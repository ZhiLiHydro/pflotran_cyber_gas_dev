#!/usr/bin/env python3
"""Behavioral checks of MODE GIRT_DAE that a gold file cannot express.

1. Model contract at the initial state, for the condensed operator and the
   assembled matrix: repeatable evaluation, history-independent
   linearization, J(a)v against finite differences (a = 0, 1/dt), exact
   preconditioner and its transpose, same run from both operators;
   Jacobian entries of activity-unit and biomass microbial rates, of
   equilibrium sorption on a mineral-hosted site density and of complexes
   of different species on one site; the deck's Newton tolerances
   converge without the model's Newton test.
2. Discrete adjoint (backward Euler): dJ/dY0 passes a Taylor test,
   including multirate and equilibrium surface complexation.
3. Forcing is independent of the waypoint intervals; a STEP velocity
   holds until its knot; snapshots inside an interval report the face
   fluxes of their time. Fixed-step backward Euler reproduces classic
   GIRT (tracer; closed calcite without depletion or activity
   coefficients).
4. A closed cell conserves calcium under BDF2; BDF1/BDF2 reach their
   orders in every reported quantity.
5. Mineral depletion is a located transition with nonnegative volume
   fractions. A BDF2 step failure redoes the interval at order 1.
6. PERIODIC TIMESTEP output, an inactive boundary cell and a negative H+
   inventory match classic GIRT; domain errors recover under FP traps;
   no PETSc leaks.
7. Unsupported physics (named in one message), classic GIRT solver
   options and coupled process models are rejected.
8. With --mpiexec, two processes match one.
"""

import argparse
import math
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile


def physical_values(path):
    values = {}
    section = None
    for line in path.read_text().splitlines():
        if line.startswith("-- "):
            kind, name = line[3:].removesuffix(" --").split(":", 1)
            section = (kind, name.strip()) if kind != "SOLUTION" else None
        elif section and ":" in line:
            label, data = line.split(":", 1)
            for component, token in enumerate(data.split()):
                token = re.sub(r'(?<=\d)([+-]\d{3})$', r'E\1', token)
                value = float(token)
                assert math.isfinite(value), (path, section, label, value)
                values[(*section, label.strip(), component)] = value
    assert values, path
    return values


def replace_once(text, old, new):
    count = text.count(old)
    assert count == 1, (old, count)
    return text.replace(old, new)


def closed_calcite(deck, final_time, activity):
    """One closed cell of calcite_dae: no flow, no boundaries, dissolving."""
    deck = replace_once(deck, "ACTIVITY_COEFFICIENTS TIMESTEP",
                        f"ACTIVITY_COEFFICIENTS {activity}")
    deck = replace_once(deck, "DATASET 1.d0 0.d0 0.d0 m/yr",
                        "DATASET 0.d0 0.d0 0.d0 m/yr")
    deck = replace_once(deck, "NXYZ 10 1 1", "NXYZ 1 1 1")
    deck = deck.replace("10.d0", "1.d0")
    deck = replace_once(deck, "    5\n    10\n", "")
    deck = replace_once(deck, "  TIMES y 0.25 0.5 0.75\n", "")
    deck = replace_once(deck, "FINAL_TIME 1.d0 y", f"FINAL_TIME {final_time} d")
    deck = replace_once(deck, "H+     1.d-8      F", "H+     5.         P")
    deck = replace_once(deck, "HCO3-  1.d-3      G  CO2(g)", "HCO3-  1.d-3      T")
    deck = replace_once(deck, "Ca++   5.d-4      M  Calcite", "Ca++   1.d-6      T")
    deck = replace_once(deck, "Calcite 1.d-5 1.d0 m^2/m^3", "Calcite 1.d-2 1.d0 m^2/m^3")
    deck = re.sub(r'BOUNDARY_CONDITION (outlet|inlet)\n.*?\nEND\n', '', deck, flags=re.S)
    return deck


def inactive_outlet(deck):
    """The last cell of a 10 m column deck inactive, under the east face."""
    deck = replace_once(deck, "STRATA\n  REGION all\n  MATERIAL soil1\nEND\n",
                        "STRATA\n  REGION all\n  MATERIAL soil1\nEND\n\n"
                        "STRATA\n  REGION last_cell\n  MATERIAL soil1\n  INACTIVE\nEND\n")
    return replace_once(deck, "\nREGION west\n",
                        "\nREGION last_cell\n  COORDINATES\n    9.d0 0.d0 0.d0\n"
                        "    10.d0 1.d0 1.d0\n  /\nEND\n\nREGION west\n")


def observation_rows(path):
    """The numeric rows of an observation file."""
    lines = path.read_text().splitlines()[1:]
    return [[float(token) for token in line.split()] for line in lines if line.strip()]


def fixed_steps(deck, dt):
    deck = re.sub(r"INITIAL_TIMESTEP_SIZE [^\n]*", f"INITIAL_TIMESTEP_SIZE {dt} d", deck)
    deck = re.sub(r"MAXIMUM_TIMESTEP_SIZE [^\n]*", f"MAXIMUM_TIMESTEP_SIZE {dt} d", deck)
    return deck


FIXED = ["-tran_ts_adapt_type", "none", "-tran_snes_atol", "1.e-30",
         "-tran_snes_rtol", "1.e-13", "-tran_snes_stol", "1.e-15",
         "-ts_model_newton_kappa", "0"]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--mpiexec")
    parser.add_argument("--remove-on-success", action="store_true")
    args = parser.parse_args()
    executable = args.executable.resolve()
    here = Path(__file__).resolve().parent
    database = (here / "../../../database/hanford.dat").resolve()
    artifacts = Path(tempfile.mkdtemp(prefix="pflotran-girt-dae-"))
    print(f"Artifacts: {artifacts}", flush=True)

    def run(name, deck, options=(), expected_error=None, np=1, timeout=300,
            iterative=False, files=()):
        case = artifacts / name
        case.mkdir()
        for f in files:
            shutil.copy(f, case / Path(f).name)
        deck = deck.replace("../../../database/hanford.dat", str(database))
        deck = deck.replace("../../ascem/batch/",
                            str(here.parent.parent / "ascem" / "batch") + "/")
        if np > 1 or iterative:
            deck = deck.replace("    SOLVER DIRECT\n", "")
            options = [*options, "-tran_ksp_type", "gmres", "-tran_pc_type",
                       "bjacobi", "-tran_ksp_rtol", "1.e-12"]
        (case / "case.in").write_text(deck)
        command = [str(executable), "-input_prefix", "case",
                   "-successful_exit_code", "0", *options]
        if np > 1:
            command = [*shlex.split(args.mpiexec), "-np", str(np), *command]
        result = subprocess.run(command, cwd=case, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=timeout)
        (case / "run.log").write_text(result.stdout)
        output = result.stdout
        if expected_error:
            assert expected_error in output, f"{name}: expected rejection missing"
            return case, output
        assert result.returncode == 0, f"{name}: exit {result.returncode}"
        assert "ERROR" not in output and "TS failed" not in output, name
        return case, output

    calcite = (here / "calcite_dae.in").read_text()
    calcite_log = (here / "calcite_log_dae.in").read_text()
    tracer = (here / "tracer_dae.in").read_text()
    microbial = (here / "microbial_dae.in").read_text()
    sorption_mr = (here / "sorption_mr_dae.in").read_text()
    sorption_vf = (here / "sorption_vf_dae.in").read_text()
    sorption_den = (here / "sorption_den_dae.in").read_text()
    srfcplx = (here / "srfcplx_negative_dae.in").read_text()
    # sorption_vf with its SITE_FRACTION/RATES skipped: equilibrium sorption
    # whose site density scales with the calcite inventory
    sorption_vf_eq = replace_once(replace_once(sorption_vf, "\n#skip\n", "\nskip\n"),
                                  "\n#noskip\n", "\nnoskip\n")
    # and a sorbed parent that decays: its rate depends on the calcite
    # inventory through the site density
    decay_vf_eq = replace_once(sorption_vf_eq, "  SORPTION\n",
                               "  RADIOACTIVE_DECAY_REACTION\n    REACTION UO2++ <-> Ca++\n"
                               "    RATE_CONSTANT 1.d-8\n  /\n  SORPTION\n")

    # 1. model contract, for the condensed operator (default; its exact
    #    preconditioner is an extra check) and for the assembled matrix
    for name, deck in (("contract_calcite", calcite), ("contract_calcite_log", calcite_log),
                       ("contract_microbial", microbial), ("contract_tracer", tracer),
                       ("contract_sorption_mr", sorption_mr),
                       ("contract_sorption_vf", sorption_vf),
                       ("contract_sorption_vf_eq", sorption_vf_eq)):
        deck = re.sub(r"FINAL_TIME [^\n]*", "FINAL_TIME 1.d-6 y", deck)
        for suffix, options, nchecks in (("", [], 6),
                                         ("_assembled", ["-rt_dae_assembled_jacobian"], 5)):
            _, out = run(name + suffix, deck, ["-ts_model_contract_test", *options])
            lines = [l for l in out.splitlines() if "TS model contract" in l]
            assert len(lines) == nchecks and \
                all(l.rstrip().endswith("PASS") for l in lines), (name + suffix, lines)
    # the two operators solve the same equations: a run agrees to the
    # Newton tolerance
    condensed, _ = run("operator_condensed", calcite, [])
    assembled, _ = run("operator_assembled", calcite, ["-rt_dae_assembled_jacobian"])
    a = physical_values(assembled / "case.regression")
    b = physical_values(condensed / "case.regression")
    err = max(abs(b[key] - value) / max(abs(value), 1.e-14) for key, value in a.items()
              if key[0] in ("CONCENTRATION", "VOLUME_FRACTION", "GENERIC"))
    assert err < 1.e-6, ("condensed vs assembled operator", err)
    print("PASS: model contract (repeatability, cold linearization, J(a)v at a = 0 and 1/dt, "
          "exact condensed preconditioner, transposed preconditioner); condensed and assembled "
          f"operators agree ({err:.1e})", flush=True)

    # entry-wise Jacobian (a norm-wise check hides wrong small entries):
    # activity-coefficient chain of an activity-unit microbial rate on charged
    # species, and the biomass column of a microbial rate in a 10 m^3 cell
    microbial_activity = replace_once(
        calcite, "  ACTIVITY_COEFFICIENTS TIMESTEP\n",
        "  MICROBIAL_REACTION\n    REACTION Ca++ <-> H+\n    RATE_CONSTANT 1.d-9\n"
        "    MONOD\n      SPECIES_NAME Ca++\n      HALF_SATURATION_CONSTANT 1.d-3\n"
        "    /\n    CONCENTRATION_UNITS ACTIVITY\n  /\n  ACTIVITY_COEFFICIENTS TIMESTEP\n")
    biomass = replace_once(microbial, "    B(aq)\n  /\n  MICROBIAL_REACTION",
                           "    B(aq)\n  /\n  IMMOBILE_SPECIES\n    D(im)\n  /\n"
                           "  MICROBIAL_REACTION")
    biomass = replace_once(biomass, "    CONCENTRATION_UNITS ACTIVITY\n  /",
                           "    CONCENTRATION_UNITS ACTIVITY\n    BIOMASS\n"
                           "      SPECIES_NAME D(im)\n      YIELD 0.01d0\n    /\n  /")
    biomass = replace_once(biomass, "    B(aq) 1.d-7 T\n  /",
                           "    B(aq) 1.d-7 T\n  /\n  IMMOBILE\n    D(im) 1.d-2\n  /")
    # grid bounds and region "all" together
    assert biomass.count("    1.d0 1.d0 1.d0\n  /") == 2
    biomass = biomass.replace("    1.d0 1.d0 1.d0\n  /", "    10.d0 1.d0 1.d0\n  /")
    biomass = replace_once(biomass, "RATE_CONSTANT 1.d-8", "RATE_CONSTANT 1.d-4")
    for name, deck in (("columns_calcite", calcite),
                       ("columns_microbial_activity", microbial_activity),
                       ("columns_biomass", biomass),
                       ("columns_sorption_vf_eq", sorption_vf_eq),
                       ("columns_decay_sorbed", decay_vf_eq),
                       ("columns_sorption_mr", sorption_mr),
                       ("columns_site_competition", srfcplx)):
        deck = re.sub(r"FINAL_TIME [^\n]*", "FINAL_TIME 1.d-6 y", deck)
        _, out = run(name, deck, ["-ts_model_column_check"])
        match = re.search(r"column check:\s*(\d+) entries off", out)
        assert match and int(match.group(1)) == 0, (name, match and match.group(0))
    # with the model's Newton test disabled the deck's tolerances converge
    run("kappa0", re.sub(r"FINAL_TIME [^\n]*", "FINAL_TIME 1.d-2 y", calcite),
        ["-ts_model_newton_kappa", "0"])
    print("PASS: Jacobian entry by entry (activity-unit and biomass microbial rates, "
          "equilibrium sorption on a mineral site, with a decaying sorbed parent, "
          "multirate sorption, complexes competing for a site); "
          "deck Newton tolerances apply without the model test", flush=True)

    # 2. discrete adjoint on backward Euler
    # (sorption decks have 15 components; PETSc 3.24's transposed block size
    # 15 ILU solve is the forward one: the condensed operator must not use it)
    for name, deck in (("adjoint_calcite", calcite), ("adjoint_calcite_log", calcite_log),
                       ("adjoint_microbial", microbial), ("adjoint_sorption_mr", sorption_mr),
                       ("adjoint_sorption_vf", sorption_vf),
                       ("adjoint_sorption_vf_eq", sorption_vf_eq),
                       ("adjoint_sorption_den", sorption_den)):
        deck = re.sub(r"FINAL_TIME [^\n]*", "FINAL_TIME 1.d-6 y", deck)
        _, out = run(name, deck, ["-ts_model_adjoint_test"])
        # eps = 1e-1 ... 1e-6: ratios near 100 until the remainder hits
        # roundoff; a wrong gradient gives ratios near 10
        ratios = [float(m) for m in re.findall(r"ratio =\s*([0-9.]+)", out)]
        assert sum(80. < r < 125. for r in ratios) >= 2, (name, ratios)
    print("PASS: discrete adjoint Taylor test (second-order remainder)", flush=True)

    # 3. backward Euler equals classic GIRT where the formulations coincide
    fixed_tracer = fixed_steps(tracer, "5.d-2")
    girt, _ = run("tracer_girt", replace_once(fixed_tracer, "      MODE GIRT_DAE\n", "      MODE GIRT\n"))
    dae, _ = run("tracer_dae_be", fixed_tracer, ["-tran_ts_type", "beuler", *FIXED])
    a = physical_values(girt / "case.regression")
    b = physical_values(dae / "case.regression")
    peak = max(abs(v) for k, v in a.items() if k[0] == "CONCENTRATION")
    for key, value in a.items():
        if key[0] == "CONCENTRATION":
            assert abs(b[key] - value) <= 1.e-12 * peak, (key, value, b[key])
    closed = fixed_steps(closed_calcite(calcite, "0.2d0", "OFF"), "1.d-2")
    girt, _ = run("closed_girt", replace_once(closed, "      MODE GIRT_DAE\n", "      MODE GIRT\n"))
    dae, _ = run("closed_dae_be", closed, ["-tran_ts_type", "beuler", *FIXED])
    a = physical_values(girt / "case.regression")
    b = physical_values(dae / "case.regression")
    for key, value in a.items():
        if key[0] in ("CONCENTRATION", "VOLUME_FRACTION", "GENERIC"):
            assert math.isclose(b[key], value, rel_tol=1.e-8, abs_tol=1.e-20), \
                (key, value, b[key])
    print("PASS: backward Euler reproduces classic GIRT (tracer; closed calcite, "
          "activity coefficients off)", flush=True)

    # 4. conservation and orders of accuracy in a closed cell
    molar_volume = 36.934e-6
    initial_ca = 0.25 * 1000. * 1.e-6 + 1.e-2 / molar_volume
    errors = {}
    finals = {}
    for order in (1, 2):
        for dt in (1.e-2, 5.e-3, 2.5e-3, 3.125e-4):
            case, _ = run(f"closed_bdf{order}_{dt}",
                          fixed_steps(closed_calcite(calcite, "0.1d0", "NEWTON_ITERATION"),
                                      f"{dt}"),
                          ["-tran_ts_type", "bdf", "-tran_ts_bdf_order", str(order), *FIXED])
            values = physical_values(case / "case.regression")
            ca = values["CONCENTRATION", "Total Ca++", "1", 0] * 250. + \
                values["VOLUME_FRACTION", "Calcite VF", "1", 0] / molar_volume
            assert abs(ca - initial_ca) <= 1.e-10 * initial_ca, (order, dt, ca, initial_ca)
            finals[order, dt] = values
    for order in (1, 2):
        reference = finals[2, 3.125e-4]
        for key in (("CONCENTRATION", "Total Ca++", "1", 0),
                    ("VOLUME_FRACTION", "Calcite VF", "1", 0),
                    ("GENERIC", "pH", "1", 0)):
            e = [abs(finals[order, dt][key] - reference[key]) for dt in (1.e-2, 5.e-3, 2.5e-3)]
            rates = [e[i] / e[i + 1] for i in range(2)]
            expected = 2 ** order
            assert all(0.75 * expected < r < 1.3 * expected for r in rates), \
                (order, key, e, rates)
    print("PASS: closed-cell calcium conserved to 1e-10; BDF1/BDF2 converge at "
          "orders 1/2 in concentration, volume fraction and pH", flush=True)

    # forcing is independent of the segmentation: a transient uniform
    # velocity with transverse dispersion, fixed backward Euler steps, one
    # waypoint interval or ten
    moving = replace_once(tracer, "  DATASET 1.d0 0.d0 0.d0 m/d\n",
                          "  DATASET LIST\n    TIME_UNITS d\n    DATA_UNITS m/d m/d m/d\n"
                          "    INTERPOLATION LINEAR\n    0.d0 1.d0 0.1d0 0.d0\n"
                          "    1.d0 3.d0 0.5d0 0.d0\n  /\n")
    moving = replace_once(moving, "  NXYZ 10 1 1\n  BOUNDS\n    0.d0 0.d0 0.d0\n"
                          "    10.d0 1.d0 1.d0", "  NXYZ 10 5 1\n  BOUNDS\n"
                          "    0.d0 0.d0 0.d0\n    10.d0 5.d0 1.d0")
    moving = replace_once(moving, "  TORTUOSITY 1.d0\nEND",
                          "  TORTUOSITY 1.d0\n  LONGITUDINAL_DISPERSIVITY 1.d0\n"
                          "  TRANSVERSE_DISPERSIVITY_H 0.5d0\n"
                          "  TRANSVERSE_DISPERSIVITY_V 0.5d0\nEND")
    moving = replace_once(moving, "    10.d0 1.d0 1.d0\n  /\nEND\n\nREGION west",
                          "    10.d0 5.d0 1.d0\n  /\nEND\n\nREGION west")
    moving = replace_once(moving, "    0.d0 0.d0 0.d0\n    0.d0 1.d0 1.d0",
                          "    0.d0 0.d0 0.d0\n    0.d0 2.d0 1.d0")
    moving = replace_once(moving, "    10.d0 0.d0 0.d0\n    10.d0 1.d0 1.d0",
                          "    10.d0 0.d0 0.d0\n    10.d0 5.d0 1.d0")
    moving = fixed_steps(moving, "1.d-2")
    segmented = replace_once(moving, "\nTIME\n",
                             "\nOUTPUT\n  TIMES d 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9\n"
                             "END\n\nTIME\n")
    one, _ = run("velocity_one_interval", moving, ["-tran_ts_type", "beuler", *FIXED])
    ten, _ = run("velocity_ten_intervals", segmented, ["-tran_ts_type", "beuler", *FIXED])
    a = physical_values(one / "case.regression")
    b = physical_values(ten / "case.regression")
    peak = max(abs(v) for k, v in a.items() if k[0] == "CONCENTRATION")
    for key, value in a.items():
        if key[0] == "CONCENTRATION":
            assert abs(b[key] - value) <= 1.e-12 * peak, (key, value, b[key])
    print("PASS: a transient velocity with transverse dispersion does not depend on "
          "the waypoint intervals", flush=True)

    # a STEP dataset holds each value until its knot, a waypoint; classic
    # GIRT reads the velocity at the step start, so it is the reference
    stepped = replace_once(tracer, "  DATASET 1.d0 0.d0 0.d0 m/d\n",
                           "  DATASET LIST\n    TIME_UNITS d\n    DATA_UNITS m/d m/d m/d\n"
                           "    INTERPOLATION STEP\n    0.d0 1.d0 0.d0 0.d0\n"
                           "    0.5d0 3.d0 0.d0 0.d0\n  /\n")
    stepped = fixed_steps(stepped, "5.d-2")
    girt, _ = run("step_velocity_girt",
                  replace_once(stepped, "      MODE GIRT_DAE\n", "      MODE GIRT\n"))
    dae, _ = run("step_velocity_dae", stepped, ["-tran_ts_type", "beuler", *FIXED])
    a = physical_values(girt / "case.regression")
    b = physical_values(dae / "case.regression")
    peak = max(abs(v) for k, v in a.items() if k[0] == "CONCENTRATION")
    for key, value in a.items():
        if key[0] == "CONCENTRATION":
            assert abs(b[key] - value) <= 1.e-12 * peak, (key, value, b[key])
    print("PASS: a STEP velocity dataset holds each value until its knot "
          "(backward Euler equals classic GIRT)", flush=True)

    # PERIODIC TIMESTEP snapshots inside a waypoint interval report the face
    # fluxes of their time, under a LINEAR dataset too
    linear = replace_once(tracer, "  DATASET 1.d0 0.d0 0.d0 m/d\n",
                          "  DATASET LIST\n    TIME_UNITS d\n    DATA_UNITS m/d m/d m/d\n"
                          "    INTERPOLATION LINEAR\n    0.d0 1.d0 0.d0 0.d0\n"
                          "    1.d0 3.d0 0.d0 0.d0\n  /\n")
    linear = replace_once(linear, "\nSTRATA\n",
                          "\nOUTPUT\n  TIME_UNITS d\n  VELOCITY_AT_CENTER\n"
                          "  SNAPSHOT_FILE\n    PERIODIC TIMESTEP 4\n"
                          "    FORMAT TECPLOT POINT\n  /\nEND\n\nSTRATA\n")
    case, _ = run("linear_velocity_snapshots", fixed_steps(linear, "5.d-2"),
                  ["-tran_ts_type", "beuler", *FIXED])
    snapshots = sorted(case.glob("case-vel-*.tec"))
    assert len(snapshots) > 3, snapshots
    for path in snapshots:
        lines = path.read_text().splitlines()
        time = float(re.search(r'"\s*(\S+)', lines[0]).group(1))
        flux = float(lines[3].split()[3])
        assert math.isclose(flux, 1.0 + 2.0*time, rel_tol=1.e-5), (path.name, time, flux)
    print(f"PASS: {len(snapshots)} snapshots inside a waypoint interval report the "
          "face fluxes of their time (LINEAR velocity dataset)", flush=True)

    # 5. mineral depletion as a regime transition
    case, out = run("depletion", calcite, ["-tran_ts_atol", "1.e-10", "-tran_ts_rtol", "1.e-6"])
    match = re.search(r"regime transitions:\s*(\d+) depletions", out)
    assert match and int(match.group(1)) > 0, "no depletion reported"
    values = physical_values(case / "case.regression")
    assert all(v >= 0. for k, v in values.items() if k[0] == "VOLUME_FRACTION"), \
        "negative volume fraction"
    print(f"PASS: {match.group(1)} depletion events, volume fractions nonnegative",
          flush=True)

    # a step failure (no rejection allowed) at BDF order 2 redoes the rest
    # of the interval at order 1, with the interval's forcing segment
    _, out = run("order_fallback", microbial, ["-tran_ts_max_reject", "0"])
    match = re.search(r"BDF order fallbacks\s+(\d+)", out)
    assert match and int(match.group(1)) > 0, "no BDF order fallback"
    print(f"PASS: {match.group(1)} BDF order fallback(s) redo the rest of an "
          "interval", flush=True)

    # 6. PERIODIC TIMESTEP observation rows and snapshots match classic GIRT
    #    one to one (the last step, a waypoint too, is written once)
    periodic = replace_once(fixed_tracer, "\nSTRATA\n",
                            "\nOBSERVATION\n  REGION obs\nEND\n\nOUTPUT\n"
                            "  OBSERVATION_FILE\n    PERIODIC TIMESTEP 1\n  /\n"
                            "  SNAPSHOT_FILE\n    PERIODIC TIMESTEP 5\n"
                            "    FORMAT TECPLOT POINT\n  /\nEND\n\n"
                            "REGION obs\n  COORDINATE 5.5d0 0.5d0 0.5d0\nEND\n\nSTRATA\n")
    girt, _ = run("periodic_girt",
                  replace_once(periodic, "      MODE GIRT_DAE\n", "      MODE GIRT\n"))
    dae, _ = run("periodic_dae", periodic, ["-tran_ts_type", "beuler", *FIXED])
    a = observation_rows(girt / "case-obs-0.pft")
    b = observation_rows(dae / "case-obs-0.pft")
    assert len(a) == 21 and len(b) == len(a), ("observation rows", len(a), len(b))
    peak = max(abs(v) for row in a for v in row[1:])
    for ra, rb in zip(a, b):
        assert math.isclose(ra[0], rb[0], rel_tol=1.e-12), (ra, rb)
        assert all(abs(x - y) <= 1.e-12 * peak for x, y in zip(ra[1:], rb[1:])), (ra, rb)
    snapshots = [sorted(p.name for p in case.glob("case-*.tec")) for case in (girt, dae)]
    assert len(snapshots[0]) == 5 and snapshots[0] == snapshots[1], snapshots
    print("PASS: PERIODIC TIMESTEP observations and snapshots match classic GIRT "
          "one to one", flush=True)

    # an inactive cell under a boundary face has no boundary chemistry (its
    # zero molalities would make the carbonate complexes infinite), and
    # backward Euler still reproduces classic GIRT around it
    run("inactive_boundary_cell", inactive_outlet(calcite))
    girt, _ = run("inactive_girt", replace_once(inactive_outlet(fixed_tracer),
                                                "      MODE GIRT_DAE\n", "      MODE GIRT\n"))
    dae, _ = run("inactive_dae_be", inactive_outlet(fixed_tracer),
                 ["-tran_ts_type", "beuler", *FIXED])
    a = physical_values(girt / "case.regression")
    b = physical_values(dae / "case.regression")
    peak = max(abs(v) for k, v in a.items() if k[0] == "CONCENTRATION")
    for key, value in a.items():
        if key[0] == "CONCENTRATION":
            assert abs(b[key] - value) <= 1.e-12 * peak, (key, value, b[key])

    # a negative inventory of a component with no negative coefficient in
    # an aqueous complex (H+, deprotonated surface sites) is admissible
    girt, _ = run("negative_inventory_girt",
                  replace_once(srfcplx, "      MODE GIRT_DAE\n", "      MODE GIRT\n"))
    dae, _ = run("negative_inventory_dae", srfcplx)
    a = physical_values(girt / "case.regression")
    b = physical_values(dae / "case.regression")
    # sorbed [mol/m^3 bulk] + porosity * 1000 L/m^3 * aqueous [mol/L]
    inventory = b["CONCENTRATION", "Total Sorbed H+", "1", 0] + \
        0.9 * 1000. * b["CONCENTRATION", "Total H+", "1", 0]
    assert inventory < 0., ("H+ inventory not negative", inventory)
    err = max(abs(b[key] - value) / max(abs(value), 1.e-14) for key, value in a.items()
              if key[0] in ("CONCENTRATION", "GENERIC"))
    assert err < 1.e-4, ("negative inventory: GIRT_DAE vs classic GIRT", err)
    print(f"PASS: an inactive cell under a boundary face; a negative H+ inventory "
          f"({inventory:.2e} mol/m^3) matches classic GIRT ({err:.1e})", flush=True)

    # failed stages are recoverable with floating-point traps on (this deck
    # has some: A vanishes, ln(c) at its floor), also a failed closure, a
    # domain error (forced: no deck reaches one); and no PETSc object or
    # allocation leaks. A linear deck does not run under -fp_trap: a vanished
    # species' log(0) in RTotalAqueous divides by zero (harmlessly)
    microbial_log = (here / "microbial_log_dae.in").read_text()
    _, out = run("fp_trap", microbial_log, ["-fp_trap", "-rt_dae_fail_closure", "50"])
    match = re.search(r"TS steps = .* cuts = +(\d+)", out)
    assert match and int(match.group(1)) > 0, ("fp_trap: no step cuts", match)
    errors = re.search(r"domain errors +(\d+)", out)
    assert errors and errors.group(1) == "1", ("fp_trap: no domain error", errors)
    _, out = run("malloc_dump", tracer, ["-malloc_debug", "-malloc_dump"])
    assert not re.search(r"^\[\s*\d+\]\s*\d+ bytes", out, re.M), "leaked memory"
    print(f"PASS: {match.group(1)} step cuts and a domain error recovered under "
          "-fp_trap; no leaked memory", flush=True)

    # 7. explicit rejection of physics the model does not define, and of
    #    classic GIRT solver options that PETSc TS replaces
    rejected = replace_once(calcite, "  ACTIVITY_COEFFICIENTS TIMESTEP\n",
                            "  ACTIVITY_COEFFICIENTS TIMESTEP\n  UPDATE_POROSITY\n")
    run("reject_porosity", rejected, expected_error="MODE GIRT_DAE does not yet support")
    # soil compressibility too, without RICHARDS_DAE flow: both named
    compressible = re.sub(r"(MATERIAL_PROPERTY \w+\n  ID \d+\n)",
                          r"\1  SOIL_COMPRESSIBILITY_FUNCTION LEIJNSE\n"
                          r"  SOIL_COMPRESSIBILITY 1.d-7\n"
                          r"  SOIL_REFERENCE_PRESSURE 101325.d0\n", rejected)
    assert compressible != rejected
    run("reject_two", compressible, expected_error="MODE GIRT_DAE does not yet support "
        "UPDATE_POROSITY/TORTUOSITY/PERMEABILITY, transient porosity (soil "
        "compressibility, evolving STRATA) without MODE RICHARDS_DAE.")
    classic_option = replace_once(calcite, "  NEWTON_SOLVER\n",
                                  "  NEWTON_SOLVER\n    NUMERICAL_JACOBIAN\n")
    run("reject_numerical_jacobian", classic_option,
        expected_error="does not support the classic GIRT solver options")
    # the Newton update limits of classic GIRT
    run("reject_max_dlnc", replace_once(calcite, "CHEMISTRY\n", "CHEMISTRY\n  MAX_DLNC 2.d0\n"),
        expected_error="classic GIRT solver options MAX_DLNC")
    # TIMESTEPPER keywords of fields the TS driver does not use
    accelerated = replace_once(calcite, "NUMERICAL_METHODS TRANSPORT\n",
                               "NUMERICAL_METHODS TRANSPORT\n  TIMESTEPPER\n"
                               "    TS_ACCELERATION 8\n  /\n")
    run("reject_ts_acceleration", accelerated,
        expected_error="does not support the TIMESTEPPER keywords TS_ACCELERATION")
    steady = replace_once(calcite, "      MODE GIRT_DAE\n",
                          "      MODE GIRT_DAE\n      OPTIONS\n        STEADY_STATE\n"
                          "      /\n")
    run("reject_steady_state", steady,
        expected_error="MODE GIRT_DAE does not yet support steady state transport.")
    ufd = here.parent.parent / "ufd"
    glass = replace_once((ufd / "glass.in").read_text(), "      MODE GIRT\n",
                         "      MODE GIRT_DAE\n")
    # GIRT_DAE runs with flow only as RICHARDS_DAE, in one model
    glass = replace_once(glass, "      MODE RICHARDS\n", "      MODE RICHARDS_DAE\n")
    glass, count = re.subn(r"(?ms)^NUMERICAL_METHODS FLOW\n.*?^END\n", "", glass,
                           count=1)
    assert count == 1, "glass.in: no NUMERICAL_METHODS FLOW block"
    run("reject_waste_form", glass, files=[ufd / "ufd-decay.dat"],
        expected_error="does not run process models coupled to it")
    print("PASS: unsupported configuration rejected (porosity update, with soil "
          "compressibility too; NUMERICAL_JACOBIAN; TS_ACCELERATION; STEADY_STATE; waste "
          "form and UFD decay coupled to transport)", flush=True)

    # 8. parallel
    if args.mpiexec:
        short = re.sub(r"FINAL_TIME [^\n]*", "FINAL_TIME 5.d-2 y", calcite)
        options = ["-tran_ts_atol", "1.e-10", "-tran_ts_rtol", "1.e-6"]
        serial, _ = run("np1", short, options, np=1, iterative=True)
        parallel, _ = run("np2", short, options, np=2, timeout=1800)
        a = physical_values(serial / "case.regression")
        b = physical_values(parallel / "case.regression")
        for key, value in a.items():
            if key[0] in ("CONCENTRATION", "VOLUME_FRACTION", "GENERIC"):
                assert math.isclose(b[key], value, rel_tol=1.e-5, abs_tol=1.e-14), \
                    (key, value, b[key])
        print("PASS: two processes match one", flush=True)
    else:
        print("SKIP: parallel check (no --mpiexec)", flush=True)

    if args.remove_on_success:
        shutil.rmtree(artifacts)
    print("All GIRT_DAE checks passed.", flush=True)


if __name__ == "__main__":
    try:
        main()
    except AssertionError as error:
        print(f"Failed : {error}", flush=True)
        sys.exit(1)
