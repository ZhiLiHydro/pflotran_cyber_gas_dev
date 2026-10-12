#!/usr/bin/env python3
"""End-to-end checks of Richards flow and reactive transport as one model
(MODE RICHARDS_DAE with MODE GIRT_DAE) and of its discrete adjoint.

Post-processor of adjoint-coupled and adjoint-coupled-np2, or run directly:
python3 check_coupled.py --executable <pflotran> [--np 2]
(parallel runs need PFLOTRAN_TEST_MPIEXEC, PETSc's mpiexec or one on PATH).

Backward Euler reproduces RICHARDS then GIRT; the Jacobian, including the
pressure derivatives of transport, matches finite differences entry by
entry (adjoint-coupled's unequal cells weight the harmonic dispersion
averages unevenly), also on an unstructured mixed-element grid; gradients
(permeability, mineral and microbial rate constants) match central
differences at two step sizes, initial pressures at one, also with soil
compressibility; CALCULATE_INITIAL_POROSITY reaches flow and
transport alike. -ts_model_newton_kappa 1e-4 lets the forward runs resolve
the differences.
"""

import argparse
import math
import os
from pathlib import Path
import re
import sys
import tempfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from adjoint_verify import deck, require, set_value
from adjoint_verify import (load_deck, parse_summary, parse_rerun_check,
                            run_deck, read_sensitivities, cleanup_run)

KAPPA = ("-ts_model_newton_kappa", "1e-4")
PARALLEL = ("-tran_ksp_type", "gmres", "-tran_ksp_rtol", "1e-13")
PERMEABILITY = r"(PERM_ISO\s+)(\S+)"
CALCITE_RATE = r"(Calcite\n\s+RATE_CONSTANT\s+)(\S+)"
MICROBIAL_RATE = r"(MICROBIAL_REACTION.*?RATE_CONSTANT\s+)(\S+)"


def close(a, b, rtol, atol=0.):
    require(math.isfinite(a) and math.isfinite(b) and
            abs(a-b) <= atol + rtol * max(abs(a), abs(b)),
            f"mismatch: {a:.16e} versus {b:.16e}")


def run(text, tag, np=1, args=(), kappa=True, **kwargs):
    args = (KAPPA if kappa else ()) + tuple(args)
    if np > 1:
        text = text.replace("    SOLVER DIRECT\n", "")
        args = PARALLEL + args
    return run_deck(text, tag, np=np, args=args, **kwargs)


def forward(text):
    """The deck without its ADJOINT block, with regression output of every
    cell (in place of the deck's own)."""
    text = re.sub(r"(?ms)^ADJOINT\n.*?^END\n", "", text)
    text = re.sub(r"(?ms)^REGRESSION\n.*?^END\n", "", text)
    cells = "".join(f"    {i}\n" for i in range(1, 21))
    return text.replace("\nSUBSURFACE\n", "\nSUBSURFACE\nREGRESSION\n"
                        f"  CELL_IDS\n{cells}  /\nEND\n", 1)


def regression_values(output):
    values = {}
    field = None
    for line in Path(output).with_suffix(".regression").read_text().splitlines():
        if line.startswith("--"):
            field = line.strip("- ")
            continue
        if field and "SOLUTION" not in field and ":" in line:
            key, value = line.split(":", 1)
            values[field, key.strip()] = float(value)
    return values


def sequential_scheme():
    # backward Euler on the steps of the classic scheme solves its equations
    coupled = forward(deck("coupled_ar"))
    classic = coupled.replace("MODE RICHARDS_DAE", "MODE RICHARDS").replace(
        "MODE GIRT_DAE", "MODE GIRT").replace(
        "NUMERICAL_METHODS TRANSPORT",
        "NUMERICAL_METHODS FLOW\n  TIMESTEPPER\n    TS_ACCELERATION 100\n"
        "  /\n  LINEAR_SOLVER\n    SOLVER DIRECT\n  /\nEND\n\n"
        "NUMERICAL_METHODS TRANSPORT", 1)
    a = regression_values(run(coupled, "sequential", require_summary=False,
                              args=("-tran_ts_type", "beuler")))
    b = regression_values(run(classic, "sequential", require_summary=False,
                              kappa=False))
    require(a.keys() == b.keys() and len(a) > 100, "regression fields differ")
    for key, value in a.items():
        close(value, b[key], rtol=1e-8, atol=1e-20)
    cleanup_run("sequential")
    print("PASS backward Euler reproduces RICHARDS then GIRT", flush=True)


def compressible(text, compressibility="1.d-6"):
    """Every material with soil compressibility about 1 atm."""
    new = re.sub(r"(MATERIAL_PROPERTY \w+\n  ID \d+\n  POROSITY \S+\n)",
                 r"\1  SOIL_COMPRESSIBILITY_FUNCTION LEIJNSE\n"
                 rf"  SOIL_COMPRESSIBILITY {compressibility}\n"
                 r"  SOIL_REFERENCE_PRESSURE 101325.d0\n", text)
    require(new.count("SOIL_COMPRESSIBILITY_FUNCTION") ==
            text.count("MATERIAL_PROPERTY"), "missing POROSITY")
    return new


def rock_surface(text):
    """Equilibrium surface complexation on ROCK_DENSITY sites."""
    new = text.replace(
        "  DATABASE ",
        "  SORPTION\n    SURFACE_COMPLEXATION_RXN\n      EQUILIBRIUM\n"
        "      ROCK_DENSITY\n      SITE >FeOH 1.d-3\n      COMPLEXES\n"
        "        >FeO-\n        >FeOH2+\n      /\n    /\n  /\n  DATABASE ", 1)
    require(new != text, "missing DATABASE")
    rock = re.sub(r"(  TORTUOSITY \S+\n)", r"\1  ROCK_DENSITY 2650.d0\n", new)
    require(rock.count("ROCK_DENSITY 2650") == text.count("MATERIAL_PROPERTY"),
            "missing TORTUOSITY")
    return rock


def millington_quirk(text):
    """Millington-Quirk tortuosity (the deck's tortuosities are 1)."""
    new = text.replace("      MODE GIRT_DAE\n    /\n",
                       "      MODE GIRT_DAE\n      OPTIONS\n"
                       "        USE_MILLINGTON_QUIRK_TORTUOSITY\n"
                       "      /\n    /\n", 1)
    require(new != text, "missing GIRT_DAE block")
    return new


def model_checks():
    """Model contract, entry-wise Jacobian and parameter Jacobians: microbial
    rates with and without biomass, a general reaction, dispersion, two
    materials, a mass-rate source, and soil compressibility with
    Millington-Quirk tortuosity."""
    ar = deck("coupled_ar")
    no_biomass = re.sub(r"    BIOMASS\n.*?    /\n", "", ar, count=1,
                        flags=re.S).replace("RATE_CONSTANT 4.d-2",
                                            "RATE_CONSTANT 1.d-9")
    general = ar.replace(
        "  IMMOBILE_DECAY_REACTION",
        "  GENERAL_REACTION\n    REACTION CH2O(aq) <-> CO2(aq)\n"
        "    FORWARD_RATE 1.d-6\n    BACKWARD_RATE 2.d-7\n  /\n"
        "  IMMOBILE_DECAY_REACTION", 1)
    for name, text in (("calcite", deck("adjoint-coupled")), ("biomass", ar),
                       ("no_biomass", no_biomass), ("general", general),
                       ("compressible", compressible(millington_quirk(ar))),
                       ("rock_compressible", compressible(rock_surface(
                           deck("adjoint-coupled"))))):
        out = Path(run(text, "model_checks", args=(
            "-ts_model_contract_test", "-ts_model_column_check",
            "-adjoint_debug_jacobianp")))
        log = out.read_text()
        contract = re.findall(r"TS model contract: .*: (PASS|FAIL)", log)
        require(len(contract) == 6 and set(contract) == {"PASS"},
                f"{name}: model contract {contract}")
        match = re.search(r"column check:\s*(\d+) entries off", log)
        require(match and int(match.group(1)) == 0,
                f"{name}: Jacobian column check failed or missing")
        errors = re.findall(r"ADJOINT debug dF/dp .*max rel diff =\s*(\S+)",
                            log)
        require(len(errors) == 2 and all(float(e) < 1e-8 for e in errors),
                f"{name}: parameter Jacobian check {errors}")
    cleanup_run("model_checks")
    print("PASS model contract, entry-wise Jacobian and parameter Jacobians",
          flush=True)


def unstructured_grid():
    """Model contract and entry-wise Jacobian on the unstructured
    mixed-element grid of girt_dae/mixed_implicit_dae, from a uniform
    initial pressure: at hydrostatic equilibrium every face has zero flux,
    where the upwinded fluxes are not differentiable."""
    source = HERE.parent / "girt_dae" / "mixed_implicit_dae.in"
    text = load_deck(source)
    for pattern, replacement in (
            # no snapshots
            (r"(?ms)^OUTPUT\n.*?^END\n", ""),
            (r"(?ms)^FLOW_CONDITION initial\n.*?^END\n",
             "FLOW_CONDITION initial\n  TYPE\n    LIQUID_PRESSURE DIRICHLET\n"
             "  /\n  LIQUID_PRESSURE 9.d4\nEND\n"),
            (r"FINAL_TIME [^\n]*", "FINAL_TIME 1.d-6 d")):
        text, count = re.subn(pattern, replacement, text, count=1)
        require(count == 1, f"unstructured deck: missing {pattern}")
    text = re.sub(r"(?mi)^(\s*(?:TYPE UNSTRUCTURED|FILE)\s+)(\S+)",
                  lambda m: m.group(1) + str((source.parent / m.group(2))
                                             .resolve()), text)
    log = Path(run(text, "unstructured", require_summary=False, args=(
        "-ts_model_contract_test", "-ts_model_column_check"))).read_text()
    contract = re.findall(r"TS model contract: .*: (PASS|FAIL)", log)
    require(len(contract) == 6 and set(contract) == {"PASS"},
            f"unstructured: model contract {contract}")
    match = re.search(r"column check:\s*(\d+) entries off", log)
    require(match and int(match.group(1)) == 0,
            "unstructured: Jacobian column check failed or missing")
    cleanup_run("unstructured")
    print("PASS unstructured grid: model contract and entry-wise Jacobian",
          flush=True)


def fd_check(text, parameters, tag, np=1, tolerance=1e-6):
    """dJ/dp of the cost against central differences of each parameter's
    deck value at relative steps 1e-4 and 1e-5."""
    costs, gradients = parse_summary(run(text, tag, np=np))
    for index, (pattern, value) in enumerate(parameters, start=1):
        for relative in (1e-4, 1e-5):
            h = value * relative
            plus = parse_summary(run(set_value(text, pattern, value+h), tag,
                                     np=np))[0]
            minus = parse_summary(run(set_value(text, pattern, value-h), tag,
                                      np=np))[0]
            fd = (plus[1] - minus[1]) / (2*h)
            close(fd, gradients[1, index], rtol=tolerance)
    cleanup_run(tag)
    print(f"PASS {tag}: FD agrees at both step sizes ({np} ranks)", flush=True)
    return gradients


def gradients():
    fd_check(deck("adjoint-coupled"),
             [(PERMEABILITY, 1e-12), (CALCITE_RATE, 1e-8)], "calcite")
    fd_check(deck("coupled_ar"),
             [(PERMEABILITY, 1e-12), (MICROBIAL_RATE, 4e-2)], "biomass")


def rerun(np):
    """-ts_model_rerun_check on the coupled model: reruns reproduce the run
    bit for bit at the deck's permeability and calcite rate constant, and
    at 1.1 times them give the cost and gradients of a run of the deck at
    those values, bit for bit."""
    text = deck("adjoint-coupled")
    parameters, costs, gradients = parse_rerun_check(
        run(text, "rerun", np=np, args=("-ts_model_rerun_check",)))
    scaled = set_value(set_value(text, PERMEABILITY, parameters[1]),
                       CALCITE_RATE, parameters[2])
    standalone = parse_summary(run(scaled, "rerun", np=np))
    require(standalone == (costs, gradients),
            f"the rerun at 1.1 p0 gives {costs} {gradients}, a run of the "
            f"deck {standalone}")
    cleanup_run("rerun")
    print(f"PASS coupled reruns reproduce the run and a run of the deck "
          f"({np} ranks)", flush=True)


def optimization(np):
    """ADJOINT OPTIMIZATION of the coupled model recovers the permeability
    and the calcite rate constant of a twin experiment, OBSERVATIONS of
    Ca++ at two cells at the final time, from 3x and 0.3x their values."""
    base = deck("adjoint-coupled")
    probe = ("  COST_FUNCTION TERMINAL_SPECIES_CONCENTRATION\n"
             "    SPECIES_NAME Ca++\n    REGION probe\n  /\n")
    require(probe in base, "adjoint-coupled: no probe cost")
    two = base.replace(probe, probe + probe.replace("probe", "well"))
    observed = parse_summary(run(two, "optimization", np=np))[0]
    with tempfile.TemporaryDirectory(
            prefix="pflotran-adjoint-optimization-") as directory:
        misfits = "  OBSERVATIONS\n"
        for i, z in ((1, "6.25d0"), (2, "3.25d0")):
            observations = Path(directory) / f"observations{i}.dat"
            observations.write_text("TIME_UNITS y\nDATA_UNITS mol/kg\n"
                                    f"5.d-2 {observed[i]:.17e}\n")
            misfits += ("    TIME_SERIES\n"
                        f"      COORDINATE 0.5d0 0.5d0 {z}\n"
                        "      OBSERVED_VARIABLE FREE_ION Ca++\n"
                        "      STANDARD_DEVIATION 1.d-6\n"
                        f"      FILE {observations}\n    /\n")
        twin = base.replace(probe, misfits + "  /\n")
        twin = set_value(set_value(twin, PERMEABILITY, 3e-12),
                         CALCITE_RATE, 3e-9)
        optimized = twin.replace("ADJOINT\n", "ADJOINT\n  OPTIMIZATION\n"
                                 "    MAXIMUM_ITERATIONS 60\n  /\n", 1)
        # serially, also through a failed evaluation: the one after it
        # equals a run of the deck at its parameters bit for bit
        for args in ((), ("-adj_opt_fail_eval", "3"))[:2 if np == 1 else 1]:
            output = run(optimized, "optimization", np=np, args=args)
            text = Path(output).read_text()
            require(re.search(r"OPTIMIZATION: \d+ iterations, \d+ "
                              r"evaluations \(\d+ failed\); converged", text),
                    f"TAO did not converge {args}")
            for name, value in (("PERMEABILITY\\[soil1\\]", 1e-12),
                                ("MINERAL_RATE_CONSTANT\\[Calcite\\]", 1e-8)):
                found = re.findall(rf"OPTIMIZATION: {name} = +(\S+)", text)
                require(found, f"no optimum of {name}")
                close(float(found[-1]), value, rtol=1e-8)
            if not args:
                continue
            rows = [line.split() for line in Path(output).with_name(
                "case-optimization.log").read_text().splitlines()
                if not line.startswith("#")]
            require(int(rows[2][1]) > 0, "evaluation 3 did not fail")
            f, p1, p2, g1, g2 = (float(rows[3][i]) for i in (2, 4, 5, 6, 7))
            costs, gradients = parse_summary(run(
                set_value(set_value(twin, PERMEABILITY, p1), CALCITE_RATE, p2),
                "optimization", np=np))
            require(costs[1] == f and gradients[1, 1] == g1 and
                    gradients[1, 2] == g2,
                    f"evaluation 4 after the failure: {rows[3]}, a run of the "
                    f"deck: {costs} {gradients}")
    cleanup_run("optimization")
    print(f"PASS coupled optimization recovers permeability and calcite rate "
          f"({np} ranks)", flush=True)


def compressible_gradients():
    """dJ/dk and dJ/dP0 with soil compressibility; the saturated-zone
    response is nonlinear, hence a small step and a tight Newton tolerance."""
    fd_check(compressible(deck("adjoint-coupled")), [(PERMEABILITY, 1e-12)],
             "compressible")
    # sorption on rock sites follows the solid fraction 1 - phi(P)
    fd_check(compressible(rock_surface(deck("adjoint-coupled"))),
             [(PERMEABILITY, 1e-12)], "rock_compressible")
    text = compressible(deck("coupled_ar"))
    args = ("-ts_model_newton_kappa", "1e-8")
    field = read_pressure_sensitivities(
        run(text, "compressible_p0", full_field=True, kappa=False, args=args))
    pattern = (r"(FLOW_CONDITION saturated\n  TYPE\n    LIQUID_PRESSURE "
               r"DIRICHLET\n  /\n  LIQUID_PRESSURE )(\S+)")
    value, h = 106000., 10.
    plus = parse_summary(run(set_value(text, pattern, value+h),
                             "compressible_p0", kappa=False, args=args))[0][1]
    minus = parse_summary(run(set_value(text, pattern, value-h),
                              "compressible_p0", kappa=False, args=args))[0][1]
    close((plus-minus)/(2*h), sum(field[i] for i in range(0, 10)), rtol=1e-6)
    cleanup_run("compressible_p0")
    print("PASS soil compressibility: initial-pressure sensitivity",
          flush=True)


def calculated_porosity():
    """CALCULATE_INITIAL_POROSITY (1 - calcite volume fraction 0.01) matches
    that porosity given in the deck."""
    text = forward(deck("adjoint-coupled"))
    calculated = text.replace(
        "  LOG_FORMULATION\n",
        "  LOG_FORMULATION\n  CALCULATE_INITIAL_POROSITY\n", 1)
    given = re.sub(r"  POROSITY \S+\n", "  POROSITY 0.99d0\n", text)
    require(calculated != text and given.count("POROSITY 0.99d0") == 2,
            "calculated porosity decks")
    a = regression_values(run(calculated, "porosity", require_summary=False))
    b = regression_values(run(given, "porosity", require_summary=False))
    require(a.keys() == b.keys() and len(a) > 100, "regression fields differ")
    for key, value in a.items():
        close(value, b[key], rtol=1e-12, atol=1e-30)
    cleanup_run("porosity")
    print("PASS calculated initial porosity in flow and transport",
          flush=True)


def initial_pressures():
    """dJ/dP0 summed over a uniform-pressure zone against central
    differences, saturated (water density, so the initial molality of a
    molarity constraint, varies with pressure) and unsaturated. dJ/dP0
    holds the initial molarity fixed. The probe is in the saturated zone,
    whose initial concentrations reach it."""
    text = deck("coupled_ar").replace("COORDINATE 0.5d0 0.5d0 5.25d0",
                                      "COORDINATE 0.5d0 0.5d0 2.25d0", 1)
    require("2.25d0" in text, "missing probe")
    out = run(text, "initial_pressure", full_field=True)
    field = read_pressure_sensitivities(out)
    # saturated: a compressibility effect, linear, so a large step;
    # unsaturated: O(h^2) difference error, so a small one
    for condition, value, cells, h in (
            ("saturated", 106000., range(0, 10), 1000.),
            ("unsaturated", 99000., range(10, 20), 1.)):
        pattern = (rf"(FLOW_CONDITION {condition}\n  TYPE\n    LIQUID_PRESSURE "
                   r"DIRICHLET\n  /\n  LIQUID_PRESSURE )(\S+)")
        expected = sum(field[i] for i in cells)
        require(abs(expected) > 0., f"{condition}: insensitive")
        plus = parse_summary(run(set_value(text, pattern, value+h),
                                 "initial_pressure"))[0][1]
        minus = parse_summary(run(set_value(text, pattern, value-h),
                                  "initial_pressure"))[0][1]
        close((plus-minus)/(2*h), expected, rtol=1e-6)
    cleanup_run("initial_pressure")
    print("PASS initial-pressure sensitivities (saturated and unsaturated)",
          flush=True)


def read_pressure_sensitivities(output):
    """dJ/dP0 of cost 1 by natural cell id - 1."""
    field = read_sensitivities(output, 1, suffix="dP0")
    return [field[i, 1] for i in range(1, len(field)+1)]


def rejections():
    ar = deck("coupled_ar")
    run(ar.replace("MODE GIRT_DAE", "MODE GIRT"), "rejections",
        expected_error="MODE RICHARDS_DAE runs with transport MODE GIRT_DAE")
    run(ar.replace("MODE RICHARDS_DAE", "MODE RICHARDS"), "rejections",
        expected_error="MODE GIRT_DAE runs without flow or with flow MODE "
        "RICHARDS_DAE")
    run(ar.replace("      MODE RICHARDS_DAE\n    /\n",
                   "      MODE RICHARDS_DAE\n      OPTIONS\n        STEADY_STATE\n"
                   "      /\n    /\n", 1), "rejections",
        expected_error="does not yet support steady state flow")
    run(ar.replace("\nSUBSURFACE\n", "\nSUBSURFACE\nONLY_VERTICAL_FLOW\n", 1),
        "rejections", expected_error="does not yet support ONLY_VERTICAL_FLOW")
    run(ar.replace("NUMERICAL_METHODS TRANSPORT",
                   "NUMERICAL_METHODS FLOW\nEND\n\nNUMERICAL_METHODS TRANSPORT"),
        "rejections", expected_error="Flow and transport are one model")
    # NUMERICAL_METHODS TRANSPORT configures transport, not flow
    for keyword, error in (
            ("ITOL_RELATIVE_UPDATE 1.d-6", "MODE GIRT_DAE does not support "
             "the classic GIRT solver options ITOL_ABSOLUTE_UPDATE/"
             "ITOL_RELATIVE_UPDATE"),
            ("PRESSURE_CHANGE_LIMIT 1.d5", 'Keyword "PRESSURE_CHANGE_LIMIT" '
             "not recognized in SUBSURFACE,NUMERICAL_METHODS,TRANSPORT,"
             "NEWTON_SOLVER")):
        run(ar.replace("NUMERICAL_METHODS TRANSPORT\n",
                       "NUMERICAL_METHODS TRANSPORT\n  NEWTON_SOLVER\n"
                       f"    {keyword}\n  /\n", 1), "rejections",
            expected_error=error)
    run(ar.replace("  LONGITUDINAL_DISPERSIVITY 0.1d0",
                   "  LONGITUDINAL_DISPERSIVITY 0.1d0\n"
                   "  TRANSVERSE_DISPERSIVITY_H 0.1d0"),
        "rejections", expected_error="does not support transverse dispersion")
    run(ar.replace("    MATERIAL soil1", "    MATERIAL soil9"), "rejections",
        expected_error='material "soil9" not found')
    # GIRT_DAE alone has no flow: no permeability derivative
    short = deck("ar_short").replace(
        "  PARAMETER RATE_CONSTANT",
        "  PARAMETER PERMEABILITY\n    MATERIAL soil1\n  /\n"
        "  PARAMETER RATE_CONSTANT", 1)
    run(short, "rejections",
        expected_error="PARAMETER PERMEABILITY requires flow in the equation")
    # dF/dk omits the derivative of a permeability-scaled source rate
    scaled = deck("adjoint-coupled").replace(
        "    RATE MASS_RATE\n", "    RATE SCALED_MASS_RATE PERM\n", 1)
    require(scaled != deck("adjoint-coupled"), "missing MASS_RATE")
    run(scaled, "rejections", expected_error="PARAMETER PERMEABILITY does "
        "not support source/sink well: its rate is scaled by permeability")
    run(scaled.replace("SCALED_MASS_RATE PERM", "SCALED_MASS_RATE VOLUME"),
        "rejections")
    cleanup_run("rejections")
    print("PASS rejections (transport mode, flow mode, steady flow, flow "
          "solver block, classic solver options, transverse dispersion, "
          "material, permeability without flow, permeability-scaled "
          "source)", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--np", type=int, default=1)
    parser.add_argument("--executable")
    options = parser.parse_args()
    if options.executable:
        os.environ["PFLOTRAN_TEST_EXECUTABLE"] = \
            str(Path(options.executable).resolve())
    test_name = os.environ.get("PFLOTRAN_TEST_NAME")
    np = 2 if test_name and test_name.endswith("-np2") else options.np
    if test_name:
        parse_summary(Path(test_name + ".out"))
    if np == 1:
        sequential_scheme()
        model_checks()
        unstructured_grid()
        gradients()
        initial_pressures()
        compressible_gradients()
        calculated_porosity()
        rerun(np)
        optimization(np)
        rejections()
    else:
        # the parallel gradients equal the serial ones
        reference = parse_summary(run(deck("adjoint-coupled"), "parallel"))[1]
        costs, parallel_gradients = parse_summary(
            run(deck("adjoint-coupled"), "parallel", np=np))
        for key, value in reference.items():
            close(value, parallel_gradients[key], rtol=1e-9)
        cleanup_run("parallel")
        print(f"PASS {np} ranks match one", flush=True)
        rerun(np)
        optimization(np)
    print("PASS coupled flow and transport adjoint suite", flush=True)


if __name__ == "__main__":
    main()
