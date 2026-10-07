#!/usr/bin/env python3
"""End-to-end checks of the MODE GIRT_DAE discrete adjoint.

Post-processor of adjoint and adjoint-np2, or run directly:
python3 check_adjoint.py --executable <pflotran> [--np 2]
(parallel runs need PFLOTRAN_TEST_MPIEXEC, PETSc's mpiexec or one on PATH).

Gradients are compared with central differences at two step sizes on fixed
backward Euler steps (-tran_ts_adapt_type none).
"""

import argparse
import math
import os
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from adjoint_verify import deck, require, set_value
from adjoint_verify import (parse_summary, parse_rerun_check,
                            run_deck, read_sensitivities, cleanup_run)

# parallel runs: an iterative solve accurate enough for finite differences
PARALLEL = ("-tran_ksp_type", "gmres", "-tran_pc_type", "bjacobi",
            "-tran_schur_sub_pc_type", "lu", "-tran_ksp_rtol", "1e-13")
MICROBIAL_RATE = r"(MICROBIAL_REACTION.*?RATE_CONSTANT\s+)\S+"
MINERAL_RATE = r"(MINERAL_KINETICS.*?RATE_CONSTANT\s+)\S+"
# J and dJ/dp of the registered deck from an independent TS formulation of
# GIRT (conserved masses as unknowns), same backward Euler steps
REFERENCE_COST = 6.466384244876979e-5
REFERENCE_GRADIENT = 1.67281389710257e-3


def close(a, b, rtol=1e-8, atol=1e-13):
    require(math.isfinite(a) and math.isfinite(b) and
            abs(a-b) <= atol + rtol * max(abs(a), abs(b)),
            f"mismatch: {a:.16e} versus {b:.16e}")


def run(text, tag, np=1, **kwargs):
    args = tuple(kwargs.pop("args", ()))
    if np > 1:
        text = text.replace("    SOLVER DIRECT\n", "")
        args = PARALLEL + args
    return run_deck(text, tag, np=np, args=args, **kwargs)


def with_parameter(text, block):
    text, count = re.subn(r"  PARAMETER RATE_CONSTANT\n.*?  /\n",
                          lambda _: block, text, count=1, flags=re.S)
    require(count == 1, "missing PARAMETER block")
    return text


MINERAL_PARAMETER = ("  PARAMETER RATE_CONSTANT\n    KINETICS MINERAL\n"
                     "    MINERAL_NAME Gypsum\n  /\n")
RATE_PARAMETER = ("  PARAMETER RATE_CONSTANT\n    KINETICS MICROBIAL\n"
                  "    REACTION_INDEX 1\n  /\n")


def fd_check(text, value, tag, np=1, pattern=MICROBIAL_RATE, scale=1.,
             tolerance=1e-6, args=()):
    """dJ/dp of every cost against central differences of the deck value
    at relative steps 1e-4 and 1e-5; scale converts the deck's units of p
    to those of the reported gradient."""
    costs, gradients = parse_summary(run(text, tag, np=np, args=args))
    for relative in (1e-4, 1e-5):
        h = abs(value) * relative
        plus = parse_summary(run(set_value(text, pattern, value+h), tag,
                                 np=np, args=args))[0]
        minus = parse_summary(run(set_value(text, pattern, value-h), tag,
                                  np=np, args=args))[0]
        for cost in costs:
            fd = (plus[cost] - minus[cost]) / (2*h) / scale
            close(fd, gradients[cost, 1], rtol=tolerance, atol=1e-13)
    print(f"PASS {tag}: FD agrees at both step sizes ({np} ranks)", flush=True)
    cleanup_run(tag)
    return costs, gradients


def registered(test_name):
    # the deck of the registered test, against the independent reference
    costs, gradients = parse_summary(Path(test_name + ".out"))
    close(costs[1], REFERENCE_COST, rtol=1e-7)
    close(gradients[1, 1], REFERENCE_GRADIENT, rtol=1e-7)
    print("PASS registered deck matches the independent reference", flush=True)


def gradients(np):
    fd_check(deck("ar_1d4"), 0.04, "one_dimensional", np=np)
    if np > 1:
        return
    short = fd_check(deck("ar_short"), 0.04, "short")
    for name in ("ar_multicost", "ar_multiwaypoint", "ar_1d4_flux",
                 "dynamic_velocity"):
        fd_check(deck(name), 0.04, name, tolerance=1e-5)
    # tighter Newton: the default test's error enters the differences;
    # LOG_FORMULATION must give the gradient of the linear state
    log = fd_check(deck("ar_log"), 0.04, "log",
                   args=("-ts_model_newton_kappa", "1e-4"))
    close(log[1][1, 1], short[1][1, 1], rtol=1e-10)
    # a cost on an immobile species, after the primaries in the state
    for name in ("ar_short", "ar_log"):
        text = deck(name)
        immobile = text.replace("SPECIES_NAME CO2(aq)", "SPECIES_NAME Biomass")
        require(immobile != text, name + ": missing CO2(aq) cost")
        fd_check(immobile, 0.04, name + "_immobile",
                 args=("-ts_model_newton_kappa", "1e-4"))
    charged = deck("charged")
    require("ACTIVITY_COEFFICIENTS OFF" in charged,
            "charged: missing ACTIVITY_COEFFICIENTS OFF")
    for mode in ("OFF", "TIMESTEP", "NEWTON_ITERATION"):
        fd_check(charged.replace("ACTIVITY_COEFFICIENTS OFF",
                                 "ACTIVITY_COEFFICIENTS " + mode),
                 2.0, "activity_" + mode)
    mineral = deck("mineral")
    varying = mineral.replace(
        "RATE_CONSTANT 1.d-6 mol/m^2-sec",
        "RATE_CONSTANT 1.d-6 mol/m^2-sec\n"
        "      SURFACE_AREA_FUNCTION VOLUME_FRACTION_RATIO\n"
        "      SURFACE_AREA_VOL_FRAC_POWER 0.6666666667")
    require(varying != mineral, "missing mineral rate constant")
    for area, text in (("constant", mineral), ("volume_fraction", varying)):
        fd_check(text, 0.5, "mineral_" + area)
        fd_check(with_parameter(text, MINERAL_PARAMETER), 1e-6,
                 "mineral_rate_" + area, pattern=MINERAL_RATE)


def failed_stage():
    # a failed closure (forced: Newton's step guard keeps every deck's
    # iterates in the domain) is a domain error; the step is redone shorter,
    # and the adjoint differentiates the steps taken
    for name, args in (("ar_short", ()),
                       ("ar_log", ("-ts_model_newton_kappa", "1e-4"))):
        args = (*args, "-rt_dae_fail_closure", "8")
        text = Path(run(deck(name), "failed_stage", args=args)).read_text()
        match = re.search(r"domain errors +(\d+)", text)
        require(match and match.group(1) == "1", f"{name}: no domain error")
        fd_check(deck(name), 0.04, "failed_stage_" + name, args=args)


def selector_and_output(np):
    text = deck("ar_1d4")
    baseline = run(text, "selectors", np=np, full_field=True)
    costs, gradients = parse_summary(baseline)
    field = read_sensitivities(baseline, 4)
    require(len(field) == 16, "wrong full-field vector size")
    # identical initial data in the four cells; only cell 2 is selected,
    # and cell 4 is owned away from rank zero
    for selector in ("COORDINATE 1.5 0.5 0.5", "CELL_IDS 2"):
        selected = text.replace("    REGION cell2", "    " + selector)
        value = run(selected, "selectors", np=np, full_field=True)
        require(parse_summary(value) == (costs, gradients),
                f"{selector} changed the objective or gradient")
        require(read_sensitivities(value, 4) == field,
                f"{selector} changed initial-condition sensitivities")
    moved = text.replace("    REGION cell2", "    COORDINATE 3.5 0.5 0.5")
    moved_out = run(moved, "selectors", np=np, full_field=True)
    moved_field = read_sensitivities(moved_out, 4)
    close(parse_summary(moved_out)[0][1], costs[1])
    cell4 = run(text.replace("    REGION cell2", "    CELL_IDS 4"),
                "selectors", np=np, full_field=True)
    require(read_sensitivities(cell4, 4) == moved_field,
            "coordinate on the last cell disagrees with CELL_IDS")
    default = run(text, "selectors", np=np)
    require(not list(Path(default).parent.glob("*-du0.bin")),
            "default output unexpectedly wrote a full field")
    if np > 1:
        serial = run(text, "selectors", full_field=True)
        serial_field = read_sensitivities(serial, 4)
        for key in field:
            close(field[key], serial_field[key])
        serial_summary = parse_summary(serial)
        for key in gradients:
            close(gradients[key], serial_summary[1][key])
    # the parameter Jacobian check must not perturb the gradient
    debug = run(text, "selectors", np=np,
                args=("-adjoint_debug_jacobianp",), full_field=True)
    require(parse_summary(debug) == (costs, gradients),
            "parameter debug check changed summary values")
    require(read_sensitivities(debug, 4) == field,
            "parameter debug check changed initial-condition sensitivities")
    cleanup_run("selectors")
    print(f"PASS selectors, compact/full output and debug restoration "
          f"({np} ranks)", flush=True)


def final_state(np):
    # the backward sweep leaves the final forward state for the output
    text = deck("ar_1d4").replace(
        "\nSUBSURFACE\n",
        "\nSUBSURFACE\nREGRESSION\n  CELL_IDS\n    1\n    4\n  /\nEND\n", 1)
    with_adjoint = run(text, "final_state", np=np)
    forward_deck = re.sub(r"(?ms)^ADJOINT\n.*?^END\n", "", text)
    require(forward_deck != text, "missing ADJOINT block")
    forward = run(forward_deck, "final_state", np=np,
                  args=("-tran_ts_type", "beuler"), require_summary=False)

    def state(output):
        return Path(output).with_suffix(".regression").read_text().split(
            "-- SOLUTION:")[0]
    require(state(with_adjoint) == state(forward),
            "adjoint did not restore the final forward state")
    cleanup_run("final_state")
    print(f"PASS final-state restoration ({np} ranks)", flush=True)


def memory_trajectory(np):
    """The memory trajectory, which would skip the last step of each TSSolve
    (one per waypoint interval), gives the gradients of the disk one."""
    # output times and the velocity knot at 2 h make five intervals
    text = deck("dynamic_velocity").replace(
        "OUTPUT\n  TIME_UNITS d\nEND",
        "OUTPUT\n  TIME_UNITS d\n  TIMES h 1.d0 3.d0 4.d0\nEND", 1)
    require("TIMES h 1.d0" in text, "dynamic_velocity: no OUTPUT block")
    basic = run(text, "memory", np=np)
    memory = run(text, "memory", np=np, args=("-ts_trajectory_type", "memory"))
    intervals = re.findall(r"== .* TS: ", Path(memory).read_text())
    require(len(intervals) > 2, f"{len(intervals)} waypoint intervals")
    for a, b in zip(parse_summary(basic), parse_summary(memory)):
        require(a.keys() == b.keys(), "different costs or parameters")
        for key in a:
            close(a[key], b[key], rtol=1e-12, atol=0.)
    cleanup_run("memory")
    print(f"PASS memory trajectory over {len(intervals)} waypoint intervals "
          f"({np} ranks)", flush=True)


def rerun(np):
    """-ts_model_rerun_check: reruns from the recorded initial state, with
    the recorded forcing replayed, reproduce the run bit for bit at the
    deck's parameter, and at 1.1 times it give the cost and gradient of a
    run of the deck at that value, bit for bit; also with STEP velocity
    datasets, the memory trajectory, LOG_FORMULATION and adaptive steps.
    The check leaves the run's results as they were; a regime transition
    at 1.1 p0 leaves none behind."""
    cases = [("ar_1d4", ()), ("dynamic_velocity", ()),
             ("dynamic_velocity", ("-ts_trajectory_type", "memory"))]
    if np == 1:
        # one cell: not on two ranks
        cases += [("ar_log", ()),
                  ("ar_short", ("-tran_ts_adapt_type", "basic"))]
    for name, args in cases:
        text = deck(name)
        parameters, costs, gradients = parse_rerun_check(
            run(text, "rerun", np=np, args=args + ("-ts_model_rerun_check",)))
        standalone = parse_summary(run(
            set_value(text, MICROBIAL_RATE, parameters[1]), "rerun", np=np,
            args=args))
        require(standalone == (costs, gradients),
                f"{name} {args}: the rerun at 1.1 p0 gives {costs} "
                f"{gradients}, a run of the deck {standalone}")
    if np == 1:
        text = deck("ar_1d4").replace(
            "\nSUBSURFACE\n",
            "\nSUBSURFACE\nREGRESSION\n  CELL_IDS\n    1\n    4\n  /\nEND\n", 1)
        require("REGRESSION" in text, "missing SUBSURFACE")

        def results(output):
            return [line for line in Path(output).with_suffix(
                ".regression").read_text().splitlines()
                if "Time (seconds)" not in line]
        require(results(run(text, "rerun", args=("-ts_model_rerun_check",)))
                == results(run(text, "rerun")),
                "-ts_model_rerun_check changed the run's regression file")
        # gypsum depletes at 1.1 p0 (no gradient) but not at p0
        depleting = with_parameter(deck("mineral"), MINERAL_PARAMETER).replace(
            "Gypsum 1.d-2 100.d0", "Gypsum 1.5d-4 100.d0", 1)
        require("Gypsum 1.5d-4" in depleting, "missing gypsum")
        log = Path(run(depleting, "rerun",
                       args=("-ts_model_rerun_check",))).read_text()
        require("TS model rerun check at 1.1 p0: no gradient, status 1" in log
                and "TS model rerun check: PASSED" in log,
                "rerun check with a regime transition at 1.1 p0")
        no_adjoint = re.sub(r"(?ms)^ADJOINT\n.*?^END\n", "", deck("ar_short"))
        run(no_adjoint, "rerun", args=("-ts_model_rerun_check",),
            expected_error="-ts_model_rerun_check requires an ADJOINT block",
            timeout=30)
    cleanup_run("rerun")
    print(f"PASS reruns reproduce the run and runs of the deck ({np} ranks)",
          flush=True)


def regime_transitions(np):
    """A mineral regime transition (depletion, dissolution resuming) is a
    TSEvent the adjoint does not differentiate: the forward run and its
    output stand, the adjoint is skipped, and the exit code says so."""
    mineral = deck("mineral")
    depleted = mineral.replace("Gypsum 1.d-2 100.d0", "Gypsum 1.d-4 100.d0")
    if np > 1:
        # only the last two cells deplete, on a rank other than the root
        depleted = mineral.replace("NXYZ 1 1 1", "NXYZ 4 1 1")
        depleted = depleted.replace("    1.d0 1.d0 1.d0", "    4.d0 1.d0 1.d0")
        initial = re.search(r"TRANSPORT_CONDITION initial\n.*?\nEND",
                            depleted, re.S).group(0)
        second = initial.replace("initial", "depleting").replace(
            "Gypsum 1.d-2 100.d0", "Gypsum 1.d-4 100.d0")
        depleted = depleted.replace("END_SUBSURFACE", second + """
REGION downstream
  COORDINATES
    2.d0 0.d0 0.d0
    4.d0 1.d0 1.d0
  /
END
INITIAL_CONDITION
  TRANSPORT_CONDITION depleting
  REGION downstream
END
END_SUBSURFACE""")
    # dissolution resumes after precipitation of an initially absent mineral
    resumes = mineral.replace("Gypsum 1.d-2 100.d0", "Gypsum 0.d0 100.d0")
    resumes = resumes.replace("SO4--  5.d-3", "SO4--  1.d-2")
    resumes = resumes.replace("Ca++   5.d-3", "Ca++   1.d-2")
    for text in (depleted,) if np > 1 else (depleted, resumes):
        text = text.replace("OUTPUT\n  PERIODIC_OBSERVATION",
                            "OUTPUT\n  FORMAT TECPLOT POINT\n  PERIODIC TIME 5 h\n"
                            "  PERIODIC_OBSERVATION")
        output = Path(run(text, "regime", np=np,
                          expected_error="ADJOINT cannot run after",
                          expected_exit=89))
        directory = output.parent
        require("forward solve completed; adjoint skipped" in
                output.read_text(), "the skip does not say that the forward "
                "run completed")
        require("Wall Clock Time:" in output.read_text(),
                "the skip lost the timing summary")
        require(not list(directory.glob("TS-data-*")),
                "the skip leaked the trajectory directory")
        forward_deck = re.sub(r"(?ms)^ADJOINT\n.*?^END\n", "", text)
        forward = Path(run(forward_deck, "regime", np=np,
                           args=("-tran_ts_type", "beuler"),
                           require_summary=False))
        snapshots = sorted(directory.glob("case-*.tec"))
        require(len(snapshots) >= 4, "the skip lost snapshots")
        for path in snapshots + [directory / "case-obs-0.pft"]:
            require(path.read_bytes() == (forward.parent / path.name).read_bytes(),
                    f"the skip changed forward output: {path.name}")
        cleanup_run("regime")
    print(f"PASS mineral regime transitions skip the adjoint ({np} ranks)",
          flush=True)


def setup_errors(np):
    text = deck("ar_1d4")
    for args, diagnostic in (
            (("-ts_trajectory_type", "invalid_trajectory"),
             "Unknown TSTrajectory type: invalid_trajectory"),
            (("-tran_ts_type", "invalid_ts"), "Unknown TS type: invalid_ts"),
            (("-ts_trajectory_type", "visualization"),
             "ADJOINT supports the TS trajectories basic"),
            (("-ts_trajectory_solution_only",),
             "ADJOINT does not support -ts_trajectory_solution_only"),
            (("-ts_trajectory_type", "memory", "-ts_trajectory_solution_only"),
             "ADJOINT does not support -ts_trajectory_solution_only"),
            (("-ts_trajectory_type", "memory", "-ts_trajectory_max_cps_ram",
              "3"), "-ts_trajectory_max_cps_ram (checkpointing) is not"),
            (("-tran_ts_type", "bdf"), "ADJOINT requires backward Euler")):
        # the messages are checked on one process: an MPI abort can lose
        # them, so in parallel only the stop is
        output = run(text, "setup_errors", np=np, args=args,
                     expected_error=diagnostic, timeout=30,
                     message_required=(np == 1))
        require("Finished Initialization" not in
                Path(output).with_suffix(".stdout").read_text(),
                "execution continued after a setup error")
    cleanup_run("setup_errors")
    print(f"PASS invalid options and integrators stop the setup ({np} ranks)",
          flush=True)


def rejections():
    short = deck("ar_short")
    run(short.replace("MODE GIRT_DAE", "MODE GIRT"), "rejections",
        expected_error="ADJOINT block requires MODE GIRT_DAE")
    run(with_parameter(short, RATE_PARAMETER*2), "rejections",
        expected_error="PARAMETER MICROBIAL_RATE_CONSTANT[1] is given twice")
    for second in ("    COORDINATE 0.5d0 0.5d0 0.5d0\n", "    REGION middle\n"):
        run(short.replace("    REGION middle\n", "    REGION middle\n" + second, 1),
            "rejections", expected_error="takes one of COORDINATE, CELL_IDS or REGION")
    run(short.replace("REACTION_INDEX 1", "REACTION_INDEX 2"), "rejections",
        expected_error="REACTION_INDEX out of range")
    run(short.replace("KINETICS MICROBIAL", "KINETICS MINERAL\n"
                      "    MINERAL_NAME Calcite").replace(
                          "    REACTION_INDEX 1\n", ""), "rejections",
        expected_error="is not a kinetic mineral")
    two = deck("mineral").replace(
        "RATE_CONSTANT 1.d-6 mol/m^2-sec",
        "PRECIPITATION_RATE_CONSTANT 1.d-6 mol/m^2-sec\n"
        "      DISSOLUTION_RATE_CONSTANT 2.d-6 mol/m^2-sec")
    run(with_parameter(two, MINERAL_PARAMETER), "rejections",
        expected_error="requires equal precipitation and dissolution")
    run(short, "rejections", args=("-adjoint_debug_jacobianp_eps", "-1"),
        expected_error="ADJOINT debug epsilon must be positive")
    # F is linear in a rate constant: a failure needs a nonlinear parameter
    half = with_parameter(short, "  PARAMETER HALF_SATURATION_CONSTANT\n"
                          "    KINETICS MICROBIAL\n    REACTION_INDEX 1\n"
                          "    SPECIES_NAME O2(aq)\n  /\n")
    run(half, "rejections", args=("-adjoint_debug_jacobianp",
                                  "-adjoint_debug_jacobianp_eps", "1e3"),
        expected_error="ADJOINT debug check FAILED:")
    # one Richardson iteration with the exact preconditioner solves each
    # system but reports failure: Newton accepts it, TSAdjointSolve must not
    weak = ("-tran_ksp_type", "richardson", "-tran_ksp_rtol", "1e-300",
            "-tran_ksp_atol", "0", "-tran_ksp_max_it", "1",
            "-tran_snes_max_linear_solve_fail", "100000")
    run(deck("ar_1d4"), "rejections", args=weak, timeout=60,
        expected_error="ADJOINT incomplete: backward solve stopped")
    stopped = short.replace("NUMERICAL_METHODS TRANSPORT\n",
                            "NUMERICAL_METHODS TRANSPORT\n  TIMESTEPPER\n"
                            "    MAXIMUM_NUMBER_OF_TIMESTEPS 2\n  /\n")
    # no gradients: the exit code says so
    run(stopped, "rejections", expected_exit=89,
        expected_error="forward run stopped before completion (stop_flag = "
        "2); adjoint solve skipped")
    cleanup_run("rejections")
    print("PASS mode, selection, debug and solver rejections; early stop",
          flush=True)


def debug_checks():
    params = (RATE_PARAMETER +
              "  PARAMETER HALF_SATURATION_CONSTANT\n    KINETICS MICROBIAL\n"
              "    REACTION_INDEX 1\n    SPECIES_NAME O2(aq)\n  /\n"
              "  PARAMETER YIELD\n    KINETICS MICROBIAL\n"
              "    REACTION_INDEX 1\n  /\n"
              "  PARAMETER DECAY_RATE_CONSTANT\n    REACTION_INDEX 1\n  /\n")
    for text, count in ((with_parameter(deck("ar_1d4"), params), 4),
                        (with_parameter(deck("mineral"), MINERAL_PARAMETER),
                         1)):
        out = run(text, "debug", args=("-adjoint_debug_jacobianp", "-fp_trap"))
        lines = [line for line in Path(out).read_text().splitlines()
                 if "ADJOINT debug dF/dp" in line]
        require(len(lines) == count and all(float(line.split("=")[-1]) < 1e-8
                                            for line in lines),
                f"parameter Jacobian check: {lines}")
    cleanup_run("debug")
    print("PASS parameter Jacobians against central differences", flush=True)


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
        registered(test_name)
    setup_errors(np)
    gradients(np)
    selector_and_output(np)
    final_state(np)
    memory_trajectory(np)
    rerun(np)
    regime_transitions(np)
    if np == 1:
        failed_stage()
        rejections()
        debug_checks()
    print("PASS GIRT_DAE adjoint integration suite", flush=True)


if __name__ == "__main__":
    main()
