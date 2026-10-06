#!/usr/bin/env python3
"""Adjoint application decks (APPLICATIONS; each deck's header describes
it) under MODE GIRT_DAE, alone and with MODE RICHARDS_DAE flow.

Post-processor of adjoint-applications and adjoint-applications-np2, or run
directly: python3 check_applications.py --executable <pflotran> [--np 2]
(parallel runs need PFLOTRAN_TEST_MPIEXEC, PETSc's mpiexec or one on PATH).

Gradients match central differences, the memory trajectory matches the disk
one, and two processes match one.
"""

import argparse
import os
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from adjoint_verify import parse_summary, cleanup_run
from check_coupled import deck, run, set_value, require, close

RATE = r"(MICROBIAL_REACTION.*?RATE_CONSTANT\s+)(\S+)"
PERMEABILITY = r"(PERM_ISO\s+)(\S+)"
# parameter, deck value, relative steps, tolerance; steps where the
# difference is resolved above the Newton tolerance and below truncation
APPLICATIONS = {
    "abcd_batch": [(RATE, 1e-6, (1e-3, 1e-4), 1e-7)],
    "cylindrical_coupled": [(RATE, 1e-6, (1e-3, 1e-4), 1e-7),
                            (PERMEABILITY, 1e-12, None, 2e-5)],
    "adjoint-applications": [(PERMEABILITY, 1e-12, (1e-4, 1e-5), 1e-7)],
}
# the radial permeability changes J only through the water density (the
# injection sets the flow), too little for small steps: Richardson-
# extrapolate two larger steps at a tighter Newton tolerance
RICHARDSON = (1e-2, 5e-3)
TIGHT = ("-ts_model_newton_kappa", "1e-8")


def central(text, pattern, value, relative, tight=False):
    h = value * relative
    kwargs = dict(args=TIGHT, kappa=False) if tight else {}
    plus = parse_summary(run(set_value(text, pattern, value + h), "fd",
                             **kwargs))[0]
    minus = parse_summary(run(set_value(text, pattern, value - h), "fd",
                              **kwargs))[0]
    return {cost: (plus[cost] - minus[cost]) / (2 * h) for cost in plus}


def gradients(name):
    text = deck(name)
    costs, grads = parse_summary(run(text, "fd"))
    for column, (pattern, value, steps, rtol) in enumerate(APPLICATIONS[name],
                                                           start=1):
        for cost in costs:
            require(abs(value * grads[cost, column]) >
                    1e-10 * abs(costs[cost]),
                    f"{name}: parameter {column} has no measurable effect")
        if steps is None:
            reference = parse_summary(run(text, "fd", args=TIGHT,
                                          kappa=False))[1]
            coarse, fine = (central(text, pattern, value, h, tight=True)
                            for h in RICHARDSON)
            for cost in costs:
                extrapolated = (4 * fine[cost] - coarse[cost]) / 3
                close(extrapolated, reference[cost, column], rtol=rtol)
            continue
        for relative in steps:
            fd = central(text, pattern, value, relative)
            for cost in costs:
                close(fd[cost], grads[cost, column], rtol=rtol)
    cleanup_run("fd")
    print(f"PASS {name}: {len(APPLICATIONS[name])} parameter(s) against "
          "central differences", flush=True)


def memory_trajectory(name, np):
    text = deck(name)
    basic = parse_summary(run(text, "memory", np=np))
    output = run(text, "memory", np=np, args=("-ts_trajectory_type", "memory"))
    intervals = len(re.findall(r"== .* TS: ", Path(output).read_text()))
    require(intervals > 1, f"{name}: {intervals} waypoint interval(s)")
    for a, b in zip(basic, parse_summary(output)):
        for key in a:
            close(a[key], b[key], rtol=1e-12, atol=0.)
    cleanup_run("memory")
    print(f"PASS {name}: memory trajectory over {intervals} intervals ({np} "
          "ranks)", flush=True)


def parallel(name, output=None):
    """Two processes against one; the radial permeability gradient, a small
    difference of large terms, depends on the linear solver at 1e-6."""
    text = deck(name)
    serial = parse_summary(run(text, "parallel"))
    if output is None:
        output = run(text, "parallel", np=2)
    for one, two in zip(serial, parse_summary(output)):
        for key in one:
            loose = name == "cylindrical_coupled" and key == (1, 2)
            close(one[key], two[key], rtol=1e-5 if loose else 1e-9)
    cleanup_run("parallel")
    print(f"PASS {name}: two ranks match one", flush=True)


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
    registered = None
    if test_name:
        registered = Path(test_name + ".out")
        parse_summary(registered)
    if np == 1:
        for name in APPLICATIONS:
            gradients(name)
        memory_trajectory("cylindrical_coupled", 1)
    else:
        # the one-cell batch cannot be partitioned
        parallel("cylindrical_coupled")
        parallel("adjoint-applications", output=registered)
        memory_trajectory("cylindrical_coupled", np)
    print("PASS adjoint applications", flush=True)


if __name__ == "__main__":
    main()
