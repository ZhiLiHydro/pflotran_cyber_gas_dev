#!/usr/bin/env python3
"""Checks of MODE RICHARDS_DAE (Richards flow as a TS equation model).

Post-processor of richards_dae, or run directly:
python3 check_richards_dae.py --executable <pflotran>

On fixed backward Euler steps the model reproduces classic RICHARDS to 1e-8
(their Newton tolerances) with every boundary and source type it accepts
but the HET ones (dataset files; only their rates' source differs), with
constant, STEP and LINEAR boundary data and with soil compressibility about
the initial pressure. Classic RICHARDS takes
each step's forcing at the step start, so LINEAR data are compared with its
knots one step earlier. Also: STEP data hold until their knot next to LINEAR
data; the Jacobian entry by entry (away from hydrostatic equilibrium, where
the upwinded relative permeability has a kink); classic Newton options are
rejected; a BDF order fallback; two processes against one.
"""

import argparse
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "adjoint"))
from adjoint_verify import mpiexec  # noqa: E402
FIXED = ("-flow_ts_type", "beuler", "-flow_ts_adapt_type", "none",
         "-ts_model_newton_kappa", "1e-4")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def run(text, args=(), np=1, expected_error=None):
    executable = Path(os.environ["PFLOTRAN_TEST_EXECUTABLE"]).resolve()
    launcher = []
    if np > 1:
        mpi = mpiexec()
        require(mpi, "MPIEXEC is required for the parallel check")
        launcher = shlex.split(mpi) + ["-n", str(np)]
        text = text.replace("    SOLVER DIRECT\n", "    SOLVER ITERATIVE\n")
    with tempfile.TemporaryDirectory(prefix="pflotran-richards-dae-") as name:
        work = Path(name)
        (work / "case.in").write_text(text)
        process = subprocess.run(launcher + [str(executable), "-pflotranin",
                                             "case.in"] + list(args),
                                 cwd=work, capture_output=True, text=True,
                                 timeout=600)
        output = process.stdout + process.stderr
        if expected_error is not None:
            require(process.returncode != 0 and expected_error in output,
                    f"expected {expected_error!r} ({process.returncode}):\n"
                    f"{output[-3000:]}")
            return None, output
        require(process.returncode == 0 and "PETSC ERROR" not in output,
                f"PFLOTRAN failed ({process.returncode}):\n{output[-3000:]}")
        values = {}
        field = None
        for line in (work / "case.regression").read_text().splitlines():
            if line.startswith("--"):
                field = line.strip("- ")
                continue
            if field and "SOLUTION" not in field and ":" in line:
                key, value = line.split(":", 1)
                values[field, key.strip()] = float(value)
    return values, output


def same(a, b, rtol, what):
    require(a.keys() == b.keys() and a, f"{what}: fields differ")
    for key, value in a.items():
        require(abs(value - b[key]) <= rtol*max(abs(value), abs(b[key])),
                f"{what}: {key} {value:.16e} versus {b[key]:.16e}")


def classic(text):
    # constant steps: no Newton-count governor
    return text.replace("MODE RICHARDS_DAE", "MODE RICHARDS").replace(
        "NUMERICAL_METHODS FLOW\n",
        "NUMERICAL_METHODS FLOW\n  TIMESTEPPER\n    TS_ACCELERATION 100\n"
        "  /\n", 1)


def transient(text, interpolation, shift):
    """The bottom pressure and the recharge as lists with knots shifted
    earlier by shift [y]; only the flux knots are classic waypoints."""
    def t(value):
        return f"{max(value - shift, 0.):.6f}d0"
    text = text.replace(
        "FLOW_CONDITION bottom\n  TYPE\n    LIQUID_PRESSURE DIRICHLET\n  /\n"
        "  LIQUID_PRESSURE 110000.d0\n",
        "FLOW_CONDITION bottom\n  TYPE\n    LIQUID_PRESSURE DIRICHLET\n  /\n"
        f"  INTERPOLATION {interpolation}\n  LIQUID_PRESSURE LIST\n"
        "    TIME_UNITS y\n    DATA_UNITS Pa\n"
        f"    0.d0 {101325. + (3675. * shift / 0.025 if interpolation == 'LINEAR' else 0.):.4f}d0\n"
        f"    {t(0.025)} 105000.d0\n    1.d0 105000.d0\n  /\n")
    text = text.replace(
        "  LIQUID_FLUX 10.d0 cm/y\n",
        f"  INTERPOLATION {interpolation}\n  LIQUID_FLUX LIST\n"
        "    TIME_UNITS y\n    DATA_UNITS m/y\n"
        f"    0.d0 {(2. * shift / 0.02 if interpolation == 'LINEAR' else 0.):.6f}d0\n"
        f"    {t(0.02)} 2.d0\n    {t(0.03)} 0.5d0\n    1.d0 0.5d0\n  /\n")
    return text


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable")
    options = parser.parse_args()
    if options.executable:
        os.environ["PFLOTRAN_TEST_EXECUTABLE"] = \
            str(Path(options.executable).resolve())
    deck = (HERE / "richards_dae.in").read_text()

    dae, _ = run(deck, FIXED)
    same(dae, run(classic(deck))[0], 1e-8, "constant forcing")
    print("PASS backward Euler reproduces classic RICHARDS", flush=True)

    # every other boundary type the model accepts, at the bottom
    bottom = ("FLOW_CONDITION bottom\n  TYPE\n    LIQUID_PRESSURE DIRICHLET\n"
              "  /\n  LIQUID_PRESSURE 110000.d0\n")
    datum = "  DATUM 0.d0 0.d0 1.d0\n  LIQUID_PRESSURE 101325.d0\n"
    boundaries = {
        "HYDROSTATIC": datum, "SEEPAGE": datum,
        "CONDUCTANCE": datum + "  CONDUCTANCE 1.d-12\n",
        "DIRICHLET_SEEPAGE": "  LIQUID_PRESSURE 110000.d0\n",
        "DIRICHLET_CONDUCTANCE": "  LIQUID_PRESSURE 110000.d0\n"
                                 "  CONDUCTANCE 1.d-12\n",
        "UNIT_GRADIENT": "  LIQUID_PRESSURE 101325.d0\n",
        "ZERO_GRADIENT": "  LIQUID_PRESSURE 101325.d0\n"}
    require(bottom in deck, "missing bottom condition")
    for kind, values in boundaries.items():
        variant = deck.replace(bottom, "FLOW_CONDITION bottom\n  TYPE\n"
                               f"    LIQUID_PRESSURE {kind}\n  /\n" + values)
        same(run(variant, FIXED)[0], run(classic(variant))[0], 1e-8, kind)
    # and every source type, at mid column
    require("\nEND_SUBSURFACE" in deck, "missing END_SUBSURFACE")
    sources = {"MASS_RATE": "1.d-6 kg/s", "VOLUMETRIC_RATE": "1.d-9 m^3/s",
               "SCALED_MASS_RATE VOLUME": "1.d-6 kg/s",
               "SCALED_VOLUMETRIC_RATE VOLUME": "1.d-9 m^3/s"}
    for kind, rate in sources.items():
        variant = deck.replace("\nEND_SUBSURFACE", f"""
REGION middle
  COORDINATE 0.5d0 0.5d0 5.1d0
END

FLOW_CONDITION well
  TYPE
    RATE {kind}
  /
  RATE {rate}
END

SOURCE_SINK well
  FLOW_CONDITION well
  REGION middle
END
END_SUBSURFACE""")
        same(run(variant, FIXED)[0], run(classic(variant))[0], 1e-8, kind)
    print(f"PASS boundary types {', '.join(boundaries)} and source types "
          f"{', '.join(sources)} reproduce classic RICHARDS", flush=True)

    # the pressure knot (0.025 y) is no flux knot, output or step boundary
    # of its own: the model adds it as a waypoint
    step = transient(deck, "STEP", 0.)
    same(run(step, FIXED)[0], run(classic(step))[0], 1e-8, "STEP forcing")
    print("PASS STEP boundary data (a knot inside a waypoint interval)",
          flush=True)

    linear = transient(deck, "LINEAR", 0.)
    shifted = transient(deck, "LINEAR", 0.005)
    same(run(linear, FIXED)[0], run(classic(shifted))[0], 1e-8,
         "LINEAR forcing")
    print("PASS LINEAR boundary data (classic RICHARDS one step earlier)",
          flush=True)

    # a LINEAR recharge list of constant value (the same physics) next to
    # a STEP pressure: the pressure still holds until its knot
    step_pressure = deck.replace(
        "  LIQUID_PRESSURE 110000.d0\n",
        "  INTERPOLATION STEP\n  LIQUID_PRESSURE LIST\n    TIME_UNITS y\n"
        "    DATA_UNITS Pa\n    0.d0 110000.d0\n    0.025d0 105000.d0\n"
        "    1.d0 105000.d0\n  /\n")
    mixed = step_pressure.replace(
        "  LIQUID_FLUX 10.d0 cm/y\n",
        "  INTERPOLATION LINEAR\n  LIQUID_FLUX LIST\n    TIME_UNITS y\n"
        "    DATA_UNITS m/y\n    0.d0 0.1d0\n    1.d0 0.1d0\n  /\n")
    require(deck != step_pressure != mixed, "missing pressure or recharge")
    same(run(mixed, FIXED)[0], run(step_pressure, FIXED)[0], 1e-12,
         "STEP pressure next to LINEAR recharge")
    print("PASS STEP data hold until their knot next to LINEAR data",
          flush=True)

    compressible = deck.replace(
        "  POROSITY 0.25d0\n", "  POROSITY 0.25d0\n"
        "  SOIL_COMPRESSIBILITY_FUNCTION LEIJNSE\n"
        "  SOIL_COMPRESSIBILITY 1.d-7\n"
        "  SOIL_REFERENCE_PRESSURE INITIAL_PRESSURE\n")
    require(compressible != deck, "missing porosity")
    same(run(compressible, FIXED)[0], run(classic(compressible))[0], 1e-8,
         "soil compressibility")
    print("PASS soil compressibility about the initial pressure", flush=True)

    _, log = run(deck, FIXED + ("-ts_model_contract_test",
                                "-ts_model_column_check"))
    contract = re.findall(r"TS model contract: .*: (PASS|FAIL)", log)
    require(len(contract) == 5 and set(contract) == {"PASS"},
            f"model contract {contract}")
    match = re.search(r"column check:\s*(\d+) entries off", log)
    require(match and int(match.group(1)) == 0, "Jacobian column check")
    print("PASS model contract and Jacobian entry by entry", flush=True)

    # classic Newton options that PETSc TS and the model replace
    newton = deck.replace(
        "NUMERICAL_METHODS FLOW\n", "NUMERICAL_METHODS FLOW\n"
        "  NEWTON_SOLVER\n    PRESSURE_CHANGE_LIMIT 1.d5\n"
        "    ITOL_RELATIVE_UPDATE 1.d-6\n  /\n", 1)
    require(newton != deck, "missing NUMERICAL_METHODS FLOW")
    run(newton, expected_error="does not support the classic RICHARDS Newton "
        "solver options PRESSURE_CHANGE_LIMIT, UPDATE_INF_TOL/")
    print("PASS classic Newton options are rejected", flush=True)

    # a step failure at BDF order 2 (two rejections allowed) redoes the
    # rest of the interval at order 1 with the interval's forcing
    _, log = run(deck, ("-flow_ts_max_reject", "2"))
    match = re.search(r"BDF order fallbacks\s+(\d+)", log)
    require(match and int(match.group(1)) > 0, "no BDF order fallback")
    print("PASS a BDF order fallback redoes the rest of an interval",
          flush=True)

    if mpiexec():
        iterative = FIXED + ("-flow_ksp_type", "gmres", "-flow_ksp_rtol",
                             "1e-12")
        same(run(deck, iterative, np=2)[0],
             run(deck.replace("    SOLVER DIRECT\n", "    SOLVER ITERATIVE\n"),
                 iterative)[0], 1e-9, "two processes")
        print("PASS two processes match one", flush=True)
    else:
        print("SKIP two processes: no MPI launcher", flush=True)
    print("PASS RICHARDS_DAE checks", flush=True)


if __name__ == "__main__":
    main()
