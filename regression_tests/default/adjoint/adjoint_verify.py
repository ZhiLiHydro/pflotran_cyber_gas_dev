"""Isolated deck runs and finite-value checks shared by the adjoint checks."""

import math
import os
from pathlib import Path
import re
import signal
import shlex
import shutil
import struct
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
_RUNS = {}


def finite(value):
    result = float(value)
    if not math.isfinite(result):
        raise RuntimeError(f"nonfinite adjoint value: {value}")
    return result


def parse_summary(path):
    text = Path(path).read_text()
    match = re.search(r"=== ADJOINT SUMMARY ===(.*?)=== END ADJOINT ===",
                      text, re.S)
    if not match:
        raise RuntimeError(f"missing complete adjoint summary: {path}")
    summary = match.group(1)
    costs = {int(i): finite(v) for i, v in
             re.findall(r"J_(\d+)\(T\)\s*=\s*(\S+)", summary)}
    gradients = {(int(i), int(j)): finite(v) for i, j, v in re.findall(
        r"dJ_(\d+)/dp_(\d+)\s*\[.*?\]\s*=\s*(\S+)", summary)}
    if not costs or not gradients:
        raise RuntimeError(f"missing costs or parameter gradients: {path}")
    for value in re.findall(r"max \|dJ_\d+/du_phys\(0\)\|\s*=\s*(\S+)", summary):
        finite(value)
    return costs, gradients


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def set_value(text, pattern, value):
    """The first match of pattern (a group, then the value) set to value."""
    text, count = re.subn(pattern, lambda m: m.group(1) + f"{value:.17e}",
                          text, count=1, flags=re.S)
    require(count == 1, f"missing {pattern}")
    return text


def deck(name):
    return load_deck(HERE / f"{name}.in")


def load_deck(path):
    path = Path(path).resolve()
    text = path.read_text()
    # Preserve the meaning of relative database paths in an isolated run.
    def database(match):
        value = Path(match.group(2))
        if not value.is_absolute():
            value = path.parent / value
        return match.group(1) + str(value.resolve())
    return re.sub(r"(?m)^(\s*DATABASE\s+)(\S+)", database, text)


def mpiexec():
    """The test manager's launcher, else PETSc's, else the one on the path
    (None if none)."""
    mpi = os.environ.get("PFLOTRAN_TEST_MPIEXEC") or os.environ.get("MPIEXEC")
    petsc = Path(os.environ.get("PETSC_DIR", "/nonexistent"),
                 os.environ.get("PETSC_ARCH", ""), "bin", "mpiexec")
    return mpi or (str(petsc) if petsc.is_file() else shutil.which("mpiexec"))


def run_deck(text, tag, args=(), np=1,
             full_field=False, expected_error=None, require_summary=True,
             expected_exit=None, timeout=180, message_required=True):
    executable = Path(os.environ.get("PFLOTRAN_TEST_EXECUTABLE",
                      os.environ.get("PFLOTRAN", ROOT / "src/pflotran/pflotran")))
    executable = executable.resolve()
    launcher = []
    if np > 1:
        mpi = mpiexec()
        if not mpi:
            raise RuntimeError("MPIEXEC is required for parallel adjoint checks")
        launcher = shlex.split(mpi) + ["-n", str(np)]
    directory = tempfile.TemporaryDirectory(prefix="pflotran-adjoint-")
    # Repeated tags are allowed, but never reuse a previous output file.
    _RUNS.setdefault(tag, []).append(directory)
    work = Path(directory.name)
    if full_field:
        text, count = re.subn(r"(?m)^ADJOINT\s*$",
                             "ADJOINT\n  OUTPUT\n    FULL_FIELD\n  /",
                             text, count=1)
        if count != 1:
            raise RuntimeError("missing ADJOINT block for full-field output")
    (work / "case.in").write_text(text)
    command = launcher + [str(executable), "-pflotranin", "case.in",
                           "-successful_exit_code", "0",
                           "-tran_ts_adapt_type", "none"] + list(args)
    process = subprocess.Popen(command, cwd=work, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True, start_new_session=True)
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        # Kill the MPI launcher and its children, not just the launcher.
        os.killpg(process.pid, signal.SIGKILL)
        process.communicate()
        raise
    combined = stdout + stderr
    (work / "case.stdout").write_text(combined)
    output = work / "case.out"
    if expected_exit is not None and process.returncode != expected_exit:
        raise RuntimeError(f"expected exit={expected_exit}, got "
                           f"exit={process.returncode}\n{combined[-5000:]}")
    if expected_error is not None:
        # an MPI abort can lose the message: message_required=False only
        # checks that the run stopped
        missing = message_required and expected_error not in combined
        if (process.returncode == 0 or missing or
                "=== ADJOINT SUMMARY ===" in combined or
                "Caught signal" in combined):
            raise RuntimeError(f"expected diagnostic {expected_error!r}, "
                               f"exit={process.returncode}\n{combined[-5000:]}")
    elif process.returncode != 0 or "PETSC ERROR" in combined:
        raise RuntimeError(f"PFLOTRAN failed: exit={process.returncode}\n"
                           f"{combined[-5000:]}")
    if require_summary and expected_error is None:
        parse_summary(output)
    return str(output)


def cleanup_run(tag):
    for directory in _RUNS.pop(tag, []):
        directory.cleanup()


def read_sensitivities(output, ncomp, cost=1, suffix="du0"):
    """Read a real-double PETSc binary Vec with 32- or 64-bit indices."""
    path = Path(output).with_name(f"case-adjoint-dJ{cost}-{suffix}.bin")
    raw = path.read_bytes()
    for integer in (">i", ">q"):
        width = struct.calcsize(integer)
        if len(raw) < 2 * width:
            continue
        classid = struct.unpack_from(integer, raw)[0]
        size = struct.unpack_from(integer, raw, width)[0]
        if classid == 1211214 and size >= 0 and len(raw) == 2 * width + 8 * size:
            values = struct.unpack_from(f">{size}d", raw, 2 * width)
            if size % ncomp:
                raise RuntimeError("sensitivity size is not a multiple of ncomp")
            return {(i // ncomp + 1, i % ncomp + 1): finite(v)
                    for i, v in enumerate(values)}
    raise RuntimeError(f"unsupported PETSc sensitivity Vec format: {path}")
