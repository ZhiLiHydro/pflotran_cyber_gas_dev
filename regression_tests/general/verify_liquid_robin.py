#!/usr/bin/env python3
"""Executable verification of a finite-volume liquid Robin boundary."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import h5py
import numpy as np


def deck(beta=1e-10, cells=1, delta=100., numerical=False, two_phase=False):
    initial = 'GAS_PRESSURE DIRICHLET\n    GAS_SATURATION DIRICHLET' if two_phase else 'LIQUID_PRESSURE DIRICHLET\n    MOLE_FRACTION DIRICHLET'
    initial_values = 'GAS_PRESSURE 101325\n  GAS_SATURATION 0.4' if two_phase else 'LIQUID_PRESSURE 200000\n  MOLE_FRACTION 0'
    return f'''SIMULATION
  SIMULATION_TYPE SUBSURFACE
  PROCESS_MODELS
    SUBSURFACE_FLOW flow
      MODE GENERAL
      OPTIONS
        ISOTHERMAL
      /
    /
  /
END
SUBSURFACE
GRID
  TYPE STRUCTURED
  NXYZ {cells} 1 1
  BOUNDS
    0 0 0
    1 1 1
  /
  GRAVITY 0 0 0
END
FLUID_PROPERTY
  PHASE LIQUID
  DIFFUSION_COEFFICIENT 1e-9
END
FLUID_PROPERTY
  PHASE GAS
  DIFFUSION_COEFFICIENT 2e-5
END
MATERIAL_PROPERTY rock
  ID 1
  POROSITY 0.25
  TORTUOSITY 1
  ROCK_DENSITY 2700
  THERMAL_CONDUCTIVITY_DRY 1
  THERMAL_CONDUCTIVITY_WET 1
  HEAT_CAPACITY 835
  PERMEABILITY
    PERM_ISO 1e-12
  /
  CHARACTERISTIC_CURVES rock
END
CHARACTERISTIC_CURVES rock
  SATURATION_FUNCTION VAN_GENUCHTEN
    LIQUID_RESIDUAL_SATURATION 0
    ALPHA 1e-4
    M 0.5
  /
  PERMEABILITY_FUNCTION MUALEM
    PHASE LIQUID
    LIQUID_RESIDUAL_SATURATION 0
    M 0.5
  /
  PERMEABILITY_FUNCTION MUALEM_VG_GAS
    PHASE GAS
    LIQUID_RESIDUAL_SATURATION 0
    GAS_RESIDUAL_SATURATION 0
    M 0.5
  /
END
NUMERICAL_METHODS FLOW
  NEWTON_SOLVER
    {'NUMERICAL_JACOBIAN' if numerical else 'ANALYTICAL_JACOBIAN'}
    ATOL 1e-14
    RTOL 1e-10
    STOL 1e-30
    MINIMUM_NEWTON_ITERATIONS 2
    MAXIMUM_NUMBER_OF_ITERATIONS 30
  /
END
TIME
  FINAL_TIME 10 s
  INITIAL_TIMESTEP_SIZE 0.01 s
  MAXIMUM_TIMESTEP_SIZE 1 s
END
OUTPUT
  ASCII_DIGITS_OF_PRECISION 15
  FORMAT HDF5
  TIMES s 0.1 1 10
  VARIABLES
    LIQUID_PRESSURE
    GAS_PRESSURE
    LIQUID_SATURATION
    GAS_SATURATION
    LIQUID_DENSITY
    GAS_DENSITY
    LIQUID_DENSITY_MOLAR
    GAS_DENSITY_MOLAR
    LIQUID_MOLE_FRACTIONS
    GAS_MOLE_FRACTIONS
    LIQUID_VISCOSITY
    TEMPERATURE
    THERMODYNAMIC_STATE
  /
  MASS_BALANCE_FILE
    PERIODIC TIMESTEP 1
  /
END
REGION all
  INFINITE
END
REGION inlet
  FACE WEST
  COORDINATES
    0 0 0
    0 1 1
  /
END
REGION outlet
  FACE EAST
  COORDINATES
    1 0 0
    1 1 1
  /
END
FLOW_CONDITION initial
  TYPE
    {initial}
    TEMPERATURE DIRICHLET
  /
  {initial_values}
  TEMPERATURE 12
END
FLOW_CONDITION inlet
  TYPE
    LIQUID_PRESSURE DIRICHLET
    MOLE_FRACTION DIRICHLET
    TEMPERATURE DIRICHLET
  /
  LIQUID_PRESSURE {200000+delta}
  MOLE_FRACTION 0
  TEMPERATURE 12
  {'' if beta is None else 'LIQUID_ROBIN_COEFFICIENT '+str(beta)}
END
FLOW_CONDITION outlet
  TYPE
    LIQUID_PRESSURE DIRICHLET
    MOLE_FRACTION DIRICHLET
    TEMPERATURE DIRICHLET
  /
  LIQUID_PRESSURE 200000
  MOLE_FRACTION 0
  TEMPERATURE 12
END
INITIAL_CONDITION
  FLOW_CONDITION initial
  REGION all
END
STRATA
  REGION all
  MATERIAL rock
END
BOUNDARY_CONDITION inlet
  FLOW_CONDITION inlet
  REGION inlet
END
BOUNDARY_CONDITION outlet
  FLOW_CONDITION outlet
  REGION outlet
END
END_SUBSURFACE
'''


def read_fields(path):
    with h5py.File(path) as h:
        groups = sorted([k for k in h if 'Time' in k], key=lambda k: float(re.search('Time\\s*:?\\s*([0-9.eE+-]+)', k).group(1)))
        return [{name.replace('_', ' ', 1) if not name.startswith('X_') else name: h[k][name][:] for name in h[k]} for k in groups]

def main():
    if not __debug__:
        raise RuntimeError('Run this verifier without Python -O: assertions are required.')
    p = argparse.ArgumentParser()
    p.add_argument('--executable', type=Path, required=True)
    p.add_argument('--reference-executable', type=Path)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    out = args.output.resolve()
    out.mkdir(exist_ok=False, parents=True)
    exe = args.executable.resolve()
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    report = {'executable_sha256': hashlib.sha256(exe.read_bytes()).hexdigest(),
              'tolerances': {'mass_balance_normalized': 1e-7,
                             'steady_rate_relative': 1e-4,
                             'active_jacobian_relative': 1e-5}}

    def run(name, text, executable=exe, invalid=False, flags=()):
        folder = out / name
        folder.mkdir()
        (folder / 'test.in').write_text(text)
        with (folder / 'run.log').open('w') as log:
            result = subprocess.run([str(executable), '-input_prefix', 'test', *flags], cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=180)
        log = (folder / 'run.log').read_text()
        success = result.returncode == 0 and 'Wall Clock Time:' in log and ('Stopping!' not in log)
        if invalid:
            assert not success, name
            assert 'Liquid Robin' in log or 'LIQUID_ROBIN_COEFFICIENT' in log, (name, log[-2000:])
            report[name] = {'rejected': True}
            return None
        assert success, (name, result.returncode, log[-2000:])
        fields = []
        for path in sorted(folder.glob('test*.h5')):
            fields.extend(read_fields(path))
        assert fields, name
        with (folder / 'test-mas.dat').open() as f:
            header = next(csv.reader(f, skipinitialspace=True))
        data = np.loadtxt(folder / 'test-mas.dat', skiprows=1, ndmin=2)
        first = fields[0]
        initial = []
        for component, mw in [('l', 18.01534), ('g', 28.96)]:
            mass = 0
            for phase, letter in [('Liquid', 'l'), ('Gas', 'g')]:
                mass += np.sum(first[f'{phase} Saturation'] * first[f'{phase} Density [kmol_m^3]'] * first[f'X_{component}^{letter}']) * mw * 0.25 / first['Liquid Saturation'].size
            initial.append(mass)
        errors = []
        for j, component in enumerate(['Water', 'Air']):
            inventory = data[:, [i for i, h in enumerate(header) if h.startswith('Global ' + component + ' Mass')]].sum(axis=1)
            flux = data[:, [i for i, h in enumerate(header) if not h.startswith('Global ') and h.endswith(component + ' Mass [kg]')]].sum(axis=1)
            error = np.abs(inventory - initial[j] - flux).max() / max(initial[j], 1.0)
            assert error < 1e-07, (name, component, error)
            errors.append(float(error))
        assert all((np.isfinite(v).all() for field in fields for v in field.values()))
        for field in fields:
            np.testing.assert_allclose(field['Liquid Saturation'] + field['Gas Saturation'],
                                       1, rtol=0, atol=1e-12)
        rate_idx = next((i for i, h in enumerate(header) if h.startswith('inlet Water Mass [kg/')))
        unit = header[rate_idx].split('[kg/')[1][:-1]
        factor = {'s': 1, 'd': 86400, 'y': 365 * 86400}[unit]
        rate = float(data[-1, rate_idx] / factor)
        report[name] = {'water_air_balance_errors': errors, 'inlet_water_rate_kg_s': rate}
        print(name, 'PASS', flush=True)
        return (fields, rate)
    results = {}
    for name, beta, n, delta in [('native', None, 1, 100), ('finite', 1e-10, 1, 100), ('small_coefficient', 1e-17, 1, 100), ('refined', 1e-10, 4, 100), ('zero', 0, 1, 100), ('equal', 1e-10, 1, 0), ('reverse', 1e-10, 1, -100), ('large', 10000000000.0, 1, 100)]:
        fields, rate = run(name, deck(beta, n, delta))
        results[name] = (fields, rate)
        if name not in ['zero', 'equal']:
            viscosity = next((v for k, v in fields[-1].items() if 'Viscosity' in k))
            mu = float(viscosity.mean())
            density = float(fields[-1]['Liquid Density [kg_m^3]'].mean())
            expected = delta / (mu / 1e-12 + (0 if beta is None else 1 / beta)) * density
            rel = abs(rate - expected) / abs(expected)
            assert rel < 0.0001, (name, rate, expected, mu, rel)
            report[name].update(expected_rate_kg_s=expected, relative_error_to_series_solution=rel)
        else:
            assert abs(rate) < 1e-14, (name, rate)
    for a, b in [('native', 'large')]:
        assert len(results[a][0]) == len(results[b][0])
        for fa, fb in zip(results[a][0], results[b][0]):
            assert fa.keys() == fb.keys()
            for field in fa:
                np.testing.assert_array_equal(fa[field], fb[field])
    assert abs(results['finite'][1] / results['refined'][1] - 1) < 0.0001
    hydro = deck().replace('FLOW_CONDITION inlet\n  TYPE\n    LIQUID_PRESSURE DIRICHLET', 'FLOW_CONDITION inlet\n  TYPE\n    LIQUID_PRESSURE HYDROSTATIC').replace('  LIQUID_PRESSURE 200100', '  DATUM 0 0 0.5\n  LIQUID_PRESSURE 200100').replace('GRAVITY 0 0 0', 'GRAVITY 0 0 -9.81')
    _, hydro_rate = run('hydrostatic', hydro)
    np.testing.assert_allclose(hydro_rate, results['finite'][1], rtol=1e-8, atol=0)
    run('negative_rejected', deck(-1), invalid=True)
    run('nan_rejected', deck('NaN'), invalid=True)
    run('infinity_rejected', deck('Infinity'), invalid=True)
    run('duplicate_rejected', deck().replace('LIQUID_ROBIN_COEFFICIENT 1e-10', 'LIQUID_ROBIN_COEFFICIENT 1e-10\n  LIQUID_ROBIN_COEFFICIENT 1e-10'), invalid=True)
    run('unsupported_type_rejected', deck().replace('FLOW_CONDITION inlet\n  TYPE\n    LIQUID_PRESSURE DIRICHLET',
                                           'FLOW_CONDITION inlet\n  TYPE\n    LIQUID_PRESSURE SEEPAGE'), invalid=True)
    run('initial_condition_rejected', deck().replace('INITIAL_CONDITION\n  FLOW_CONDITION initial', 'INITIAL_CONDITION\n  FLOW_CONDITION inlet'), invalid=True)
    run('non_darcy_rejected', deck().replace('        ISOTHERMAL', '        ISOTHERMAL\n        NON_DARCY_FLOW'), invalid=True)
    ta, _ = run('two_phase_analytic', deck(two_phase=True), flags=('-flow_snes_test_jacobian', '-flow_snes_test_jacobian_view'))
    tn, _ = run('two_phase_numeric', deck(two_phase=True, numerical=True))
    for key in ['Liquid Saturation', 'Gas Saturation', 'Liquid Pressure [Pa]', 'Gas Pressure [Pa]']:
        np.testing.assert_allclose(ta[-1][key], tn[-1][key], rtol=1e-07, atol=1e-08)
    draining = lambda numerical: deck(two_phase=True, numerical=numerical).replace('LIQUID_PRESSURE 200100', 'LIQUID_PRESSURE 80000').replace('LIQUID_PRESSURE 200000', 'LIQUID_PRESSURE 80000')
    da, _ = run('draining_analytic', draining(False), flags=('-flow_snes_test_jacobian', '-flow_snes_test_jacobian_view'))
    dn, _ = run('draining_numeric', draining(True))
    for key in ['Liquid Saturation', 'Gas Saturation', 'Liquid Pressure [Pa]', 'Gas Pressure [Pa]']:
        np.testing.assert_allclose(da[-1][key], dn[-1][key], rtol=1e-07, atol=1e-08)
    # Exercise heat advection and temperature-dependent mobility as well.
    thermal = deck().replace('        ISOTHERMAL', '')
    inlet_start = thermal.index('FLOW_CONDITION inlet\n')
    outlet_start = thermal.index('FLOW_CONDITION outlet\n')
    thermal = (thermal[:inlet_start]
               + thermal[inlet_start:outlet_start].replace('TEMPERATURE 12', 'TEMPERATURE 18')
               + thermal[outlet_start:])
    thermal_analytic, _ = run('nonisothermal_analytic', thermal)
    thermal_numeric, _ = run('nonisothermal_numeric',
                              thermal.replace('ANALYTICAL_JACOBIAN', 'NUMERICAL_JACOBIAN'))
    assert thermal_analytic[-1].keys() == thermal_numeric[-1].keys()
    for key in thermal_analytic[-1]:
        np.testing.assert_allclose(thermal_analytic[-1][key],
                                   thermal_numeric[-1][key], rtol=1e-7, atol=1e-8)
    for name in ['two_phase_analytic', 'draining_analytic']:
        log = (out / name / 'run.log').read_text()
        matrices = []
        for title in ['Hand-coded Jacobian', 'Finite difference Jacobian']:
            blocks = re.findall(re.escape(title) + ' -+\\n(?:[^\\n]*\\n){2}row 0:([^\\n]*)\\nrow 1:([^\\n]*)', log)
            assert blocks, (name, title)
            matrices.append([np.array([[float(v) for col, v in re.findall('\\((\\d+), ([^)]*)\\)', row) if int(col) < 2] for row in block]) for block in blocks])
        assert len(matrices[0]) == len(matrices[1])
        error = max((np.linalg.norm(a - b) / max(np.linalg.norm(a), 1e-30) for a, b in zip(*matrices)))
        assert error < 1e-05, (name, error)
        report[name]['max_active_jacobian_relative_difference'] = float(error)
    drain = np.loadtxt(out / 'draining_analytic/test-mas.dat', skiprows=1)
    assert drain[0, 8] < 0 and drain[-1, 8] > 0, 'Must exercise both upwind directions'
    if args.reference_executable:
        report['reference_executable_sha256'] = hashlib.sha256(
            args.reference_executable.resolve().read_bytes()).hexdigest()
        old, _ = run('old_native', deck(None), args.reference_executable.resolve())
        assert len(results['native'][0]) == len(old)
        for fa, fb in zip(results['native'][0], old):
            assert fa.keys() == fb.keys()
            for field in fa:
                np.testing.assert_array_equal(fa[field], fb[field])
        report['unchanged_native_fields_bitwise'] = True
    report['all_passed'] = True
    (out / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
if __name__ == '__main__':
    main()
