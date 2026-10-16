#!/usr/bin/env python3
"""Finite-difference checks of the MODE GIRT_DAE adjoint for several
microbial reactions, Monod, inhibition, yield and decay parameters, and the
OBSERVATIONS objective. Post-processor of adjoint-extensions and
adjoint-extensions-np2; see check_adjoint.py for direct use."""

import argparse
import os
from pathlib import Path
import re
import tempfile

from check_adjoint import (close, deck, parse_summary,
                           require, cleanup_run, read_sensitivities)
from adjoint_verify import parse_rerun_check
from check_adjoint import run as run_adjoint


def run(text, tag, np=1, **kwargs):
    # converge Newton below the finite-difference and formulation
    # comparisons; the model's default test is too loose for them
    kwargs["args"] = ("-ts_model_newton_kappa", "1e-4") + \
        tuple(kwargs.get("args", ()))
    return run_adjoint(text, tag, np=np, **kwargs)


def example(name):
    return deck("ar_" + name)


def log_formulation(text):
    # the first Newton update of a trace species in log coordinates raises
    # the residual by orders of magnitude before it converges
    text = text.replace("CHEMISTRY\n", "CHEMISTRY\n  LOG_FORMULATION\n")
    return text.replace("NUMERICAL_METHODS TRANSPORT\n",
                        "NUMERICAL_METHODS TRANSPORT\n  NEWTON_SOLVER\n"
                        "    DTOL 1.d20\n  /\n", 1)


def replace_nth(text, pattern, value, index=0):
    matches = list(re.finditer(pattern, text, re.M))
    require(len(matches) > index, f"missing occurrence {index} of {pattern}")
    match = matches[index]
    return text[:match.start(1)] + str(value) + text[match.end(1):]


def parameter(kind, reaction=1, species=None):
    kinetics = '' if kind == 'DECAY_RATE_CONSTANT' else '    KINETICS MICROBIAL\n'
    species = f'    SPECIES_NAME {species}\n' if species else ''
    return (f'  PARAMETER {kind}\n{kinetics}    REACTION_INDEX {reaction}\n'
            f'{species}  /\n')


def with_parameters(text, parameters):
    return re.sub(r'  PARAMETER RATE_CONSTANT\n.*?  /\n',
                  lambda _: ''.join(parameters), text, count=1, flags=re.S)


def gradient_check(text, parameters, tag, np=1, tolerance=1e-6):
    """Check every supplied parameter against both adjacent FD step sizes."""
    costs, gradients = parse_summary(run(text, tag, np=np))
    require(len(gradients) == len(costs)*len(parameters), 'wrong gradient count')
    for column, (pattern, occurrence, value) in enumerate(parameters, 1):
        require(any(abs(gradients[cost, column]) > 1e-14 for cost in costs),
                f'{tag}: parameter {column} has only zero gradients')
        for relative in (1e-4, 1e-5):
            h = abs(value)*relative
            plus = replace_nth(text, pattern, f'{value+h:.17e}', occurrence)
            minus = replace_nth(text, pattern, f'{value-h:.17e}', occurrence)
            jp = parse_summary(run(plus, tag, np=np))[0]
            jm = parse_summary(run(minus, tag, np=np))[0]
            for cost in costs:
                close((jp[cost]-jm[cost])/(2*h), gradients[cost, column],
                      rtol=tolerance, atol=1e-11)
    cleanup_run(tag)
    print(f'PASS {tag}: {len(parameters)} parameters, both FD step sizes ({np} ranks)', flush=True)
    return costs, gradients


RATE = r'^    RATE_CONSTANT (\S+)'
HALF = r'^      HALF_SATURATION_CONSTANT (\S+)'


def multireaction(np):
    text = example('1d4')
    start = text.index('  MICROBIAL_REACTION\n')
    end = text.index('  IMMOBILE_DECAY_REACTION\n', start)
    second = text[start:end].replace('4.d-2', '1.d-2').replace('2.d-4', '4.d-4')
    text = text[:end] + second + text[end:]
    text = with_parameters(text, [parameter('RATE_CONSTANT', 1),
                                 parameter('RATE_CONSTANT', 2),
                                 parameter('HALF_SATURATION_CONSTANT', 1, 'O2(aq)'),
                                 parameter('HALF_SATURATION_CONSTANT', 2, 'O2(aq)')])
    parameters = [(RATE, 0, .04), (RATE, 1, .01),
                  (HALF, 0, 2e-4), (HALF, 2, 4e-4)]
    gradient_check(text, parameters, 'multireaction', np)
    if np == 1:
        for units in ('MOLALITY', 'MOLARITY', 'ACTIVITY'):
            units_text = text.replace('  MICROBIAL_REACTION\n',
                                     '  MICROBIAL_REACTION\n'
                                     f'    CONCENTRATION_UNITS {units}\n')
            units_text = log_formulation(units_text)
            units_text = units_text.replace('REFERENCE_LIQUID_DENSITY 1000.d0',
                                            'REFERENCE_LIQUID_DENSITY 1200.d0')
            gradient_check(units_text, parameters, 'multireaction_log_' + units)
        run(replace_nth(text, RATE, '-1.d0'), 'extension_rejections',
            expected_error='RATE_CONSTANT adjoint requires a finite, nonnegative')
        run(replace_nth(text, HALF, '0.d0'), 'extension_rejections',
            expected_error='HALF_SATURATION_CONSTANT adjoint requires')
        invalid = text.replace('SPECIES_NAME O2(aq)\n  /\n',
                               'SPECIES_NAME CO2(aq)\n  /\n')
        run(invalid, 'extension_rejections', expected_error='does not match any MONOD')
        cleanup_run('extension_rejections')



def inhibition(np):
    threshold = r'^      THRESHOLD_CONCENTRATION (\S+)'
    width = r'^      SMOOTHSTEP_INTERVAL (\S+)'
    for kind in ('MONOD', 'SMOOTHSTEP'):
        for direction in ('BELOW', 'ABOVE'):
            block = ('    INHIBITION\n      SPECIES_NAME O2(aq)\n'
                     f'      TYPE {kind}\n'
                     '      THRESHOLD_CONCENTRATION 1.d-3\n'
                     f'      INHIBIT_{direction}_THRESHOLD\n')
            parameters = [parameter('INHIBITION_CONSTANT', 1, 'O2(aq)')]
            perturbations = [(threshold, 0, 1e-3)]
            if kind == 'SMOOTHSTEP':
                block += '      SMOOTHSTEP_INTERVAL 3.d0\n'
                parameters.append(parameter('SMOOTHSTEP_INTERVAL', 1, 'O2(aq)'))
                perturbations.append((width, 0, 3.0))
            block += '    /\n'
            text = example('1d4').replace('    BIOMASS\n', block + '    BIOMASS\n')
            text = with_parameters(text, parameters)
            tag = f'inhibition_{kind}_{direction}'
            reference = gradient_check(text, perturbations, tag, np)
            if np == 1:
                log = log_formulation(text)
                gradient_check(log, perturbations, tag + '_log')
                debug = parse_summary(run(text, tag, args=('-adjoint_debug_jacobianp',)))
                for key in reference[1]:
                    close(reference[1][key], debug[1][key], rtol=1e-10)
                cleanup_run(tag)
                if kind == 'SMOOTHSTEP':
                    for flat_threshold in ('1.d-10', '1.d3'):
                        flat = replace_nth(text, threshold, flat_threshold)
                        _, gradients = parse_summary(run(flat, tag))
                        require(all(abs(g) < 1e-14 for g in gradients.values()),
                                'flat SMOOTHSTEP region has nonzero derivatives')
                    cleanup_run(tag)
    if np == 1:
        run(replace_nth(text, width, '0.d0'), 'inhibition_rejections',
            expected_error='SMOOTHSTEP adjoint requires a finite')
        run(replace_nth(text, threshold, '0.d0'), 'inhibition_rejections',
            expected_error='INHIBITION_CONSTANT adjoint requires a finite')
        run(text.replace('TYPE SMOOTHSTEP', 'TYPE THRESHOLD'),
            'inhibition_rejections', expected_error='supports TYPE MONOD and SMOOTHSTEP')
        cleanup_run('inhibition_rejections')


def yield_decay(np):
    yield_pattern = r'^      YIELD (\S+)'
    base = with_parameters(example('1d4'), [parameter('YIELD'),
                                           parameter('DECAY_RATE_CONSTANT')])
    perturbations = [(yield_pattern, 0, 1e-4), (RATE, 1, 1e-6)]
    gradient_check(base, perturbations, 'yield_decay', np)
    # Distinct biomass yields in two reactions must populate distinct columns.
    start = base.index('  MICROBIAL_REACTION\n')
    end = base.index('  IMMOBILE_DECAY_REACTION\n', start)
    second = base[start:end].replace('4.d-2', '1.d-2').replace('YIELD 1.d-4', 'YIELD 2.d-4')
    multiple = base[:end] + second + base[end:]
    multiple = multiple.replace(parameter('YIELD'),
                                parameter('YIELD') + parameter('YIELD', 2), 1)
    gradient_check(multiple, [(yield_pattern, 0, 1e-4),
                              (yield_pattern, 1, 2e-4), (RATE, 2, 1e-6)],
                   'multireaction_yield_decay', np)
    if np == 1:
        log = log_formulation(base)
        gradient_check(log, perturbations, 'yield_decay_log')
        for column, (pattern, occurrence, value) in enumerate(perturbations, 1):
            # zero is physical: a second-order one-sided difference keeps
            # the parameter nonnegative
            zero = replace_nth(base, pattern, '0.d0', occurrence)
            costs, gradients = parse_summary(run(zero, 'zero_parameters'))
            require(abs(gradients[1, column]) > 1e-10, 'zero parameter lost its derivative')
            for relative in (1e-4, 1e-5):
                h = value*relative
                jp = parse_summary(run(replace_nth(zero, pattern, f'{h:.17e}', occurrence),
                                       'zero_parameters'))[0]
                jpp = parse_summary(run(replace_nth(zero, pattern, f'{2*h:.17e}', occurrence),
                                        'zero_parameters'))[0]
                for cost in costs:
                    fd = (-3*costs[cost]+4*jp[cost]-jpp[cost])/(2*h)
                    close(fd, gradients[cost, column], rtol=1e-6, atol=1e-10)
            cleanup_run('zero_parameters')
        no_biomass = re.sub(r'    BIOMASS\n.*?    /\n', '', base, count=1, flags=re.S)
        run(no_biomass, 'yield_rejections', expected_error='YIELD parameter requires')
        run(base.replace('PARAMETER DECAY_RATE_CONSTANT\n    REACTION_INDEX 1',
                         'PARAMETER DECAY_RATE_CONSTANT\n    REACTION_INDEX 2'),
            'yield_rejections', expected_error='REACTION_INDEX out of range')
        cleanup_run('yield_rejections')
        print('PASS zero yield and decay derivatives, plus invalid selections', flush=True)


def observation_block(*entries):
    """OBSERVATIONS of TIME_SERIES (cell, variable, species, error, file)."""
    text = '  OBSERVATIONS\n'
    for cell, variable, species, error, path in entries:
        text += (f'    TIME_SERIES\n      CELL_ID {cell}\n'
                 f'      OBSERVED_VARIABLE {variable} {species}\n'
                 f'      {error}\n      FILE {path}\n    /\n')
    return text + '  /\n'


def with_observations(text, block):
    text, count = re.subn(r'(?m)^ADJOINT\n', 'ADJOINT\n' + block, text, count=1)
    require(count == 1, 'no ADJOINT block')
    return text


def write_series(path, rows, time_units='h', data_units='mol/kg'):
    path.write_text(f'TIME_UNITS {time_units}\nDATA_UNITS {data_units}\n' +
                    ''.join(f'{t:.17g} {v:.17e}\n' for t, v in rows))


def observed(output):
    """(time [s], simulated, sigma, residual) of each row of the
    observations file."""
    rows = []
    path = Path(output).with_name('case-observations.dat')
    for line in path.read_text().splitlines():
        if not line.startswith('#'):
            fields = line.split()
            rows.append(tuple(float(fields[k]) for k in (1, 6, 7, 8)))
    return rows


def run_to(text, final, times):
    """text run to final [h], with waypoints at times [h] before it."""
    text, count = re.subn(r'(?m)^(\s*FINAL_TIME\s+)\S+\s+\S+',
                          lambda m: m.group(1) + f'{final:.17g} h', text, count=1)
    require(count == 1, 'no FINAL_TIME')
    earlier = [t for t in times if t < final]
    if earlier:
        text, count = re.subn(r'(?m)^OUTPUT\n', 'OUTPUT\n  TIMES h ' +
                              ' '.join(f'{t:.17g}' for t in earlier) + '\n', text, count=1)
        require(count == 1, 'no OUTPUT block')
    return text


def observations(np):
    """OBSERVATIONS: at the final time the objective is the terminal
    concentration's 1/2((C-d)/sigma)^2, also with LOG_FORMULATION and in the
    initial-state sensitivities; a twin over time series at three cells
    (free ion and immobile) matches finite
    differences and the sum of terminal costs at each time (one sweep for
    the M runs it replaces); the sweep in segments, the memory trajectory
    and reruns reproduce it bit for bit; units, a COORDINATE and bad input
    (also a COORDINATE between processes) are checked."""
    base = example('1d4')
    terminal = re.search(r'  COST_FUNCTION TERMINAL_SPECIES_CONCENTRATION\n.*?  /\n',
                         base, re.S).group(0)
    with tempfile.TemporaryDirectory(prefix='pflotran-adjoint-observations-') as directory:
        directory = Path(directory)
        final = directory / 'final.dat'
        write_series(final, [(5, 5e-5)])
        one = observation_block((2, 'FREE_ION', 'CO2(aq)', 'STANDARD_DEVIATION 1.d-4', final))
        deck = with_observations(base, one)
        reference = None
        for log in (False, True):
            text = log_formulation(deck) if log else deck
            costs, gradients = gradient_check(text, [(RATE, 0, .04)], 'observations', np)
            factor = (costs[1]-5e-5)/1e-8
            close(costs[2], 0.5*((costs[1]-5e-5)/1e-4)**2, rtol=1e-10)
            close(gradients[2, 1], factor*gradients[1, 1], rtol=1e-8)
            output = run(text, 'observations_output', np=np, full_field=True)
            first = read_sensitivities(output, 4, cost=1)
            second = read_sensitivities(output, 4, cost=2)
            require(first.keys() == second.keys(), 'sensitivity layout differs')
            for key in first:
                close(second[key], factor*first[key], rtol=1e-8)
            cleanup_run('observations_output')
            if reference:
                for cost in costs:
                    close(costs[cost], reference[0][cost], rtol=1e-7)
                    close(gradients[cost, 1], reference[1][cost, 1], rtol=1e-7)
            reference = costs, gradients
        plain = parse_summary(run(deck, 'observations_forms', np=np))
        # the same observation in sec and mol/g (sigma in the data's units),
        # and as a MEASUREMENT
        write_series(final, [(18000, 5e-8)], 'sec', 'mol/g')
        in_grams = deck.replace('STANDARD_DEVIATION 1.d-4', 'STANDARD_DEVIATION 1.d-7')
        close(parse_summary(run(in_grams, 'observations_forms', np=np))[0][2], plain[0][2],
              rtol=1e-12)
        write_series(final, [(5, 5e-5)])
        measurement = deck.replace(one, '  OBSERVATIONS\n    MEASUREMENT\n      TIME 5 h\n'
                                   '      CELL_ID 2\n      OBSERVED_VARIABLE FREE_ION CO2(aq)\n'
                                   '      VALUE 5.d-5 mol/kg\n      STANDARD_DEVIATION 1.d-4\n'
                                   '    /\n  /\n')
        require(parse_summary(run(measurement, 'observations_forms', np=np)) == plain,
                'a MEASUREMENT differs from the same TIME_SERIES row')
        # the file's keywords in any case, and the cell at a COORDINATE
        final.write_text('time_units h\ndata_units mol/kg\n5 5.d-5\n')
        require(parse_summary(run(deck, 'observations_forms', np=np)) == plain,
                'lowercase TIME_UNITS and DATA_UNITS change the objective')
        write_series(final, [(5, 5e-5)])
        at_point = deck.replace('      CELL_ID 2\n', '      COORDINATE 1.5d0 0.5d0 0.5d0\n')
        require(parse_summary(run(at_point, 'observations_forms', np=np)) == plain,
                'a COORDINATE in cell 2 differs from CELL_ID 2')
        cleanup_run('observations_forms')

        # a twin: data from a run at the true rate constant, fitted at 0.06
        series = {
            'co2': (2, 'FREE_ION', 'CO2(aq)', 'STANDARD_DEVIATION 1.d-6', [1, 2, 2.5, 4, 5]),
            'o2': (3, 'FREE_ION', 'O2(aq)', 'ERROR ABSOLUTE 1.d-6 RELATIVE 1.d-2', [1.5, 3]),
            'biomass': (4, 'IMMOBILE', 'Biomass', 'ERROR ABSOLUTE 1.d-6 RELATIVE 1.d-2', [3, 5]),
        }
        units = {'co2': 'mol/kg', 'o2': 'mol/kg', 'biomass': 'mol/m^3'}
        for name, (_, _, _, _, times) in series.items():
            write_series(directory / f'{name}.dat', [(t, 1e-4) for t in times],
                         data_units=units[name])
        block = observation_block(*[(cell, variable, species, error, directory / f'{name}.dat')
                                    for name, (cell, variable, species, error, _)
                                    in series.items()])
        twin = with_observations(base, block)
        rows = iter(observed(run(twin, 'observations_twin', np=np)))
        for name, (_, _, _, _, times) in series.items():
            write_series(directory / f'{name}.dat', [(t, next(rows)[1]) for t in times],
                         data_units=units[name])
        require(parse_summary(run(twin, 'observations_twin', np=np))[0][2] == 0.,
                'the twin does not fit its own data')
        fitted = replace_nth(twin, RATE, '6.d-2')
        costs, gradients = gradient_check(fitted, [(RATE, 0, .06)], 'observations_twin', np)
        require(costs[2] > 0., 'the fitted twin has no misfit')
        whole = (costs, gradients)
        for args in (('-adjoint_debug_split',), ('-ts_trajectory_type', 'memory')):
            require(parse_summary(run(fitted, 'observations_twin', np=np, args=args)) == whole,
                    f'{args} changes the objective or its gradient')
        parse_rerun_check(run(fitted, 'observations_twin', np=np,
                              args=('-ts_model_rerun_check',)))
        # one sweep = the sum over observations of r/sigma times the terminal
        # concentration's gradient in a run to its time
        output = run(fitted, 'observations_twin', np=np)
        all_times = sorted({t for *_, times in series.values() for t in times})
        total = 0.
        rows = iter(observed(output))
        for cell, variable, species, _, times in series.values():
            for t in times:
                _, _, sigma, r = next(rows)
                cost = (f'  COST_FUNCTION TERMINAL_SPECIES_CONCENTRATION\n'
                        f'    SPECIES_NAME {species}\n    CELL_IDS {cell}\n  /\n')
                short = run_to(replace_nth(base, RATE, '6.d-2').replace(terminal, cost),
                               t, all_times)
                total += r/sigma*parse_summary(run(short, 'observations_sum', np=np))[1][1, 1]
        close(total, gradients[2, 1], rtol=1e-10)
        cleanup_run('observations_twin')
        cleanup_run('observations_sum')

        if np > 1:
            # a COORDINATE on the face between the cells of two processes,
            # of an unstructured grid (a structured one shifts it into one)
            grid = directory / 'line.ugi'
            grid.write_text('4 20\n' + ''.join(
                f'H {i+1} {i+2} {i+7} {i+6} {i+11} {i+12} {i+17} {i+16}\n'
                for i in range(4)) + ''.join(
                f'{x}. {y}. {z}.\n' for z in (0, 1) for y in (0, 1) for x in range(5)))
            unstructured, count = re.subn(r'(?ms)^GRID\n.*?^END\n',
                                          f'GRID\n  TYPE UNSTRUCTURED {grid}\nEND\n',
                                          deck, count=1)
            require(count == 1, 'no GRID block')
            run(unstructured.replace('      CELL_ID 2\n', '      COORDINATE 2.d0 0.5d0 0.5d0\n'),
                'observations_face', np=np, expected_error='is in cells of 2 processes')
            cleanup_run('observations_face')

        if np == 1:
            write_series(final, [(5, 5e-5)])
            for content, error in (
                    ('TIME_UNITS h\nDATA_UNITS M\n5 5.d-5\n', 'Units'),
                    ('DATA_UNITS mol/kg\n5 5.d-5\n', 'requires TIME_UNITS'),
                    ('TIME_UNITS h\n5 5.d-5\n', 'requires DATA_UNITS'),
                    ('TIME_UNITS h\nDATA_UNITS mol/kg\n', 'has no data'),
                    ('TIME_UNITS h\nDATA_UNITS mol/kg\n5 5.d-5\n4 5.d-5\n',
                     'times must increase'),
                    ('TIME_UNITS h\nDATA_UNITS mol/kg\n0 5.d-5\n', 'is not after the start'),
                    ('TIME_UNITS h\nDATA_UNITS mol/kg\n6 5.d-5\n', 'is after the final time'),
                    ('TIME_UNITS h\nDATA_UNITS mol/kg\n5 5.d-5 1.d-6\n',
                     'a row is a time and one value'),
                    ('TIME_UNITS h\nDATA_UNITS mol/kg\n5 inf\n', 'is not finite')):
                (directory / 'bad.dat').write_text(content)
                bad = deck.replace(str(final), str(directory / 'bad.dat'))
                run(bad, 'observations_rejections', expected_error=error)
            for old, new, error in (
                    ('      CELL_ID 2\n', '', 'requires CELL_ID or COORDINATE'),
                    ('      CELL_ID 2\n', '      REGION cell2\n', 'is at one cell'),
                    ('      CELL_ID 2\n', '      CELL_ID 9\n', 'not a cell of the grid'),
                    ('FREE_ION CO2(aq)', 'TOTAL Biomass', 'TOTAL "Biomass" is not a primary'),
                    ('FREE_ION CO2(aq)', 'FREE_ION Biomass', 'is not a primary species'),
                    ('      OBSERVED_VARIABLE FREE_ION CO2(aq)\n', '',
                     'requires OBSERVED_VARIABLE'),
                    ('      STANDARD_DEVIATION 1.d-4\n', '',
                     'requires ERROR or STANDARD_DEVIATION'),
                    ('      STANDARD_DEVIATION 1.d-4\n',
                     '      STANDARD_DEVIATION 1.d-4\n      ERROR ABSOLUTE 0 RELATIVE 0\n',
                     'gives the error twice'),
                    ('STANDARD_DEVIATION 1.d-4', 'STANDARD_DEVIATION Infinity',
                     'STANDARD_DEVIATION must be positive and finite'),
                    ('STANDARD_DEVIATION 1.d-4', 'STANDARD_DEVIATION -1.d-4',
                     'STANDARD_DEVIATION must be positive and finite'),
                    ('STANDARD_DEVIATION 1.d-4', 'ERROR ABSOLUTE Inf RELATIVE 0',
                     'ERROR must be finite and not negative'),
                    (f'      FILE {final}\n', '', 'requires FILE'),
                    (str(final), str(final) + '-missing', 'not found'),
                    ('  OBSERVATIONS\n', '  OBSERVATIONS\n  /\n  OBSERVATIONS\n',
                     'has no TIME_SERIES or MEASUREMENT'),
                    (one, one + one, 'takes one OBSERVATIONS block')):
                require(old in deck, f'missing {old!r}')
                run(deck.replace(old, new, 1), 'observations_rejections', expected_error=error)
            zero = deck.replace('STANDARD_DEVIATION 1.d-4', 'ERROR ABSOLUTE 0 RELATIVE 1.d-2')
            write_series(final, [(5, 0.)])
            run(zero, 'observations_rejections', expected_error='is zero')
            write_series(final, [(5, 5e-5)])
            for old, new, error in (
                    ('      VALUE 5.d-5 mol/kg\n', '      VALUE 5.d-5\n', 'VALUE units'),
                    ('      TIME 5 h\n', '', 'requires TIME'),
                    ('      TIME 5 h\n', '      TIME 5\n', 'TIME units'),
                    ('FREE_ION CO2(aq)', 'IMMOBILE CO2(aq)', 'Units of [volume] were expected')):
                run(measurement.replace(old, new), 'observations_rejections',
                    expected_error=error)
            run(measurement.replace('FREE_ION CO2(aq)', 'IMMOBILE CO2(aq)').replace(
                'mol/kg', 'mol/m^3'), 'observations_rejections',
                expected_error='is not an immobile species')
            run(measurement.replace('FREE_ION CO2(aq)', 'PH').replace(
                'VALUE 5.d-5 mol/kg', 'VALUE 7.d0'), 'observations_rejections',
                expected_error='PH requires H+ as a primary species')
            cleanup_run('observations_rejections')
    print(f'PASS observations: final-time formula, time series vs FD and terminal sums, '
          f'segments, units and rejections ({np} ranks)', flush=True)


CALCITE_RATE = r'^      RATE_CONSTANT (\S+) mol/m\^2-sec'


def chemistry_observations(np):
    """TOTAL and PH observations (ar_carbonate: carbonate speciation with
    activity coefficients): the objective's gradient against finite
    differences, also with LOG_FORMULATION; the simulated values against
    PFLOTRAN's Total [mol/L] (over the density) and pH output at the same
    times; recording them leaves the forward run as it is; segments, the
    memory trajectory and reruns reproduce it; without complexes TOTAL is
    FREE_ION; bad input."""
    base = example('carbonate')
    with tempfile.TemporaryDirectory(prefix='pflotran-adjoint-chemistry-') as directory:
        directory = Path(directory)
        dic, ph = directory / 'dic.dat', directory / 'ph.dat'
        write_series(dic, [(36, 1.5e-3), (72, 1.8e-3), (120, 2.e-3)])
        ph.write_text('TIME_UNITS d\n2.5 6.3\n5 6.6\n')
        block = ('  OBSERVATIONS\n'
                 '    TIME_SERIES\n      CELL_ID 1\n      OBSERVED_VARIABLE TOTAL HCO3-\n'
                 f'      ERROR ABSOLUTE 1.d-5 RELATIVE 2.d-2\n      FILE {dic}\n    /\n'
                 '    TIME_SERIES\n      CELL_ID 4\n      OBSERVED_VARIABLE PH\n'
                 f'      STANDARD_DEVIATION 5.d-2\n      FILE {ph}\n    /\n'
                 '    MEASUREMENT\n      TIME 2 d\n      CELL_ID 2\n      OBSERVED_VARIABLE PH\n'
                 '      VALUE 7.d0\n      STANDARD_DEVIATION 5.d-2\n    /\n'
                 '    MEASUREMENT\n      TIME 5 d\n      CELL_ID 3\n'
                 '      OBSERVED_VARIABLE TOTAL Ca++\n      VALUE 4.d-4 mol/kg\n'
                 '      STANDARD_DEVIATION 1.d-5\n    /\n'
                 '  /\n')
        deck = with_observations(base, block)
        reference = None
        for log in (False, True):
            if log and np > 1:
                break
            text = log_formulation(deck) if log else deck
            tag = f'chemistry_observations_{log}'
            costs, gradients = gradient_check(text, [(CALCITE_RATE, 0, 1e-7)], tag, np,
                                              tolerance=1e-7)
            require(costs[2] > 0. and abs(gradients[2, 1]) > 1., 'no misfit gradient')
            if reference:
                for key in gradients:
                    close(gradients[key], reference[1][key], rtol=1e-6)
            reference = costs, gradients

        # the simulated values are PFLOTRAN's (mol/L at 1.2 kg/L for TOTAL)
        points = ''.join(f'OBSERVATION\n  REGION cell{i}\nEND\n\n' for i in (1, 2))
        heavy = deck.replace('FLUID_PROPERTY\n', 'REFERENCE_LIQUID_DENSITY 1200.d0\n\n'
                             'FLUID_PROPERTY\n').replace('STRATA\n', points + 'STRATA\n')
        heavy = heavy.replace('OUTPUT\n  TIME_UNITS d\n',
                              'OUTPUT\n  TIME_UNITS d\n  PERIODIC_OBSERVATION TIMESTEP 1\n')
        heavy = heavy.replace('REGION cell2\n  COORDINATE 1.5d0 0.5d0 0.5d0\nEND\n',
                              'REGION cell2\n  COORDINATE 1.5d0 0.5d0 0.5d0\nEND\n\n'
                              'REGION cell3\n  COORDINATE 2.5d0 0.5d0 0.5d0\nEND\n\n'
                              'REGION cell4\n  COORDINATE 3.5d0 0.5d0 0.5d0\nEND\n')
        heavy = heavy.replace('OBSERVATION\n  REGION cell2\nEND\n\n',
                              'OBSERVATION\n  REGION cell2\nEND\n\n'
                              'OBSERVATION\n  REGION cell3\nEND\n\n'
                              'OBSERVATION\n  REGION cell4\nEND\n\n')
        output = Path(run(heavy, 'chemistry_values', np=np))
        # the observation points' output, one file per process with points
        columns = {}
        for path in sorted(output.parent.glob('case-obs-*.pft')):
            lines = path.read_text().splitlines()
            header = [c.strip().strip('"') for c in lines[0].split(',')]
            table = [[float(v) for v in line.split()] for line in lines[1:]]
            for i, name in enumerate(header[1:], 1):
                columns[name] = [(row[0], row[i]) for row in table]
        compared = 0
        for line in output.with_name('case-observations.dat').read_text().splitlines():
            if line.startswith('#'):
                continue
            fields = line.split()
            time, cell, variable, species, simulated = (float(fields[1])/86400.,
                                                        int(fields[2]), fields[3],
                                                        fields[4], float(fields[6]))
            name = 'pH' if variable == 'PH' else f'Total {species} [M]'
            matches = [c for c in columns if c.startswith(name) and f'cell{cell} ' in c]
            require(len(matches) == 1, f'missing output column {name} cell{cell}')
            at, value = min(columns[matches[0]], key=lambda r: abs(r[0]-time))
            require(abs(at-time) < 1e-9, f'no output at {time} d')
            if variable == 'PH':
                close(simulated, value, rtol=0., atol=2e-6)
            else:
                close(1.2*simulated, value, rtol=2e-6)
            compared += 1
        require(compared == 7, 'not every observation was compared with the output')
        cleanup_run('chemistry_values')

        # TOTAL and PH evaluate the chemistry at the recorded state: the
        # forward run and the terminal cost are those with FREE_ION at the
        # same cells and times (the same waypoints)
        free = block.replace('TOTAL HCO3-', 'FREE_ION HCO3-').replace(
            'OBSERVED_VARIABLE PH\n', 'OBSERVED_VARIABLE FREE_ION H+\n').replace(
            'TOTAL Ca++', 'FREE_ION Ca++').replace('      VALUE 7.d0\n',
                                                   '      VALUE 1.d-7 mol/kg\n')
        ph.write_text('TIME_UNITS d\nDATA_UNITS mol/kg\n2.5 1.d-6\n5 1.d-6\n')
        plain = parse_summary(run(with_observations(base, free), 'chemistry_forms', np=np))
        ph.write_text('TIME_UNITS d\n2.5 6.3\n5 6.6\n')
        whole = parse_summary(run(deck, 'chemistry_forms', np=np))
        require(whole[0][1] == plain[0][1] and whole[1][1, 1] == plain[1][1, 1],
                'recording TOTAL or PH changed the forward run')
        for args in (('-adjoint_debug_split',), ('-ts_trajectory_type', 'memory')):
            require(parse_summary(run(deck, 'chemistry_forms', np=np, args=args)) == whole,
                    f'{args} changes the objective or its gradient')
        parse_rerun_check(run(deck, 'chemistry_forms', np=np,
                              args=('-ts_model_rerun_check',)))
        cleanup_run('chemistry_forms')

        # without aqueous complexes the total is the free ion
        ar = example('1d4')
        final = directory / 'final.dat'
        write_series(final, [(2, 4.e-4), (5, 1.e-3)])
        one = observation_block((2, 'FREE_ION', 'CO2(aq)', 'STANDARD_DEVIATION 1.d-4', final))
        free_ion = parse_summary(run(with_observations(ar, one), 'chemistry_free', np=np))
        total = parse_summary(run(with_observations(ar, one.replace('FREE_ION', 'TOTAL')),
                                  'chemistry_free', np=np))
        close(total[0][2], free_ion[0][2], rtol=1e-13)
        close(total[1][2, 1], free_ion[1][2, 1], rtol=1e-12)
        cleanup_run('chemistry_free')

        if np == 1:
            for old, new, error in (
                    ('      VALUE 7.d0\n', '      VALUE 7.d0 mol/kg\n', 'PH is dimensionless'),
                    ('      VALUE 4.d-4 mol/kg\n', '      VALUE 4.d-4\n', 'VALUE units are missing'),
                    ('OBSERVED_VARIABLE TOTAL Ca++', 'OBSERVED_VARIABLE TOTAL Calcite',
                     'TOTAL "Calcite" is not a primary species')):
                require(old in deck, f'missing {old!r}')
                run(deck.replace(old, new, 1), 'chemistry_rejections', expected_error=error)
            ph.write_text('TIME_UNITS d\nDATA_UNITS mol/kg\n2.5 6.3\n5 6.6\n')
            run(deck, 'chemistry_rejections', expected_error='PH is dimensionless')
            ph.write_text('TIME_UNITS d\n2.5 6.3\n5 6.6\n')
            dic.write_text('TIME_UNITS h\n36 1.5d-3\n')
            run(deck, 'chemistry_rejections', expected_error='requires DATA_UNITS')
            cleanup_run('chemistry_rejections')
    print(f'PASS TOTAL and PH observations: FD, output values, forward run unchanged, '
          f'segments, reruns, free ion and rejections ({np} ranks)', flush=True)


BIOMASS = r'^      Biomass\s+(\S+)'
CH2O = r'^      CH2O\(aq\)\s+(\S+)'
O2 = r'^      O2\(aq\)\s+(\S+)'


def initial_parameter(species, region):
    return (f'  PARAMETER INITIAL_CONCENTRATION\n    SPECIES_NAME {species}\n'
            f'    REGION {region}\n  /\n')


def with_cell2_constraint(text, ch2o):
    """ar_1d4 with cell2 initialized by a copy of the initial constraint
    with CH2O(aq) at ch2o."""
    initial = re.search(r'TRANSPORT_CONDITION initial\n.*?\nEND\n', text, re.S).group(0)
    second = initial.replace('initial', 'second').replace('CH2O(aq)  2.d-3', f'CH2O(aq)  {ch2o}')
    require(second.count(ch2o) == 1, 'no CH2O(aq) in the initial constraint')
    coupler = 'INITIAL_CONDITION\n  TRANSPORT_CONDITION initial\n  REGION all\nEND\n'
    require(coupler in text, 'no initial condition coupler')
    return text.replace(coupler, coupler + '\n' + second + '\nINITIAL_CONDITION\n'
                        '  TRANSPORT_CONDITION second\n  REGION cell2\nEND\n')


def initial_fd(text, gradients, column, pattern, occurrences, value, tag, np, scale=1.):
    """dJ/dp of an INITIAL_CONCENTRATION against central differences of the
    deck values it is (every occurrence), scale = dc0/d(deck value)."""
    for relative in (1e-4, 1e-5):
        h = value*relative
        jp, jm = text, text
        for occurrence in occurrences:
            jp = replace_nth(jp, pattern, f'{value+h:.17e}', occurrence)
            jm = replace_nth(jm, pattern, f'{value-h:.17e}', occurrence)
        plus = parse_summary(run(jp, tag, np=np))[0]
        minus = parse_summary(run(jm, tag, np=np))[0]
        for cost in plus:
            close((plus[cost]-minus[cost])/(2*h)/scale, gradients[cost, column],
                  rtol=1e-6, atol=1e-13)


def initial_concentrations(np):
    """PARAMETER INITIAL_CONCENTRATION: immobile biomass in all cells and
    the free CH2O(aq) of cell2 (its own constraint) against central
    differences of the deck's values, also with LOG_FORMULATION; a total
    at density 1.2 kg/L is 1.2 times the free molality; reruns at 1.1 p0
    equal a run of the deck there, bit for bit; TAO's gradient test (through
    the reruns' initial states); the objective of OBSERVATIONS; bad input."""
    text = with_cell2_constraint(example('1d4'), '3.d-3')
    text = text.replace('  PARAMETER RATE_CONSTANT\n',
                        initial_parameter('Biomass', 'all') +
                        initial_parameter('CH2O(aq)', 'cell2') +
                        '  PARAMETER RATE_CONSTANT\n', 1)
    for log in (False, True):
        if log and np > 1:
            break
        deck = log_formulation(text) if log else text
        tag = f'initial_{log}'
        output = run(deck, tag, np=np)
        costs, gradients = parse_summary(output)
        names = re.findall(r'dJ_1/dp_\d+ \[(\S+)\] =', Path(output).read_text())
        require(names == ['INITIAL_CONCENTRATION[Biomass,all]',
                          'INITIAL_CONCENTRATION[CH2O(aq),cell2]',
                          'MICROBIAL_RATE_CONSTANT[1]'], f'parameters {names}')
        # the biomass is in both constraints, the CH2O(aq) of cell2 second
        initial_fd(deck, gradients, 1, BIOMASS, (0, 1), 1e-4, tag, np)
        initial_fd(deck, gradients, 2, CH2O, (1,), 3e-3, tag, np)
        cleanup_run(tag)

    # the parameter is the free molality: a total at 1.2 kg/L is 1.2 times it
    dense = text.replace('REFERENCE_LIQUID_DENSITY 1000.d0', 'REFERENCE_LIQUID_DENSITY 1200.d0')
    dense = dense.replace('O2(aq)    1.5d-3   F', 'O2(aq)    1.5d-3   T').replace(
        initial_parameter('Biomass', 'all') + initial_parameter('CH2O(aq)', 'cell2'),
        initial_parameter('O2(aq)', 'all'))
    gradients = parse_summary(run(dense, 'initial_dense', np=np))[1]
    initial_fd(dense, gradients, 1, O2, (0, 1), 1.5e-3, 'initial_dense', np, scale=1/1.2)
    cleanup_run('initial_dense')

    # reruns set the initial state: at 1.1 p0 a run of the deck there
    parameters, costs, gradients = parse_rerun_check(
        run(text, 'initial_rerun', np=np, args=('-ts_model_rerun_check',)))
    scaled = replace_nth(replace_nth(text, BIOMASS, f'{parameters[1]:.17e}', 0),
                         BIOMASS, f'{parameters[1]:.17e}', 1)
    scaled = replace_nth(scaled, CH2O, f'{parameters[2]:.17e}', 1)
    scaled = replace_nth(scaled, RATE, f'{parameters[3]:.17e}', 0)
    require(parse_summary(run(scaled, 'initial_rerun', np=np)) == (costs, gradients),
            'a rerun at 1.1 p0 differs from a run of the deck there')
    cleanup_run('initial_rerun')

    with tempfile.TemporaryDirectory(prefix='pflotran-adjoint-initial-') as directory:
        series = Path(directory) / 'co2.dat'
        write_series(series, [(3, 1.e-3), (5, 1.5e-3)])
        terminal = re.search(r'  COST_FUNCTION TERMINAL_SPECIES_CONCENTRATION\n.*?  /\n',
                             text, re.S).group(0)
        fitted = with_observations(text.replace(terminal, ''), observation_block(
            (2, 'FREE_ION', 'CO2(aq)', 'STANDARD_DEVIATION 1.d-5', series)))
        costs, gradients = parse_summary(run(fitted, 'initial_observations', np=np))
        initial_fd(fitted, gradients, 1, BIOMASS, (0, 1), 1e-4, 'initial_observations', np)
        if np == 1:
            # at the start; TAO's first step reaches a flat region
            optimized = fitted.replace('ADJOINT\n', 'ADJOINT\n  OPTIMIZATION\n'
                                       '    MAXIMUM_ITERATIONS 1\n  /\n', 1)
            output = run(optimized, 'initial_observations', args=('-adj_tao_test_gradient',))
            ratios = [float(v) for v in re.findall(
                r'2-norm \|\|G - Gfd\|\|/\|\|G\|\| = (\S+),',
                Path(output).with_suffix('.stdout').read_text())]
            require(ratios and ratios[0] < 1e-6, f'TAO gradient test: {ratios}')
        cleanup_run('initial_observations')

    if np == 1:
        block = initial_parameter('Biomass', 'all')
        for old, new, error in (
                (block, block.replace('    REGION all\n', ''), 'requires SPECIES_NAME and REGION'),
                (block, block.replace('Biomass', 'Calcite'), 'is not a primary or immobile'),
                (block, block.replace('REGION all', 'REGION nowhere'), 'region "nowhere" not found'),
                (block, block.replace('    REGION all\n', '    REGION all\n    KINETICS MICROBIAL\n'),
                 'Keyword "KINETICS" not recognized'),
                (block, block + initial_parameter('Biomass', 'cell2'), 'their regions overlap'),
                (initial_parameter('CH2O(aq)', 'cell2'), initial_parameter('CH2O(aq)', 'all'),
                 'the parameter is one value')):
            require(old in text, f'missing {old!r}')
            run(text.replace(old, new, 1), 'initial_rejections', expected_error=error)
        cleanup_run('initial_rejections')
    print(f'PASS initial concentrations: FD of the deck values, free molality at a density, '
          f'reruns, TAO gradient test, observations and rejections ({np} ranks)', flush=True)


def nonunit_example(volume=2.0):
    text = example('1d4')
    for old, new in (('4.d0 1.d0 1.d0', f'{4*volume} 1.d0 1.d0'),
                     ('1.5d0 0.5d0 0.5d0', f'{1.5*volume} 0.5d0 0.5d0'),
                     ('POROSITY 1.d0', 'POROSITY 0.35d0')):
        require(old in text, f'missing {old}')
        text = text.replace(old, new)
    return text


def biomass_geometry(np):
    for aqueous in (False, True):
        for volume in (0.1, 2.0):
            text = nonunit_example(volume)
            biomass_pattern, biomass_value, yield_value = r'^      Biomass\s+(\S+)', 1e-4, 1e-4
            if aqueous:
                text = text.replace('REFERENCE_LIQUID_DENSITY 1000.d0',
                                    'REFERENCE_LIQUID_DENSITY 1200.d0')
                text = text.replace('    CO2(aq)\n', '    CO2(aq)\n    Bio(aq)\n', 1)
                text = re.sub(r'  IMMOBILE_SPECIES\n.*?  /\n', '', text, count=1, flags=re.S)
                text = re.sub(r'  IMMOBILE_DECAY_REACTION\n.*?  /\n', '', text, count=1, flags=re.S)
                text = re.sub(r'    IMMOBILE\n.*?    /\n', '', text, count=1, flags=re.S)
                database = Path(__file__).with_name('aqueous_biomass.dat').resolve()
                text = re.sub(r'  DATABASE .*', f'  DATABASE {database}', text)
                text = text.replace('SPECIES_NAME Biomass', 'SPECIES_NAME Bio(aq)')
                text = text.replace('YIELD 1.d-4', 'YIELD 1.d-3')
                text = text.replace('      O2(aq)    1.5d-3   F',
                                    '      Bio(aq)   1.d-7 F\n      O2(aq)    1.5d-3   F')
                biomass_pattern, biomass_value, yield_value = r'^      Bio\(aq\)\s+(\S+)', 1e-7, 1e-3
            text = with_parameters(text, [parameter('RATE_CONSTANT'), parameter('YIELD')])
            for log in (False, True):
                deck = log_formulation(text) if log else text
                tag = f'biomass_{aqueous}_{volume}_{log}'
                gradient_check(deck, [(RATE, 0, .04), (r'^      YIELD (\S+)', 0, yield_value)], tag, np)
                out = run(deck, tag, np=np, full_field=True,
                          args=('-ts_model_column_check', '-adjoint_debug_jacobianp'))
                debug_output = Path(out).read_text()
                # the model's column check runs on one process only
                match = re.search(r'column check:\s*(\d+) entries off', debug_output)
                require(np > 1 or (match and int(match.group(1)) == 0),
                        f'{tag}: state Jacobian column check failed or missing')
                errors = re.findall(r'ADJOINT debug dF/dp .*max rel diff =\s*(\S+)',
                                    debug_output)
                require(len(errors) == 2 and all(float(error) < 1e-8 for error in errors),
                        f'{tag}: missing or failing parameter Jacobian checks')
                sensitivities = read_sensitivities(out, 4)
                # Uniform initial perturbation: sum sensitivities over all cells.
                expected = sum(value for (cell, component), value in sensitivities.items()
                               if component == 4)
                if aqueous:
                    # constraints are in mol/L, FULL_FIELD derivatives per
                    # mol/kg (conversion as in ReactionEquilibrateConstraint)
                    expected *= 1000.0/1200.0
                require(abs(expected) > 1e-8, 'initial biomass test is insensitive')
                for relative in (1e-4, 1e-5):
                    h = biomass_value*relative
                    jp = parse_summary(run(replace_nth(deck, biomass_pattern, f'{biomass_value+h:.17e}'), tag, np=np))[0][1]
                    jm = parse_summary(run(replace_nth(deck, biomass_pattern, f'{biomass_value-h:.17e}'), tag, np=np))[0][1]
                    close((jp-jm)/(2*h), expected, rtol=1e-6, atol=1e-10)
                cleanup_run(tag)
    print(f'PASS biomass volume, porosity, state Jacobian and initial sensitivities ({np} ranks)', flush=True)


def selection_and_domain(np):
    base = nonunit_example()
    monod = ('    INHIBITION\n      SPECIES_NAME O2(aq)\n      TYPE MONOD\n'
             '      THRESHOLD_CONCENTRATION 1.d-3\n      INHIBIT_BELOW_THRESHOLD\n    /\n')
    smooth = monod.replace('TYPE MONOD', 'TYPE SMOOTHSTEP').replace(
        '1.d-3', '2.d-3').replace('INHIBIT_BELOW_THRESHOLD',
                                 'INHIBIT_ABOVE_THRESHOLD\n      SMOOTHSTEP_INTERVAL 3.d0')
    for blocks in (monod+smooth, smooth+monod, smooth+smooth):
        deck = base.replace('    BIOMASS\n', blocks+'    BIOMASS\n')
        selected = with_parameters(deck, [parameter('INHIBITION_CONSTANT', species='O2(aq)')])
        run(selected, 'ambiguous_terms', np=np, expected_error='Ambiguous adjoint term', timeout=30)
        selected = with_parameters(deck, [parameter('SMOOTHSTEP_INTERVAL', species='O2(aq)')])
        if blocks == smooth+smooth:
            run(selected, 'ambiguous_terms', np=np, expected_error='Ambiguous adjoint term', timeout=30)
        else:
            gradient_check(selected, [(r'^      SMOOTHSTEP_INTERVAL (\S+)', 0, 3.)], 'mixed_types', np)
    first_monod = re.search(r'    MONOD\n.*?    /\n', base, re.S).group(0)
    duplicate = base.replace(first_monod, first_monod*2, 1)
    duplicate = with_parameters(duplicate, [parameter('HALF_SATURATION_CONSTANT', species='O2(aq)')])
    run(duplicate, 'ambiguous_terms', np=np, expected_error='Ambiguous adjoint term', timeout=30)
    cleanup_run('ambiguous_terms')
    # the positive-denominator violation lives on one rank only
    invalid = base.replace('      HALF_SATURATION_CONSTANT 2.d-4',
                           '      HALF_SATURATION_CONSTANT 2.d-4\n      THRESHOLD_CONCENTRATION 1.d-3')
    invalid = with_parameters(invalid, [parameter('HALF_SATURATION_CONSTANT', species='O2(aq)')])
    # reach parameter assembly before the unphysical rate exhausts oxygen
    # and trips the forward positivity safeguard
    invalid = invalid.replace('FINAL_TIME 5.d0 h', 'FINAL_TIME 0.5d0 h')
    invalid = invalid.replace('TIMESTEP_SIZE 1.d0 h', 'TIMESTEP_SIZE 0.1d0 h')
    initial = re.search(r'TRANSPORT_CONDITION initial\n.*?\nEND', invalid, re.S).group(0)
    second = initial.replace('initial', 'violating').replace('1.5d-3', '1.d-4')
    for bounds in ('0 0 0\n    4 1 1', '4 0 0\n    8 1 1'):
        deck = invalid.replace('END_SUBSURFACE', second + '\nREGION half\n  COORDINATES\n    ' + bounds +
                               '\n  /\nEND\nINITIAL_CONDITION\n  TRANSPORT_CONDITION violating\n'
                               '  REGION half\nEND\nEND_SUBSURFACE')
        run(deck, 'denominator_rejection', np=np, expected_exit=87, timeout=30,
            expected_error='K + concentration - threshold > 0')
    cleanup_run('denominator_rejection')
    print(f'PASS unambiguous selection and collective domain rejection ({np} ranks)', flush=True)


def extra_branches(np):
    base = nonunit_example()
    threshold = base.replace('      HALF_SATURATION_CONSTANT 2.d-4',
                             '      HALF_SATURATION_CONSTANT 2.d-4\n      THRESHOLD_CONCENTRATION 1.d-4')
    threshold = with_parameters(threshold, [parameter('HALF_SATURATION_CONSTANT', species='O2(aq)')])
    gradient_check(threshold, [(HALF, 0, 2e-4)], 'monod_threshold', np)
    arrhenius = base.replace('    RATE_CONSTANT 4.d-2',
                            '    RATE_CONSTANT 4.d-2\n    ACTIVATION_ENERGY 5.d4')
    arrhenius = arrhenius.replace('REFERENCE_LIQUID_DENSITY', 'REFERENCE_TEMPERATURE 35.d0\nREFERENCE_LIQUID_DENSITY')
    require('ACTIVATION_ENERGY' in arrhenius and 'REFERENCE_TEMPERATURE' in arrhenius,
            'arrhenius deck')
    reference = gradient_check(arrhenius, [(RATE, 0, .04)], 'arrhenius', np)
    # the rate constant's units convert the rate constant, also with the
    # activation energy before it (it was scaled instead)
    for text in (arrhenius.replace('    RATE_CONSTANT 4.d-2\n', '    RATE_CONSTANT 3.456d3 1/d\n'),
                 arrhenius.replace('    RATE_CONSTANT 4.d-2\n    ACTIVATION_ENERGY 5.d4\n',
                                   '    ACTIVATION_ENERGY 5.d1 kJ/mol\n'
                                   '    RATE_CONSTANT 3.456d3 1/d\n')):
        require(text != arrhenius, 'units deck unchanged')
        costs, gradients = parse_summary(run(text, 'arrhenius_units', np=np))
        for key in costs:
            close(costs[key], reference[0][key], rtol=1e-13)
        for key in gradients:
            close(gradients[key], reference[1][key], rtol=1e-12)
    cleanup_run('arrhenius_units')
    # a disabled reaction (k = 0, reaction 2) still has a nonzero dJ/dk
    start = arrhenius.index('  MICROBIAL_REACTION\n')
    end = arrhenius.index('  IMMOBILE_DECAY_REACTION\n', start)
    disabled = replace_nth(arrhenius[start:end], RATE, '0.d0')
    zero = arrhenius[:end] + disabled + arrhenius[end:]
    zero = with_parameters(zero, [parameter('RATE_CONSTANT', 2)])
    costs, gradients = parse_summary(run(zero, 'zero_rate', np=np,
                                         args=('-adjoint_debug_jacobianp',)))
    require(abs(gradients[1, 1]) > 1e-8, 'zero k lost its derivative')
    for h in (4e-6, 4e-7):
        jp = parse_summary(run(replace_nth(zero, RATE, f'{h:.17e}', 1), 'zero_rate', np=np))[0][1]
        jpp = parse_summary(run(replace_nth(zero, RATE, f'{2*h:.17e}', 1), 'zero_rate', np=np))[0][1]
        close((-3*costs[1]+4*jp-jpp)/(2*h), gradients[1, 1], rtol=1e-6, atol=1e-11)
    cleanup_run('zero_rate')
    print(f'PASS Monod threshold, Arrhenius and zero rate constant ({np} ranks)', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--np', type=int, default=1)
    parser.add_argument('--executable')
    parser.add_argument('--feature', default='all',
                        choices=('all', 'multireaction', 'inhibition', 'yield_decay',
                                 'observations', 'chemistry_observations',
                                 'initial_concentrations',
                                 'geometry', 'selection', 'branches'))
    args = parser.parse_args()
    if args.executable:
        os.environ['PFLOTRAN_TEST_EXECUTABLE'] = str(Path(args.executable).resolve())
    test_name = os.environ.get('PFLOTRAN_TEST_NAME', '')
    np = 2 if test_name.endswith('-np2') else args.np
    if args.feature in ('all', 'multireaction'):
        multireaction(np)
    if args.feature in ('all', 'inhibition'):
        inhibition(np)
    if args.feature in ('all', 'yield_decay'):
        yield_decay(np)
    if args.feature in ('all', 'observations'):
        observations(np)
    if args.feature in ('all', 'chemistry_observations'):
        chemistry_observations(np)
    if args.feature in ('all', 'initial_concentrations'):
        initial_concentrations(np)
    if args.feature in ('all', 'geometry'):
        biomass_geometry(np)
    if args.feature in ('all', 'selection'):
        selection_and_domain(np)
    if args.feature in ('all', 'branches'):
        extra_branches(np)
    print('PASS RT adjoint extensions', flush=True)


if __name__ == '__main__':
    main()
