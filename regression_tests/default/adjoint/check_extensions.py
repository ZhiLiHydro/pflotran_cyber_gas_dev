#!/usr/bin/env python3
"""Finite-difference checks of the MODE GIRT_DAE adjoint for several
microbial reactions, Monod, inhibition, yield and decay parameters, and the
least-squares misfit cost. Post-processor of adjoint-extensions and
adjoint-extensions-np2; see check_adjoint.py for direct use."""

import argparse
import os
from pathlib import Path
import re
import tempfile

from check_adjoint import (close, deck, parse_summary,
                           require, cleanup_run, read_sensitivities)
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


def misfit(np):
    # the terminal concentration cost is the reference for J and for
    # dJmisfit/dp = 2*w*(C-obs)/sigma**2 * dC/dp
    base = example('1d4')
    terminal = re.search(r'  COST_FUNCTION TERMINAL_SPECIES_CONCENTRATION\n.*?  /\n',
                         base, re.S).group(0)
    with tempfile.TemporaryDirectory(prefix='pflotran-adjoint-observations-') as directory:
        observations = Path(directory) / 'observations.dat'
        observations.write_text('TIME_UNITS h\nDATA_UNITS mol/kg\n0 0\n5 5.d-5\n')
        misfit_cost = terminal.replace('TERMINAL_SPECIES_CONCENTRATION',
                                       'MISFIT_LEAST_SQUARES').replace(
            '    SPECIES_NAME', f'    OBSERVATION_FILE {observations}\n'
            '    SIGMA 1.d-4\n    WEIGHT 2.d0\n    SPECIES_NAME')
        deck = base.replace(terminal, terminal + misfit_cost)
        reference = None
        for log in (False, True):
            text = log_formulation(deck) if log else deck
            costs, gradients = gradient_check(text, [(RATE, 0, .04)], 'misfit', np)
            factor = 4*(costs[1]-5e-5)/1e-8
            close(costs[2], 2*((costs[1]-5e-5)/1e-4)**2, rtol=1e-10)
            close(gradients[2, 1], factor*gradients[1, 1], rtol=1e-8)
            output = run(text, 'misfit_output', np=np, full_field=True)
            first = read_sensitivities(output, 4, cost=1)
            second = read_sensitivities(output, 4, cost=2)
            require(first.keys() == second.keys(), 'misfit sensitivity layout differs')
            for key in first:
                close(second[key], factor*first[key], rtol=1e-8)
            cleanup_run('misfit_output')
            default = run(text, 'misfit_output', np=np)
            require(not list(Path(default).parent.glob('*-du0.bin')),
                    'misfit wrote FULL_FIELD by default')
            cleanup_run('misfit_output')
            if reference:
                for cost in costs:
                    close(costs[cost], reference[0][cost], rtol=1e-7)
                    close(gradients[cost, 1], reference[1][cost, 1], rtol=1e-7)
            reference = costs, gradients
        # four identical cells, summed across ranks: 4x the one-cell misfit
        all_cells = deck.replace(misfit_cost, misfit_cost.replace('REGION cell2', 'REGION all'))
        aggregate = parse_summary(run(all_cells, 'misfit_cells', np=np))
        close(aggregate[0][2], 4*reference[0][2], rtol=1e-7)
        close(aggregate[1][2, 1], 4*reference[1][2, 1], rtol=1e-7)
        cleanup_run('misfit_cells')
        # the same observation in s and mol/g gives the same cost
        observations.write_text('TIME_UNITS sec\nDATA_UNITS mol/g\n0 0\n18000 5.d-8\n')
        converted = parse_summary(run(deck, 'misfit_units', np=np))
        close(converted[0][2], reference[0][2], rtol=1e-7)
        cleanup_run('misfit_units')
        observations.write_text('TIME_UNITS h\n4 5.d-5\n')
        rejected = Path(run(deck, 'misfit_time_rejection', np=np,
                            expected_error='final observation time must',
                            expected_exit=87, timeout=30))
        require(not list(rejected.parent.glob('TS-data-*')),
                'early observation-time rejection leaked a trajectory')
        require('Finished Initialization' not in rejected.read_text(),
                'observation-time rejection ran the forward solve')
        cleanup_run('misfit_time_rejection')
        if np == 1:
            for content, error in (
                ('TIME_UNITS h\nDATA_UNITS M\n5 5.d-5\n', 'Units'),
                ('TIME_UNITS h\n4 5.d-5\n', 'final observation time must'),
                ('TIME_UNITS h\n5 5.d-5\n4 5.d-5\n', 'observation times must increase'),
                ('TIME_UNITS h\n5 nan\n', 'While reading "column data"'),
                ('TIME_UNITS h\n5 5.d-5 6.d-5\n', 'exactly one data column'),
                ('CYCLIC\nTIME_UNITS h\n0 0\n5 5.d-5\n', 'does not support CYCLIC'),
            ):
                observations.write_text(content)
                run(deck, 'misfit_rejections', expected_error=error)
            observations.write_text('TIME_UNITS h\n5 5.d-5\n')
            run(deck.replace('SIGMA 1.d-4', 'SIGMA 0.d0'), 'misfit_rejections',
                expected_error='SIGMA must be finite and positive')
            run(deck.replace('WEIGHT 2.d0', 'WEIGHT -1.d0'), 'misfit_rejections',
                expected_error='requires finite, nonnegative WEIGHT')
            run(deck.replace(str(observations), str(observations)+'-missing'),
                'misfit_rejections', expected_error='not found.')
            cleanup_run('misfit_rejections')
    print(f'PASS terminal misfit formula, units, LOG and distributed output ({np} ranks)', flush=True)


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
                        choices=('all', 'multireaction', 'inhibition', 'yield_decay', 'misfit',
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
    if args.feature in ('all', 'misfit'):
        misfit(np)
    if args.feature in ('all', 'geometry'):
        biomass_geometry(np)
    if args.feature in ('all', 'selection'):
        selection_and_domain(np)
    if args.feature in ('all', 'branches'):
        extra_branches(np)
    print('PASS RT adjoint extensions', flush=True)


if __name__ == '__main__':
    main()
