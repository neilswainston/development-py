'''
Created on 20 Oct 2015

@author: neilswainston
'''
import math
import random
import sys
import RBS_Calculator
import RBS_MC_Design

_RBS_CALC = RBS_Calculator.RBS_Calculator('A', [0, 0])


def optimise(sequences, tir_target=None,
             dg_target=None, max_iter=10000):
    '''Simulated annealing method for optimising RBS.'''

    # Initialization
    counter = 0
    accepts = 0
    rejects = 0
    r_t = 0.6

    # Define the energy/cost function based on the dg_target and the other,
    # optional targets
    dg_target, energy, calc = _init(sequences, tir_target, dg_target)

    while energy > 0.25 and counter < max_iter:
        counter += 1
        accepted = False

        rbs_new, energy_new = _mutate(sequences, dg_target)

        if RBS_MC_Design.calc_constraints(rbs_new, calc):
            energy_new = float('inf')

        print 'New energy = ', energy_new

        if energy_new < energy:
            # Accept move immediately
            sequences[1] = rbs_new
            energy = energy_new
            accepted = True
            print 'Move immediately accepted'
        elif math.exp((energy - energy_new) / r_t) > random.random():
            # Accept move based on conditional probability
            sequences[1] = rbs_new
            energy = energy_new
            accepts += 1
            accepted = True
            print 'Move conditionally accepted'
        else:
            # Reject move
            rejects += 1

        if accepted:
            calc.print_dG(calc.infinity)

        # Simulated annealing control
        if accepts + rejects > 50:
            ratio = float(accepts) / float(accepts + rejects)

            if ratio > 0.2:
                # Too many accepts, reduce rt
                r_t /= 2.0
                accepts = 0
                rejects = 0

                print 'Accepting too many conditional moves, reducing ' + \
                    'temperature'

            elif ratio < 0.01:
                # Too many rejects, increase rt
                r_t *= 2.0
                accepts = 0
                rejects = 0

                print 'Rejecting too many conditional moves, increasing ' + \
                    'temperature'

    if tir_target is not None:
        # Convert back to TIR (Translation Initiation Rate)
        return (_RBS_CALC.K * math.exp(-calc.dG_total_list[0] /
                                       _RBS_CALC.RT_eff), sequences[1],
                counter)
    else:
        return (calc.dG_total_list[0], sequences[1], counter)


def _init(sequences, tir_target, dg_target):
    '''Initialises the simulated annealing job.'''
    # Check if dg_total or TIR (translation initiation rate) was specified.
    # If TIR, then convert to dg_total.
    dg_target = _RBS_CALC.RT_eff * \
        (_RBS_CALC.logK - math.log(float(tir_target))) \
        if tir_target is not None else dg_target

    # If rbs_Init is given, use it. Otherwise, randomly choose one that is a
    # decent starting point.
    if sequences[1] is None:
        (sequences[1], calc) = \
            RBS_MC_Design.GetInitialRBS(sequences[0], sequences[2], dg_target)
    else:
        calc = RBS_MC_Design.Run_RBS_Calculator(sequences[0], sequences[2],
                                                sequences[1])

    energy = calc.dG_total_list[0] - dg_target
    calc.print_dG(calc.infinity)

    return dg_target, energy, calc


def _mutate(sequences, dg_target):
    '''Mutates and scores RBS.'''
    weighted_moves = [('insert', 0.1), ('delete', 0.1), ('replace', 0.8)]
    move = RBS_MC_Design.weighted_choice(weighted_moves)
    pos = int(random.random() * len(sequences[1]))

    if move == 'insert':
        letter = random.choice(['A', 'T', 'G', 'C'])
        rbs_new = sequences[1][0:pos] + letter + \
            sequences[1][pos:len(sequences[1])]
    if move == 'delete' and len(sequences[1]) > 1:
        rbs_new = sequences[1][0:pos] + sequences[1][pos + 1:len(sequences[1])]
    if move == 'replace':
        letter = random.choice(['A', 'T', 'G', 'C'])
        rbs_new = sequences[1][0:pos] + letter + \
            sequences[1][pos + 1:len(sequences[1])]

    rbs_new = RBS_MC_Design.RemoveStartCodons(rbs_new)

    if len(rbs_new) > RBS_MC_Design.Max_RBS_Length:
        rbs_new = rbs_new[
            len(rbs_new) - RBS_MC_Design.Max_RBS_Length:len(rbs_new) + 1]

    calc = RBS_MC_Design.Run_RBS_Calculator(
        sequences[0], sequences[2], rbs_new, True)

    return rbs_new, calc.dG_total_list[0] - dg_target


def main(argv):
    '''main method.'''
    print optimise([argv[1], None, argv[2]], tir_target=float(argv[3]))


if __name__ == '__main__':
    main(sys.argv)
