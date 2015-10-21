'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import math
import random


def optimise(solution, max_iter=10000, verbose=False):
    '''Optmises a solution with simulated annealing.'''

    # Initialization:
    counter = 0
    accepts = 0
    rejects = 0
    r_temp = 0.6

    energy = solution.get_energy()

    print str(counter) + '\t' + str(solution) + '\t' + \
        str(energy)

    while energy > 0.25 and counter < max_iter:
        counter += 1
        energy_new = solution.mutate(verbose=verbose)

        if energy_new < energy:
            # Accept move immediately:
            solution.accept()
            energy = energy_new
            print str(counter) + '\t' + str(solution) + '\t' + \
                str(energy)
        elif energy == energy_new:
            # Take no action:
            continue
        elif math.exp((energy - energy_new) / r_temp) > random.random():
            # Accept move based on conditional probability:
            solution.accept()
            energy = energy_new
            print str(counter) + '\t' + str(solution) + '\t' + \
                str(energy)
            accepts += 1
        else:
            # Reject move:
            rejects += 1

        # Simulated annealing control:
        if accepts + rejects > 50:
            if float(accepts) / float(accepts + rejects) > 0.2:
                # Too many accepts, reduce r_temp:
                r_temp /= 2.0
                accepts = 0
                rejects = 0
            elif float(accepts) / float(accepts + rejects) < 0.01:
                # Too many rejects, increase r_temp:
                r_temp *= 2.0
                accepts = 0
                rejects = 0

    return (solution, counter)
