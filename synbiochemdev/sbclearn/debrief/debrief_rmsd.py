'''
sbclearn (c) University of Manchester 2017

sbclearn is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import json
import sys
import urllib

import matplotlib.pyplot as plt
import numpy as np


def get_data(project_id):
    '''Gets data.'''
    url = 'http://debrief.synbiochem.co.uk/data/' + project_id
    response = urllib.urlopen(url)
    data = json.loads(response.read())

    data = sorted([(np.mean(mutation['active_site_rmsd']),
                    np.std(mutation['active_site_rmsd']),
                    mutation['active'])
                   for mutation in data['mutations']
                   if 'active_site_rmsd' in mutation])

    return [np.array(arr) for arr in zip(*data)]


def plot(means, stds, activities):
    '''Plots data.'''
    y_pos = np.arange(len(means))
    plt.bar(y_pos, means, align='center', yerr=stds,
            color=['green' if active else 'red' for active in activities])
    plt.xlabel('Variant')
    plt.ylabel('Mean active site rmsd')
    plt.title('Active site rmsd')

    plt.show()


def main(args):
    '''main method.'''
    means, stds, activities = get_data(args[0])
    plot(means, stds, activities)


if __name__ == '__main__':
    main(sys.argv[1:])
