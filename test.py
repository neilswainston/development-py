'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys
import urllib


def main(args):
    '''main method.'''
    fle = urllib.urlopen(args[0])
    print fle.read()

if __name__ == '__main__':
    main(sys.argv[1:])
