'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import gzip
import os
import re


def unzip_file(input_filepath, destination):
    '''Unzips files to target destination.'''
    input_filename = os.path.split(input_filepath)[-1:][0]
    target_filename = os.path.join(destination, input_filename[3:-7] + '.pdb')

    input_file = gzip.open(input_filepath, 'rb')
    output_file = open(target_filename, 'wb')

    for line in input_file:
        output_file.write(line)

    input_file.close()
    output_file.close()


def unzip_dir(input_directory, pattern, destination):
    '''Recursively unzips all files in directory matching the pattern.'''
    for root, subdirs, files in os.walk(input_directory):
        for filename in files:
            if re.match(pattern, filename):
                unzip_file(os.path.join(root, filename), destination)

        for subdir in subdirs:
            unzip_dir(os.path.join(input_directory, subdir), pattern,
                      destination)


def main():
    '''main method.'''
    destination = os.path.join(os.path.expanduser('~'), 'structure_utils',
                               'pdb')

    input_directory = os.path.join(os.path.expanduser('~'), 'structure_utils',
                                   'pdb')

    unzip_dir(input_directory, 'pdb.*\\.ent.gz', destination)


if __name__ == '__main__':
    main()
