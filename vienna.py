'''
Created on 17 May 2017

@author: neilswainston
'''
import subprocess


def write_csv_to_fasta(csv_filename, fasta_filename):
    '''Writes a csv file of name, sequence to a fasta file.'''

    # Open csv file in read mode, fasta file in write mode...
    with open(csv_filename, 'r') as csv_file, open(fasta_filename, 'w') as fasta_file:
        # Read csv file line by line...
        for line in csv_file:
            # Split line into "tokens", separated by commas...
            tokens = line.split(',')

            # Write header to fasta...
            fasta_file.write('>' + tokens[0] + '\n')

            # Write sequence to fasta...
            fasta_file.write(tokens[1])


def run_rnafold(fasta_filename):
    '''Run RNAfold from fasta file.'''
    with open(fasta_filename) as fasta_file:
        process = subprocess.Popen(['RNAfold', '--noPS'], stdin=fasta_file)
        process.wait()


def main():
    '''main method to start the program.'''
    fasta_filename = 'all_seq.fasta.txt'

    # Write fasta file:
    write_csv_to_fasta('all_seq.csv', fasta_filename)

    # Runs the RNAfold program with input from fasta file:
    run_rnafold(fasta_filename)

if __name__ == '__main__':
    main()
