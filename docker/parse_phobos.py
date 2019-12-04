#!/usr/bin/env python3

import argparse
import re
import sys


def reverse_complement(seq):
    '''Reverse complements sequence string.'''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))


def parse_fastq_phobos(filename, minflanking, flt_only):
    '''
    Loops sequentially through PHOBOS output for updating each read sequence
    with it's largest microsatellite. Assumes that the PHOBOS output is
    sequentially ordered by seq-name.

    Prints to stdout an interleaved paired FASTA.
    '''
    current_seq = None
    phobos_repeat = set([])
    with open(filename, 'r') as handle:
        for line in handle:
            # skip header and column definition line
            if line[0] == '#' or re.match("seq-name", line):
                continue
            cols = line.strip('\n').split('\t')
            # skip if flanking sequence is smaller than minimum
            if len(cols[14]) < minflanking or len(cols[15]) < minflanking:
                continue
            try:
                ms_len = int(cols[5])
            except ValueError as e:
                continue
            if flt_only:
                out = line
            else:
                out = ">%saaa%iaaa%i_%s_%s_%s_%s_%i_1\n%s\n>%saaa%iaaa%i_%s_%s_%s_%s_%i_2\n%s\n" % \
                    (
                    cols[0], len(cols[14]), len(cols[15]), cols[2], 
                    cols[3], cols[4], cols[5], len(cols[6]), 
                    cols[14], 
                    cols[0], len(cols[14]), len(cols[15]), cols[2], 
                    cols[3], cols[4], cols[5], len(cols[6]),
                    reverse_complement(cols[15])
                    )
            # initiate current seq if first line
            if not current_seq:
                current_seq = {
                    'name' : cols[0],
                    'length' : ms_len,
                    'out' : out,
                    'repeat' : cols[13]
                    }
            # checks if new seq, where initiates new current_seq
            # and writes old current_seq
            if current_seq['name'] != cols[0]:
                sys.stdout.write(current_seq['out'])
                phobos_repeat.add(current_seq['repeat'])
                current_seq = {
                    'name' : cols[0],
                    'length' : ms_len,
                    'out' : out,
                    'repeat' : cols[13]
                    }
            else:
                # Updates to new current_seq if larger than the current_seq
                if ms_len > current_seq['length']:
                    current_seq = {
                    'name' : cols[0],
                    'length' : ms_len,
                    'out' : out,
                    'repeat' : cols[13]
                    }
    sys.stdout.write(current_seq['out'])


def write_map(phobos_map):
    '''Writes map to file with key as filename.'''
    for key, val in phobos_map.items():
        filename = key + ".flt.phobos"
        with open(filename, 'w') as handle:
            for cols in val:
                out = ">%saaa%iaaa%i_%s_%s_%s_%s_%i\n%s %s\n" % \
                    (cols[0], len(cols[14]), len(cols[15]), cols[2], 
                    cols[3], cols[4], cols[5], len(cols[6]), cols[14], 
                    cols[15])
                handle.write(out)


def main():
    desc = "Parses PHOBOS output for largest microsatellite per read"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--minflanking", metavar="INT", 
                        default=10, type=int,
                        help="minimum flanking sequence size")
    parser.add_argument("--fltonly", action='store_true',
                        help="filters PHOBOS output only for largest MS")
    parser.add_argument("phobos", metavar="PHOBOS",
                        help="PHOBOS file")
    args = parser.parse_args()
    parse_fastq_phobos(args.phobos, args.minflanking, args.fltonly)


if __name__ == "__main__":
    main()
