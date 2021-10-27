#!/usr/bin/env python3
"""
Author : Ken Youens-Clark <kyclark@gmail.com>
Date   : 2021-10-20
Purpose: Join sequences/SNPs
"""

import argparse
import re
import sys
from collections import defaultdict
from Bio import SeqIO
from typing import List, NamedTuple, TextIO


class Args(NamedTuple):
    """ Command-line arguments """
    files: List[TextIO]
    outfile: TextIO


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Join sequences/SNPs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('files',
                        help='Input file(s)',
                        metavar='FILE',
                        nargs='+',
                        type=argparse.FileType('rt'))
    
    parser.add_argument('-s',
                        '--samples',
                        help='File containing list of all samples',
                        metavar='Samples_FILE',
                        type=argparse.FileType('rt'))

    parser.add_argument('-o',
                        '--outfile',
                        help='Output file',
                        metavar='FILE',
                        type=argparse.FileType('wt'),
                        default='out.fa')

    args = parser.parse_args()

    if not args.samples:
        parser.error('Please provide list of all sample names with -s for fasta headings.')

    return args


# --------------------------------------------------
def main() -> None:
    """ Make a jazz noise here """

    args = get_args()

    samples = []
    for line in args.samples:
        samples.append(line.rstrip())

    seqs = defaultdict(list)
    
    for fh in args.files:
        seq_ids = []
        for rec in SeqIO.parse(fh, 'fasta'):
            seq_id = re.sub('^uce-\d+_', '', rec.id)
            if seq_id not in samples:
                sys.exit(f'Error: Sample "{seq_id}"" is not in {args.samples.name}.')
            seq_ids.append(seq_id)
            seqs[seq_id].append(str(rec.seq))
        for sample in samples:
            if sample not in seq_ids:
                seqs[sample].append('-')

    for seq_id, recs in sorted(seqs.items()):
        print('>{}'.format(seq_id), file=args.outfile)
        print("".join(recs), file=args.outfile)

    print(f'Done, see output in "{args.outfile.name}".')


# --------------------------------------------------
if __name__ == '__main__':
    main()
