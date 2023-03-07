#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GFF SQUISHER!

@author: Oh Purin
@email: ablaglou@cestcertain.org
@project: https://github.com/ab-n
"""
import argparse
import logging
import sys

logging.basicConfig(format="%(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)


def analyze_contig(gff_in):
    """
    """
    logger.info("Analyzing contigs...")
    contig_dict = dict()
    add_to_seq = 0
    for linus in gff_in:
        if linus == "##FASTA\n":
            break
    for linus in gff_in:
        linus = linus.strip('\n')
        if linus.startswith('>'):
            # add cumsum of the contig size
            contig_dict[linus.lstrip('>')] = add_to_seq
        else:
            add_to_seq += len(linus)
    return(contig_dict)


def squish_contigs(contig_dict, gff_in, gff_out):
    """
    """
    logger.info("Squishing time!")
    contig_new_name = "squish"
    gff_in.seek(0)
    # Deal with the header
    for linus in gff_in:
        if linus == "##FASTA\n":
            break
        elif not linus.startswith("##"):
            annot = linus.split('\t')
            curr_contig = annot[0]
            # Shift the annotation positions
            annot[3] = str(int(annot[3]) + contig_dict[curr_contig])
            annot[4] = str(int(annot[4]) + contig_dict[curr_contig])
            # Homogenize the contig names
            annot[0] = contig_new_name
            annot = "\t".join(annot)
            gff_out.write(annot)

    # We're past header & annotation, deal with the sequences
    gff_out.write(f"##FASTA\n>{contig_new_name}\n")
    for linus in gff_in:
        if not linus.startswith(">"):
            gff_out.write(linus)


def main():
    # Get command-line arguments
    parser = argparse.ArgumentParser(
        description="Squishing a GFF like there's no tomorrow")
    parser.add_argument('file', type=argparse.FileType('r'), nargs=1,
                        help="input file")
    parser.add_argument('-o', '--gff_out', default=sys.stdout,
                        type=argparse.FileType('w'), help='comes out')
    args = parser.parse_args()

    if args.gff_out == sys.stdout:
        logger.info("Output file not specified, printing to stdout")

    contig_dict = analyze_contig(args.file[0])
    squish_contigs(contig_dict, args.file[0], args.gff_out)

    args.file[0].close()
    args.gff_out.close()


if __name__ == '__main__':
    main()
