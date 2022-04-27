#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 17:56:44 2022.

@author: kalavattam
"""

from __future__ import division, absolute_import, print_function
import argparse
import pysam
from tqdm import tqdm


def detect_filetype_from_path(path):
    """
    Based on its extension, determine what type of file a given path is.

    Parameters
    ----------
    path : str
        Path to check
    """
    #  Convert to lowercase for simpler comparisons
    path = path.tolower()
    #  Remove '.gz' at the end of a pathname if it is present
    if path.endswith('.gz'):
        path = path[:-4]
    #  Parse .sam and .bam extensions
    if path.endswith('.sam'):
        return 'sam'
    elif path.endswith('.bam'):
        return 'bam'
    else:
        return None


def samprint(alignment, output=None):
    """
    Print to output file or STDOUT if none is supplied.

    Parameters
    ----------
    alignment : pysam.AlignedSegment
        Alignment to be written
    output : str
        Output file
    """
    if output is None:
        print(alignment)
    else:
        output.write(alignment)


def filter_qname(file_bam, file_qnames, file_out=None):
    """Filter reads in a SAM/BAM file by their query names."""
    #  Read in QNAME-sorted bam file
    bam = pysam.AlignmentFile(file_bam, "rb")
    #  Output bam file
    ftype = detect_filetype_from_path(file_out)
    if ftype is None:
        output = None
    elif ftype == 'SAM':
        output = pysam.AlignmentFile(file_out, "w", template=bam)
    elif ftype == 'BAM':
        output = pysam.AlignmentFile(file_out, "wb", template=bam)
    else:
        raise ValueError(
            "Unknown output file format for `file_out`: {}".format(file_out)
        )
    #  Read in IDs to be removed (list(set(...)) to only keep unique IDs)
    ids = list(set([l.rstrip() for l in open(file_qnames, 'r').readlines()]))

    #  Sort IDs to match BAM for efficient processing (destructive function)
    ids.sort(key=str.lower, reverse=True)
    n_ids = len(ids)

    #  Progress bar using total alignment counts
    pbar = tqdm()
    reads = bam.fetch(until_eof=True)
    read = next(reads)

    #  Variable for ensuring BAM is sorted by query name
    last_q = reads.query_name
    while True:
        if read.query_name < last_q:
            raise ValueError(
                'Alignment file is not sorted. {} < {}'.format(
                    read.query_name, last_q
                )
            )
        #  If read name is greater than current top of stack
        if read.query_name > ids[-1]:
            #  If this is the last ID
            if n_ids == 1:
                #  Write remaining reads to new file, break out of while loop
                samprint(read, output)
                for r in bam:
                    samprint(read, output)
                break
            #  Otherwise pop that ID, try again
            else:
                ids.pop()
                n_ids -= 1
        #  Skip reads that match the top of the stack
        elif read.query_name == ids[-1]:
            read = next(reads)
        #  If read name is less that top of stack, write it and move on
        else:
            samprint(read, output)
            last_q = read.query_name
            read = next(reads)
        pbar.update()
    if output is not None:
        output.close()


def main():
    """
    #TODO Some description of the script.

    Parameters
    ----------
    #TODO Some description of the parameters.
    """
    ap = argparse.ArgumentParser(
        description="#TODO Write a description for the parser.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    ap.add_argument(
        "-i",
        "--bam_in",
        type=str,
        help="QNAME-sorted bam infile to be filtered, including path"
    )
    ap.add_argument(
        "-q",
        "--qname",
        type=str,
        help="Text infile containing QNAMEs to be removed, including path"
    )
    ap.add_argument(
        "-o",
        "--bam_out",
        type=str,
        help="Bam outfile, including path"
    )

    arguments = ap.parse_args()

    filter_qname(
        file_bam=arguments.bam_in,
        file_qnames=arguments.qname,
        file_out=arguments.bam_out
    )
