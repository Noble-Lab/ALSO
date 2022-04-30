#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 17:56:44 2022.

@author: kalavattam
"""

from __future__ import division, absolute_import, print_function
import argparse
import gzip
import pysam
from tqdm import tqdm

# import os
# os.chdir("/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/files_bam")


def detect_filetype_from_path(file):
    """
    Determine type of file (bam, sam, txt, txt.gz) based on extension.

    Parameters
    ----------
    file : str
        Infile, including path

    """
    #  Convert to lowercase for simpler comparisons
    file = file.lower()

    #  Parse sam, bam, txt, txt.gz file extensions
    if file.endswith('.sam'):
        return 'sam'
    elif file.endswith('.bam'):
        return 'bam'
    elif file.endswith('.txt'):
        return 'txt'
    elif file.endswith('.txt.gz'):
        return 'txt.gz'
    else:
        return None


def write_alignment(alignment, outfile=None):
    """
    Print alignment to outfile or, if none is supplied, STDOUT .

    Parameters
    ----------
    alignment : pysam.AlignedSegment
        Alignment to be written to bam outfile

    outfile: str
        Outfile

    """
    if outfile is None:
        print(alignment)
    else:
        outfile.write(alignment)


def filter_qname(infile, txt, outfile=None):
    """
    Filter out reads in a sam/bamfile based on QNAME.

    Parameters
    ----------
    infile: QNAME-sorted bam file, including path (str)

    txt: text file containing QNAME(s) to exclude, including path (str)

    outfile: bam outfile, including path (str)
    '''
    """
    # infile = "Disteche_sample_13.dedup.CAST.sort-c.rm.chr19.sort-n.bam"
    # txt = "Disteche_sample_13.dedup.CAST.sort-c.rm.chr19.singleton.txt.gz"
    # outfile = "Disteche_sample_13.dedup.CAST.sort-c.rm.chr19.rm-singleton.bam"

    #  Read in QNAME-sorted bam file
    bam = pysam.AlignmentFile(infile, "rb")

    #  Set up outfile
    ftype = detect_filetype_from_path(outfile)
    if ftype is None:
        outfile = None
    elif ftype == 'sam':
        outfile = pysam.AlignmentFile(outfile, "w", template=bam)
    elif ftype == 'bam':
        outfile = pysam.AlignmentFile(outfile, "wb", template=bam)
    else:
        raise ValueError(
            "Unknown format for `outfile`: {}".format(outfile)
        )

    # Evaluate QNAME file, read in QNAMEs to be removed, then sort QNAMEs to
    # match BAM for efficient processing
    ftype = detect_filetype_from_path(txt)
    if ftype == 'txt':
        ids = list(set(
            [i.rstrip() for i in open(txt, "r").readlines()]
        ))
        ids.sort(key=str.lower, reverse=True)
    elif ftype == 'txt.gz':
        ids = list(set(
            [i.decode().rstrip() for i in gzip.open(txt, "r").readlines()]
        ))
        ids.sort(key=str.lower, reverse=True)
    else:
        raise ValueError(
            "Unknown format for `txt`: {}".format(txt)
        )
    n_ids = len(ids)

    #  Set up progress bar that draws on total alignment counts
    pbar = tqdm()

    # #FIXME Loop will run on non-QNAME-sorted bam files, but this is incorrect
    reads = bam.fetch(until_eof=True)
    read = next(reads)
    while True:
        #  If QNAME is greater than current top of stack
        if read.query_name > ids[-1]:
            # if this is the last ID
            if n_ids == 1:
                #  Write remaining reads to new file, break out of while loop
                outfile.write(read)
                for r in bam:
                    outfile.write(read)
                break
            #  Otherwise, pop that ID and try again
            else:
                ids.pop()
                n_ids -= 1
        #  Skip reads that match the top of the stack
        elif read.query_name == ids[-1]:
            read = next(reads)
        #  If QNAME is less than top of stack, write it and move on
        else:
            outfile.write(read)
            read = next(reads)
        pbar.update()
    outfile.close()


def main():
    """
    #TODO Some description of the script.

    Parameters
    ----------
    #TODO Some description of the parameters.
    """
    ap = argparse.ArgumentParser(
        description="""
        Filter bam infile to exclude reads based on list of QNAMEs (supplied in
        txt file); write to bam outfile or STDOUT.
        """,
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
        "--txt",
        type=str,
        help="text infile containing QNAMEs to be removed, including path"
    )
    ap.add_argument(
        "-o",
        "--bam_out",
        type=str,
        help="bam outfile, including path"
    )

    arguments = ap.parse_args()
    bam_in = arguments.bam_in
    txt = arguments.txt
    bam_out = arguments.bam_out

    filter_qname(bam_in, txt, bam_out)


main()
