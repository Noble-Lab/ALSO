#!/usr/bin/env python
#
# 201902019 Giancarlo

import sys
import string
import gzip
import re
import pysam
from subprocess import Popen, PIPE
# import os
# import signal

USAGE = """USAGE: ~.py <*.bam> ...

  This program takes as input a file of filtered sciRNA mapped reads
  with clipped read sequences (and quality scores).

  The SNP file contains the following five columns:
    - chromosome number
    - position (indexed from 1)
    - strand (ignored)
    - reference base
    - SNP base

"""


#############################################################################
# Used to take the complement of a string.
complement = string.maketrans('ATCGN', 'TAGCN')

def reverseComplement(sequence):
    return sequence.upper().translate(complement)[::-1]


#############################################################################
def Illumina18_Qscore(Qcode):
    """
    Converts Illumina1.8 (PhredILL+33) to value
    """
    return ord(Qcode)-33


# #############################################################################
# def os_system(command):
#     """
#     Run system command
#     """
#     process = Popen(command, stdout=PIPE, shell=True)
#     while True:
#         line = process.stdout.readline()
#         if not line:
#             break
#         yield line


#############################################################################
# MAIN
#############################################################################

# Parse the command line.
if len(sys.argv) != 7:
    sys.stderr.write(USAGE)
    sys.exit(1)
sciATACbam = sys.argv[1]
refsnpFileName = sys.argv[2] # 129
altsnpFileName = sys.argv[3] # CAST
root = sys.argv[4]
SNPthresh = int(sys.argv[5])
Qthresh = int(sys.argv[6])

# Read the ref SNP file into a dictionary.
refsnps = {}  # key = (chrom, position), value = (ref allele, alt allele)
refsnpFile = gzip.open(refsnpFileName, "r")
# tally=0
for line in refsnpFile:
    if (line[0] == '#'):
        continue
    words = line.rstrip().split("\t")
    if (len(words) < 5):
        sys.stderr.write("Error reading from %s.\n%s" % (refsnpFileName, line))
        sys.exit(1)
    chrom = words[0]
    position = int(words[1])
    refsnps[(chrom, position)] = (words[3], words[4])
    # tally += 1
    # if tally >= 10000:
    #     break
refsnpFile.close()
sys.stderr.write("Read %d SNPs from %s.\n" % (len(refsnps), refsnpFileName))

# Read the alt SNP file into a dictionary.
altsnps = {}  # key = (chrom, position), value = (ref allele, alt allele)
altsnpFile = gzip.open(altsnpFileName, "r")
# tally=0
for line in altsnpFile:
    if (line[0] == '#'):
        continue
    words = line.rstrip().split("\t")
    if (len(words) < 5):
        sys.stderr.write("Error reading from %s.\n%s" % (altsnpFileName, line))
        sys.exit(1)
    chrom = words[0]
    position = int(words[1])
    altsnps[(chrom, position)] = (words[3], words[4])
    # tally += 1
    # if tally >= 10000:
    #     break
altsnpFile.close()
sys.stderr.write("Read %d SNPs from %s.\n" % (len(altsnps), altsnpFileName))

# Open the BAM files.

# THIS DOESN'T WORK SO USE SYSTEM COMMAND
# inputSAM = pysam.view(sciATACbam)
# sys.stderr.write("Read %d records from %s.\n" % (len(inputSAM), sciATACbam))
# sciATACbam='/net/noble/vol2/home/gbonora/proj/2019_sciATAC_analysis/results/gbonora/20190208_sciATAC_mouseDiff/vijayWork/sciATAC_EB_MKII.split.q10.sort.bam'
# sciATACbam='/net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20190208_sciATAC_mouseDiff/allelicWork/sciATAC_EB_MKII.split.q10.sort.nsorted.bam'
samtoolsCmd = "samtools view " + sciATACbam
process = Popen(samtoolsCmd, stdout=PIPE, shell=True)
# # https://stackoverflow.com/questions/6549669/how-to-kill-process-and-child-processes-from-python/36970370
# # os.kill(process.pid, signal.TERM)
# https://stackoverflow.com/questions/4084322/killing-a-process-created-with-pythons-subprocess-popen
# process.kill()
inputStream = pysam.AlignmentFile(sciATACbam, "rb")
refFile = pysam.AlignmentFile("%s.ref.bam" % root, "wb", template=inputStream)
altFile = pysam.AlignmentFile("%s.alt.bam" % root, "wb", template=inputStream)
ambigFile = pysam.AlignmentFile("%s.ambig.bam" % root, "wb", template=inputStream)
contraFile = pysam.AlignmentFile("%s.contra.bam" % root, "wb", template=inputStream)

# Initialize counters.
lineNumber = 0
readPairTally = 0

numRef = 0
numAlt = 0
numContra = 0
numAmbig = 0


numE1truAltSNP = 0 # Listed Cast SNPs
numE1infAltSNP = 0 # Assumed Cast SNPs
numE1truRefSNP = 0 # Listed 129 SNPs
numE1infRefSNP = 0 # Assumed 129 SNPs
numE1refSNP = 0 # 129 SNPs
numE1altSNP = 0 # Cast SNPs
numE1huhAltSNP = 0 # Base not a listed SNP
numE1huhRefSNP = 0 # Base not a listed SNP
numE1AltNs = 0 # N
numE1RefNs = 0 # N
numE1LowQualAltBases = 0
numE1LowQualRefBases = 0

numE1CoveredBases = 0
numE1AmbigBases = 0


numE2truAltSNP = 0 # Listed Cast SNPs
numE2infAltSNP = 0 # Assumed Cast SNPs
numE2truRefSNP = 0 # Listed 129 SNPs
numE2infRefSNP = 0 # Assumed 129 SNPs
numE2refSNP = 0 # 129 SNPs
numE2altSNP = 0 # Cast SNPs
numE2huhAltSNP = 0 # Base not a listed SNP
numE2huhRefSNP = 0 # Base not a listed SNP
numE2AltNs = 0 # N
numE2RefNs = 0 # N
numE2LowQualAltBases = 0
numE2LowQualRefBases = 0

numE2CoveredBases = 0
numE2AmbigBases = 0


# Traverse the reads file.
# for readSAMLine in process.stdout:
# https://stackoverflow.com/questions/803265/getting-realtime-output-using-subprocess
while True:
    readSAMLine = process.stdout.readline()
    if not readSAMLine: break
    readSAMLine = readSAMLine.strip()
    readLine = inputStream.next()
    lineNumber += 1

    # Skip header lines.
    if (readSAMLine[0] == '@'):
        refFile.write(readSAMLine)
        altFile.write(readSAMLine)
        ambigFile.write(readSAMLine)
        contraFile.write(readSAMLine)
        continue

    # Read in mate pair
    readSAMLinePE = process.stdout.readline()
    readSAMLinePE = readSAMLinePE.strip()
    readLinePE = inputStream.next()
    lineNumber += 1
    readPairTally += 1
    # sys.stderr.write("\n~~~\n%d:\t%s\n"% (lineNumber, readSAMLine))
    # sys.stderr.write("%d:\t%s\n"% (lineNumber, readSAMLinePE))
    # sys.stderr.write("~~~\n%d:\t%s\n"% (lineNumber, readLine))
    # sys.stderr.write("%d:\t%s\n"% (lineNumber, readLinePE))

    # Verify that we have the right number of fields.
    readWords = readSAMLine.split("\t")
    if (len(readWords) < 11):
        sys.stderr.write("Error parsing line %d from %s.\n%s"
                         % (lineNumber, sciATACbam, readSAMLine))
        sys.exit(1)

    # Verify that we have the right number of fields.
    readWordsPE = readSAMLinePE.split("\t")
    if (len(readWords) < 11):
        sys.stderr.write("Error parsing line %d from %s.\n%s"
                         % (lineNumber, sciATACbam, readSAMLine))
        sys.exit(1)

    # Verify that the two lines are for the same read.
    if (readWords[0] != readWordsPE[0]):
        sys.stderr.write("Read mismatch at line %d. readWords != readWordsPE. (%s != %s).\n"
                         % (readPair, readWords[0], readWordsPE[0]))
        sys.exit(1)

    # # For testing
    # if readPairTally > 5:
    #     break

    # sys.stderr.write("%d read pair.\n" % readPairTally)

    # Get the read and its coordinates.
    E1chrom = readWords[2][3:]
    E1startPosition = int(readWords[3])
    E1Read = readWords[9]
    E1Qual = readWords[10]
    E1ReadLength = len(E1Read)
    E1cigar = readWords[5]

    E2chrom = readWordsPE[2][3:]
    E2startPosition = int(readWordsPE[3])
    E2Read = readWordsPE[9]
    E2Qual = readWordsPE[10]
    E2ReadLength = len(E2Read)
    E2cigar = readWordsPE[5]

    pysamE1chrom = readLine.reference_name[3:]
    pysamE1startPosition = readLine.reference_start
    pysamE1Read = readLine.query_sequence
    pysamE1Qual = readLine.query_qualities
    pysamE1ReadLength = readLine.query_length
    pysamE1cigar = readLine.cigarstring

    pysamE2chrom = readLinePE.reference_name[3:]
    pysamE2startPosition = readLinePE.reference_start
    pysamE2Read = readLinePE.query_sequence
    pysamE2Qual = readLinePE.query_qualities
    pysamE2ReadLength = readLinePE.query_length
    pysamE2cigar = readLinePE.cigarstring


    # Parse ref CIGAR 1
    # http://bioinformatics.cvr.ac.uk/blog/tag/cigar-string/
    # E1cigar = '42S108M'; E1Read='AACCTTGTTGTCTATCAACTGATGAATGGAAAGCTGAGGGATAGTCTCTATAGGGTTAGGCATATCCTCTCCCACGAGGCCAAAAAAGGCAGTTCTCTGCTACATGTTCCTAGGGCCTTGGACCAGCCTATGTGTGCTTTTTGATTGGTG'
    # E1cigar = '98M52S'; E1Read='AAAACAACCAGAATGGCATAACGTGCATTTTTTCTCTCTGCTCAAGAGGACTGACTGAGTTGACCCTGTGTGATTTAAATTTGTGATGGAAGAATTTAAGCTGAGGGATAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAA'
    # E1cigar = '1S1D50M'; E1Read='AAAACAACCAGAATGGCATAACGTGCATTTTTTCTCTCTGCTCAAGAGGACTGACTGAGTTGACCCTGTGTGATTTAAATTTGTGATGGAAGAATTTAAGCTGAGGGATAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAA'
    # E1cigar = '1S3I1D50M'; E1Read='AAAACAACCAGAATGGCATAACGTGCATTTTTTCTCTCTGCTCAAGAGGACTGACTGAGTTGACCCTGTGTGATTTAAATTTGTGATGGAAGAATTTAAGCTGAGGGATAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAA'
    E1ReadCIGARsed = ''
    E1QualCIGARsed = ''
    E1startPositionCIGARsed = E1startPosition
    origE1cigar = E1cigar
    prevPos = 0
    if E1cigar != '*':
        firstFlag = True
        while len(E1cigar)>0:
            m = re.search('^(\d+)(\D)', E1cigar)
            mcount = int(m.group(1))
            mtype = m.group(2)
            if mtype in ['M', '=', 'X']:
                E1ReadCIGARsed = E1ReadCIGARsed + E1Read[prevPos:(prevPos+mcount)]
                E1QualCIGARsed = E1QualCIGARsed + E1Qual[prevPos:(prevPos+mcount)]
                prevPos += mcount
            elif mtype in ['S', 'H']:
                prevPos += mcount
                if firstFlag:
                    E1startPositionCIGARsed += mcount
            elif mtype in ['I', 'P']:
                prevPos += mcount
                if firstFlag:
                    E1startPositionCIGARsed += mcount
            elif mtype in ['D', 'N']:
                E1ReadCIGARsed = E1ReadCIGARsed + 'N' * mcount
                E1QualCIGARsed = E1QualCIGARsed + 'N' * mcount
            E1cigar = E1cigar.replace(m.group(0), '', 1)
            firstFlag = False
    else:
        E1ReadCIGARsed = E1Read
        E1QualCIGARsed = E1Qual

    E1ReadCIGARsedLength = len(E1ReadCIGARsed)
    pysamE1ReadCIGARsed = readLine.query_alignment_sequence
    pysamE1ReadCIGARsedLength = readLine.query_alignment_length
    pysamE1startPositionCIGARsed = readLine.query_alignment_start + pysamE1startPosition

    if E1ReadCIGARsedLength != pysamE1ReadCIGARsedLength or E1startPosition != E1startPositionCIGARsed:
        sys.stderr.write("\n~~~\n%s\t%d\t%d\t%s\t%d\t%d\n"% (str(E1chrom), E1startPosition, E1startPositionCIGARsed, origE1cigar, E1ReadLength, E1ReadCIGARsedLength))
        sys.stderr.write("\n%s\t%d\t%d\t%s\t%d\t%d\n"% (str(pysamE1chrom), pysamE1startPosition, pysamE1startPositionCIGARsed, pysamE1cigar, pysamE1ReadLength, pysamE1ReadCIGARsedLength))
        #
        # sys.stderr.write("\n~~~\n%s\t%d\t%d\t%s\t%d\t%d\n%s\n"% (str(E1chrom), E1startPosition, E1startPositionCIGARsed, origE1cigar, E1ReadLength, E1ReadCIGARsedLength, E1ReadCIGARsed))
        # sys.stderr.write("\n%s\t%d\t%d\t%s\t%d\t%d\n%s\n"% (str(pysamE1chrom), pysamE1startPosition, pysamE1startPositionCIGARsed, pysamE1cigar, pysamE1ReadLength, pysamE1ReadCIGARsedLength, pysamE1ReadCIGARsed))
        #
        # sys.stderr.write("\n~~~\n%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\n"% (E1chrom, E1startPosition, E1startPositionCIGARsed, origE1cigar, E1Read, E1ReadCIGARsed, E1ReadLength, E1ReadCIGARsedLength))
        # sys.stderr.write("\n%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\n"% (pysamE1chrom, pysamE1startPosition, pysamE1startPositionCIGARsed, pysamE1cigar, pysamE1Read, pysamE1ReadCIGARsed, pysamE1ReadLength, pysamE1ReadCIGARsedLength))


    # Parse ref CIGAR 2
    # http://bioinformatics.cvr.ac.uk/blog/tag/cigar-string/
    # E2cigar = '42S108M'; E2Read='AACCTTGTTGTCTATCAACTGATGAATGGAAAGCTGAGGGATAGTCTCTATAGGGTTAGGCATATCCTCTCCCACGAGGCCAAAAAAGGCAGTTCTCTGCTACATGTTCCTAGGGCCTTGGACCAGCCTATGTGTGCTTTTTGATTGGTG'
    # E2cigar = '98M52S'; E2Read='AAAACAACCAGAATGGCATAACGTGCATTTTTTCTCTCTGCTCAAGAGGACTGACTGAGTTGACCCTGTGTGATTTAAATTTGTGATGGAAGAATTTAAGCTGAGGGATAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAA'
    # E2cigar = '1S1D50M'; E2Read='AAAACAACCAGAATGGCATAACGTGCATTTTTTCTCTCTGCTCAAGAGGACTGACTGAGTTGACCCTGTGTGATTTAAATTTGTGATGGAAGAATTTAAGCTGAGGGATAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAA'
    # E2cigar = '1S3I1D50M'; E2Read='AAAACAACCAGAATGGCATAACGTGCATTTTTTCTCTCTGCTCAAGAGGACTGACTGAGTTGACCCTGTGTGATTTAAATTTGTGATGGAAGAATTTAAGCTGAGGGATAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAA'
    E2ReadCIGARsed = ''
    E2QualCIGARsed = ''
    E2startPositionCIGARsed = E2startPosition
    origE2cigar = E2cigar
    prevPos = 0
    if E2cigar != '*':
        firstFlag = True
        while len(E2cigar)>0:
            m = re.search('^(\d+)(\D)', E2cigar)
            mcount = int(m.group(1))
            mtype = m.group(2)
            if mtype in ['M', '=', 'X']:
                E2ReadCIGARsed = E2ReadCIGARsed + E2Read[prevPos:(prevPos+mcount)]
                E2QualCIGARsed = E2QualCIGARsed + E2Qual[prevPos:(prevPos+mcount)]
                prevPos += mcount
            elif mtype in ['S', 'H']:
                prevPos += mcount
                if firstFlag:
                    E2startPositionCIGARsed += mcount
            elif mtype in ['I', 'P']:
                prevPos += mcount
                if firstFlag:
                    E2startPositionCIGARsed += mcount
            elif mtype in ['D', 'N']:
                E2ReadCIGARsed = E2ReadCIGARsed + 'N' * mcount
                E2QualCIGARsed = E2QualCIGARsed + 'N' * mcount
            E2cigar = E2cigar.replace(m.group(0), '', 1)
            firstFlag = False
    else:
        E2ReadCIGARsed = E2Read
        E2QualCIGARsed = E2Qual

    E2ReadCIGARsedLength = len(E2ReadCIGARsed)
    pysamE2ReadCIGARsed = readLinePE.query_alignment_sequence
    pysamE2ReadCIGARsedLength = readLinePE.query_alignment_length
    pysamE2startPositionCIGARsed = readLinePE.query_alignment_start + pysamE2startPosition

    if E2ReadCIGARsedLength != pysamE2ReadCIGARsedLength or E2startPosition != E2startPositionCIGARsed:
        sys.stderr.write("\n~~~\n%s\t%d\t%d\t%s\t%d\t%d\n"% (str(E2chrom), E2startPosition, E2startPositionCIGARsed, origE2cigar, E2ReadLength, E2ReadCIGARsedLength))
        sys.stderr.write("\n%s\t%d\t%d\t%s\t%d\t%d\n"% (str(pysamE2chrom), pysamE2startPosition, pysamE2startPositionCIGARsed, pysamE2cigar, pysamE2ReadLength, pysamE2ReadCIGARsedLength))
        #
        # sys.stderr.write("\n~~~\n%s\t%d\t%d\t%s\t%d\t%d\n%s\n"% (str(E2chrom), E2startPosition, E2startPositionCIGARsed, origE2cigar, E2ReadLength, E2ReadCIGARsedLength, E2ReadCIGARsed))
        # sys.stderr.write("\n%s\t%d\t%d\t%s\t%d\t%d\n%s\n"% (str(pysamE2chrom), pysamE2startPosition, pysamE2startPositionCIGARsed, pysamE2cigar, pysamE2ReadLength, pysamE2ReadCIGARsedLength, pysamE2ReadCIGARsed))
        #
        # sys.stderr.write("\n~~~\n%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\n"% (E2chrom, E2startPosition, E2startPositionCIGARsed, origE2cigar, E2Read, E2ReadCIGARsed, E2ReadLength, E2ReadCIGARsedLength))
        # sys.stderr.write("\n%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\n"% (pysamE2chrom, pysamE2startPosition, pysamE2startPositionCIGARsed, pysamE2cigar, pysamE2Read, pysamE2ReadCIGARsed, pysamE2ReadLength, pysamE2ReadCIGARsedLength))


    # check for E1 SNPs
    E1truAltSNP = 0 # Listed Cast SNPs
    E1infAltSNP = 0 # Assumed Cast SNPs
    E1truRefSNP = 0 # Listed 129 SNPs
    E1infRefSNP = 0 # Assumed 129 SNPs
    E1huhAltSNP = 0 # Not a Cast/B6 or 129 SNP
    E1huhRefSNP = 0 # Not a 129/B6 or Cast SNP

    for offset in range(0, E1ReadCIGARsedLength):
        key = (E1chrom, E1startPositionCIGARsed + offset)
        if altsnps.has_key(key):
            numE1CoveredBases += 1
            if Illumina18_Qscore(E1QualCIGARsed[offset]) <= Qthresh:
                numE1LowQualAltBases += 1
            elif E1ReadCIGARsed[offset] == 'N':
                numE1AltNs += 1
            elif altsnps[key][1] == E1ReadCIGARsed[offset]:
                E1truAltSNP += 1 # True CAST SNP
                # sys.stderr.write("E1 alt SNP %s at offset %d, position %d.\n"
                #                  % (E1ReadCIGARsed[offset], offset, E1startPositionCIGARsed + offset))
            # Check to see if a ref SNP occurs at the same position.
            elif refsnps.has_key(key) and (refsnps[key][1] == E1ReadCIGARsed[offset]):
                E1truRefSNP += 1 # True 129 SNP
            # Neither Cast nor 129 SNP
            elif altsnps[key][0] == E1ReadCIGARsed[offset]:
                E1infRefSNP += 1 # Assume B6 base at this locus pertains to 129 too.
            else:
                E1huhAltSNP += 1 # This should not happen as it should either be a Cast or 129/B6 base.
        # If not an alt SNP, check whether its a ref SNP.
        elif refsnps.has_key(key):
            numE1CoveredBases += 1
            if Illumina18_Qscore(E1QualCIGARsed[offset]) <= Qthresh:
                numE1LowQualRefBases += 1
            elif E1ReadCIGARsed[offset] == 'N':
                numE1RefNs += 1
            elif refsnps[key][1] == E1ReadCIGARsed[offset]:
                E1truRefSNP += 1 # True CAST SNP
                # sys.stderr.write("E1 ref SNP %s at offset %d, position %d.\n"
                #                  % (E1ReadCIGARsed[offset], offset, E1startPositionCIGARsed + offset))
            # This would arleady be detected above:
            # # Check to see if a ref SNP occurs at the same position.
            # elif refsnps.has_key(key) and (refsnps[key][1] == E1ReadCIGARsed[offset]):
            #     E1truRefSNP += 1 # True 129 SNP
            # Neither Cast nor 129 SNP
            elif refsnps[key][0] == E1ReadCIGARsed[offset]:
                E1infAltSNP += 1 # Assume B6 base at this locus pertains to Cast too.
            else:
                E1huhRefSNP += 1 # This should not happen as it should either be a 129 or Cast/B6 base.
        else:
            numE1AmbigBases += 1

    E1refSNP = E1truRefSNP + E1infRefSNP
    E1altSNP = E1truAltSNP + E1infAltSNP

    numE1truAltSNP += E1truAltSNP # Listed Cast SNPs
    numE1infAltSNP += E1infAltSNP # Assumed Cast SNPs
    numE1truRefSNP += E1truRefSNP # Listed 129 SNPs
    numE1infRefSNP += E1infRefSNP # Assumed 129 SNPs
    numE1refSNP += E1refSNP # 129 SNPs
    numE1altSNP += E1altSNP # Cast SNPs
    numE1huhAltSNP += E1huhAltSNP
    numE1huhRefSNP += E1huhRefSNP


    # check for E2 SNPs
    E2truAltSNP = 0 # Listed Cast SNPs
    E2infAltSNP = 0 # Assumed Cast SNPs
    E2truRefSNP = 0 # Listed 129 SNPs
    E2infRefSNP = 0 # Assumed 129 SNPs
    E2huhAltSNP = 0 # Not a Cast/B6 or 129 SNP
    E2huhRefSNP = 0 # Not a 129/B6 or Cast SNP

    for offset in range(0, E2ReadCIGARsedLength):
        key = (E2chrom, E2startPositionCIGARsed + offset)
        if altsnps.has_key(key):
            numE2CoveredBases += 1
            if Illumina18_Qscore(E2QualCIGARsed[offset]) <= Qthresh:
                numE2LowQualAltBases += 1
            elif E2ReadCIGARsed[offset] == 'N':
                numE2AltNs += 1
            elif altsnps[key][1] == E2ReadCIGARsed[offset]:
                E2truAltSNP += 1 # True CAST SNP
                # sys.stderr.write("E2 alt SNP %s at offset %d, position %d.\n"
                #                  % (E2ReadCIGARsed[offset], offset, E2startPositionCIGARsed + offset))
            # Check to see if a ref SNP occurs at the same position.
            elif refsnps.has_key(key) and (refsnps[key][1] == E2ReadCIGARsed[offset]):
                E2truRefSNP += 1 # True 129 SNP
            # Neither Cast nor 129 SNP
            elif altsnps[key][0] == E2ReadCIGARsed[offset]:
                E2infRefSNP += 1 # Assume B6 base at this locus pertains to 129 too.
            else:
                E2huhAltSNP += 1 # This should not happen as it should either be a Cast or 129/B6 base.
        # If not an alt SNP, check whether its a ref SNP.
        elif refsnps.has_key(key):
            numE2CoveredBases += 1
            if Illumina18_Qscore(E2QualCIGARsed[offset]) <= Qthresh:
                numE2LowQualRefBases += 1
            elif E2ReadCIGARsed[offset] == 'N':
                numE2RefNs += 1
            elif refsnps[key][1] == E2ReadCIGARsed[offset]:
                E2truRefSNP += 1 # True CAST SNP
                # sys.stderr.write("E2 ref SNP %s at offset %d, position %d.\n"
                #                  % (E2ReadCIGARsed[offset], offset, E2startPositionCIGARsed + offset))
            # This would arleady be detected above:
            # # Check to see if a ref SNP occurs at the same position.
            # elif refsnps.has_key(key) and (refsnps[key][1] == E2ReadCIGARsed[offset]):
            #     E2truRefSNP += 1 # True 129 SNP
            # Neither Cast nor 129 SNP
            elif refsnps[key][0] == E2ReadCIGARsed[offset]:
                E2infAltSNP += 1 # Assume B6 base at this locus pertains to Cast too.
            else:
                E2huhRefSNP += 1 # This should not happen as it should either be a 129 or Cast/B6 base.
        else:
            numE2AmbigBases += 1

    E2refSNP = E2truRefSNP + E2infRefSNP
    E2altSNP = E2truAltSNP + E2infAltSNP

    numE2truAltSNP += E2truAltSNP # Listed Cast SNPs
    numE2infAltSNP += E2infAltSNP # Assumed Cast SNPs
    numE2truRefSNP += E2truRefSNP # Listed 129 SNPs
    numE2infRefSNP += E2infRefSNP # Assumed 129 SNPs
    numE2refSNP += E2refSNP # 129 SNPs
    numE2altSNP += E2altSNP # Cast SNPs
    numE2huhAltSNP += E2huhAltSNP
    numE2huhRefSNP += E2huhRefSNP


    # Decide where to output this read.
    if (E1refSNP >= SNPthresh or E2refSNP >= SNPthresh) and (E1altSNP < SNPthresh and E2altSNP < SNPthresh):
        refFile.write(readLine)
        refFile.write(readLinePE)
        numRef += 1
    elif (E1altSNP >= SNPthresh or E2altSNP >= SNPthresh) and (E1refSNP < SNPthresh and E2refSNP < SNPthresh):
        altFile.write(readLine)
        altFile.write(readLinePE)
        numAlt += 1
    else:
        if E1refSNP < SNPthresh and E1altSNP < SNPthresh and E2refSNP < SNPthresh and E2altSNP < SNPthresh:
            ambigFile.write(readLine)
            ambigFile.write(readLinePE)
            numAmbig += 1
        elif (E1refSNP >= SNPthresh and E1altSNP >= SNPthresh) or (E2refSNP >= SNPthresh and E2altSNP >= SNPthresh):
            # NOTE: No longer includes refalt and altref,
            # Only instances where SNPs contradict one another in same read
            contraFile.write(readLine)
            contraFile.write(readLinePE)
            numContra += 1


inputStream.close()
refFile.close()
altFile.close()
ambigFile.close()
contraFile.close()


# Print a summary.
sys.stdout.write("\n~~~\nSummary:\n")
sys.stdout.write("\n%d total reads.\n" % readPairTally)
sys.stdout.write("%d reads assigned as ref based on one end.\n" % numRef)
sys.stdout.write("%d reads assigned as alt based on one end.\n" % numAlt)
sys.stdout.write("%d reads contradictory or unresolvable.\n" % numContra)
sys.stdout.write("%d reads unassigned because ambiguous.\n" % numAmbig)


sys.stdout.write("\n%d (%.2f%%) SNP-associated bases in E1s.\n" % ( numE1CoveredBases, float(numE1CoveredBases)/(numE1CoveredBases+numE1AmbigBases)*100 ) )
sys.stdout.write("%d (%.2f%%) SNP-unassociated bases in E1s.\n" % ( numE1AmbigBases, float(numE1AmbigBases)/(numE1CoveredBases+numE1AmbigBases)*100 ) )
sys.stdout.write("%d total SNP base calls assigned to 129 allele in E1s.\n" % numE1refSNP)
sys.stdout.write("%d total SNP base calls assigned to Cast in E1s.\n" % numE1altSNP)

totRefSNPsConsideredInE1s = numE1truRefSNP + numE1infAltSNP + numE1huhRefSNP + numE1LowQualRefBases + numE1RefNs
sys.stdout.write("\n%d total # of Sanger 129 SNP base calls considered in E1s.\n" % totRefSNPsConsideredInE1s)
sys.stdout.write("%d 129 base calls at 129 SNP loci in E1s.\n" % numE1truRefSNP)
sys.stdout.write("%d B6 base calls (inferred as Cast-eminating) at 129 SNP loci in E1s.\n" % numE1infAltSNP)
sys.stdout.write("%d non-129/B6/Cast base calls at 129 SNP loci in E1s.\n" % numE1huhRefSNP)
sys.stdout.write("%d Low quality base calls at 129 SNP sites in E1s.\n" % numE1LowQualRefBases)
sys.stdout.write("%d Ns at 129 SNP sites in E1s.\n" % numE1RefNs)

totAltSNPsConsideredInE1s = numE1truAltSNP + numE1infRefSNP + numE1huhAltSNP + numE1LowQualAltBases + numE1AltNs
sys.stdout.write("\n%d total # of Sanger Cast SNPs considered in E1s.\n" % totAltSNPsConsideredInE1s)
sys.stdout.write("%d Cast base calls at Cast SNP loci in E1s.\n" % numE1truAltSNP)
sys.stdout.write("%d B6 base calls (inferred as 129-eminating) at Cast SNP loci in E1s.\n" % numE1infRefSNP)
sys.stdout.write("%d non-129/B6/Cast bases at Cast SNP loci in E1s.\n" % numE1huhAltSNP)
sys.stdout.write("%d Low quality base calls at Cast SNP sites in E1s.\n" % numE1LowQualAltBases)
sys.stdout.write("%d Ns at Cast SNP sites in E1s.\n" % numE1AltNs)


sys.stdout.write("\n%d (%.2f%%) SNP-associated bases in E2s.\n" % ( numE2CoveredBases, float(numE2CoveredBases)/(numE2CoveredBases+numE2AmbigBases)*100 ) )
sys.stdout.write("%d (%.2f%%) SNP-unassociated bases in E2s.\n" % ( numE2AmbigBases, float(numE2AmbigBases)/(numE2CoveredBases+numE2AmbigBases)*100 ) )
sys.stdout.write("%d total SNP base calls assigned to 129 allele in E2s.\n" % numE2refSNP)
sys.stdout.write("%d total SNP base calls assigned to Cast in E2s.\n" % numE2altSNP)

totRefSNPsConsideredInE2s = numE2truRefSNP + numE2infAltSNP + numE2huhRefSNP + numE2LowQualRefBases + numE2RefNs
sys.stdout.write("\n%d total # of Sanger 129 SNP base calls considered in E2s.\n" % totRefSNPsConsideredInE2s)
sys.stdout.write("%d 129 base calls at 129 SNP loci in E2s.\n" % numE2truRefSNP)
sys.stdout.write("%d B6 base calls (inferred as Cast-eminating) at 129 SNP loci in E2s.\n" % numE2infAltSNP)
sys.stdout.write("%d non-129/B6/Cast base calls at 129 SNP loci in E2s.\n" % numE2huhRefSNP)
sys.stdout.write("%d Low quality base calls at 129 SNP sites in E2s.\n" % numE2LowQualRefBases)
sys.stdout.write("%d Ns at 129 SNP sites in E2s.\n" % numE2RefNs)

totAltSNPsConsideredInE2s = numE2truAltSNP + numE2infRefSNP + numE2huhAltSNP + numE2LowQualAltBases + numE2AltNs
sys.stdout.write("\n%d total # of Sanger Cast SNPs considered in E2s.\n" % totAltSNPsConsideredInE2s)
sys.stdout.write("%d Cast base calls at Cast SNP loci in E2s.\n" % numE2truAltSNP)
sys.stdout.write("%d B6 base calls (inferred as 129-eminating) at Cast SNP loci in E2s.\n" % numE2infRefSNP)
sys.stdout.write("%d non-129/B6/Cast bases at Cast SNP loci in E2s.\n" % numE2huhAltSNP)
sys.stdout.write("%d Low quality base calls at Cast SNP sites in E2s.\n" % numE2LowQualAltBases)
sys.stdout.write("%d Ns at Cast SNP sites in E2s.\n" % numE2AltNs)
