#!/bin/bash -x
#$ -cwd
#$ -j y
#$ -pe serial 10
#$ -l mfree=4G

#  Load modules
# module load samtools/latest  # KA, 2021-0928: No longer works
# module load python/2.7.3  # KA, 2021-0928: No longer works
# GB, 2019-1112: julia/latest not defined & julia/1.1.1 gives errors:
# GB, 2019-1112: ERROR: LoadError: Failed to precompile Levenshtein [cb002d51-b850-5775-a4f8-1db47be676b8] to /net/noble/vol2/home/gbonora/.julia/compiled/v1.1/Levenshtein/NjaMZ.ji.
# module load julia/0.6.0 # julia/1.1.1 # julia/latest
# module load pysam/0.8.1
# module load coreutils/8.24
# GB, 2019-1113
# Actually may not need to do this but doesn't hurt
# GB, 2019-1112
# Must load java to avoid errors:
# subprocess.CalledProcessError: Command 'java -Xmx1G -jar ...'
# module load java/8u151

module load samtools/1.9
module load python/2.7.13
module load julia/1.1.1
module load java/1.8.0
#  KA, 2021-0928: modules coreutils and pysam no longer available on the cluster

#  What's available via conda?
#+ - bowtie2
#+ - trimmomatic
#+ - samtools
#+ - julia
#+ - pysam

main="/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data"
directory_bcl="/net/shendure/vol9/seq/NEXTSEQ/180822_NS500488_0676_AH5FCGAFXY"
directory_out="/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021"
directory_genome="/net/noble/vol1/home/kga0/genomes/UCSC.hg38/bowtie2/indices"
prefix="test"
threads=10
python ./runall.py \
-R "${directory_bcl}" \
-O "${main}/${directory_out}" \
-P "${prefix}" \
-G "${directory_genome}" \
-p "${threads}" \
--keep_intermediates

#  Example call from Giancarlo
# python ./runall.py \
# -R "/net/shendure/vol9/seq/NEXTSEQ/180822_NS500488_0676_AH5FCGAFXY" \
# -O "/net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouseDiff_Nmasked/readProcessing" \
# -P "sciATAC_mouseDiff" \
# -G "/net/noble/vol2/home/gbonora/refdata/mm10.Nmasked.Cast129/bowtie2-index/mm10.Cast129_Nmasked" \
# -p 10 \
# --keep_intermediates
# # --force_overwrite_all

#  runall.py help
# python2.7 runall.py -h
# usage: runall.py [-h] -R RUNDIR -O OUTDIR [-P PREFIX] [--miseq] [-E MAXEDIT]
#                  [-G GENOME] [-p NTHREADS] [--force_overwrite_all]
#                  [--keep_intermediates]
#
# A program to convert BCL files to cleaned, corrected and mapped BAM files for
# sci-ATAC-seq analysis.
#
# optional arguments:
#   -h, --help            show this help message and exit
#   -R RUNDIR, --rundir RUNDIR
#                         Run directory containing BCL files
#   -O OUTDIR, --outdir OUTDIR
#                         Output directory
#   -P PREFIX, --prefix PREFIX
#                         Output file prefix, otherwise default(out)
#   --miseq               Run on MiSeq, otherwise assumes NextSeq
#   -E MAXEDIT, --maxedit MAXEDIT
#                         Maximum allowed edit distance (default = 3)
#   -G GENOME, --genome GENOME
#                         Path to genome annotation (default=/net/trapnell/vol1/
#                         genomes/human/hg19/bowtie2/)
#   -p NTHREADS, --nthreads NTHREADS
#                         Number of cores available (default=1)
#   --force_overwrite_all
#                         Force overwrite of all steps of pipeline regardless of
#                         files already present.
#   --keep_intermediates  Skip clean up steps to keep intermediate files.
