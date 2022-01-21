#!/bin/bash -x
#$ -cwd
#$ -j y
#$ -pe serial 10
#$ -l mfree=4G

module load samtools/latest
module load python/2.7.3
# 20191112 julia/latest not defined & julia/1.1.1 gives errors:
# ERROR: LoadError: Failed to precompile Levenshtein [cb002d51-b850-5775-a4f8-1db47be676b8] to /net/noble/vol2/home/gbonora/.julia/compiled/v1.1/Levenshtein/NjaMZ.ji.
module load julia/0.6.0 # julia/1.1.1 # julia/latest
module load pysam/0.8.1
module load coreutils/8.24
# 20191113
# Actually may not need to do this but doesn't hurt
# 20191112
# Must load java to avoid errors:
# subprocess.CalledProcessError: Command 'java -Xmx1G -jar ...'
module load java/8u151
python ./runall.py \
-R "/net/shendure/vol9/seq/NEXTSEQ/180822_NS500488_0676_AH5FCGAFXY" \
-O "/net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouseDiff_Nmasked/readProcessing" \
-P "sciATAC_mouseDiff" \
-G "/net/noble/vol2/home/gbonora/refdata/mm10.Nmasked.Cast129/bowtie2-index/mm10.Cast129_Nmasked" \
-p 10 --keep_intermediates
# --force_overwrite_all
