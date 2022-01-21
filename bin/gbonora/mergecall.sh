#!/bin/bash
#$ -l mfree=64G
#$ -cwd
#$ -j y

source /net/noble/vol2/home/gbonora/.bashrc
source /net/noble/vol2/home/gbonora/.bash_profile

# export PATH=/net/noble/vol2/home/gbonora/anaconda2/bin:/net/gs/vol3/software/modules-sw/samtools/1.9/Linux/RHEL6/x86_64/bin/:/net/gs/vol3/software/modules-sw/bowtie2/2.2.3/Linux/RHEL6/x86_64/scripts:/net/gs/vol3/software/modules-sw/bowtie2/2.2.3/Linux/RHEL6/x86_64/:/net/gs/vol3/software/modules-sw/bcl2fastq/2.16/Linux/RHEL6/x86_64/bin:/net/gs/vol3/software/modules-sw/julia/1.1.1/Linux/RHEL6/x86_64/bin:/net/noble/vol2/home/gbonora/arch/Linux/RHEL6/x86_64/bin:/net/noble/vol2/home/gbonora/bin:/net/noble/vol2/home/gbonora/.local/bin:/noble/bin:/net/noble/vol2/home/gbonora/bin/UCSC.utilities.20160601:/net/noble/vol2/home/gbonora/bin/UCSC.utilities:/net/noble/vol2/home/gbonora/arch/Linux/RHEL6/x86_64/bin:/net/noble/vol2/home/gbonora/bin:/net/noble/vol2/home/gbonora/.local/bin:/noble/bin:/net/noble/vol2/home/gbonora/bin/UCSC.utilities.20160601:/net/noble/vol2/home/gbonora/bin/UCSC.utilities:/net/noble/vol2/home/gbonora/arch/Linux/RHEL6/x86_64/bin:/net/noble/vol2/home/gbonora/bin:/net/noble/vol2/home/gbonora/.local/bin:/noble/bin:/net/noble/vol2/home/gbonora/bin/UCSC.utilities.20160601:/net/noble/vol2/home/gbonora/bin/UCSC.utilities:/net/gs/vol3/software/modules-sw/R/3.5.0/Linux/RHEL6/x86_64/bin/:/net/gs/vol3/software/modules-sw/gm4/1.4.16/Linux/RHEL6/x86_64/bin/:/net/gs/vol3/software/modules-sw/hdf5/1.8.4/Linux/RHEL6/x86_64/bin:/net/gs/vol3/software/modules-sw/autoconf/2.69/Linux/RHEL6/x86_64/bin/:/net/gs/vol3/software/modules-sw/python/2.7.3/Linux/RHEL6/x86_64/bin/:/net/gs/vol3/software/modules-sw/gcc/4.9.1/Linux/RHEL6/x86_64/bin/:/net/gs/vol3/software/modules-sw/gmp/5.0.2/Linux/RHEL6/x86_64/bin:/net/gs/vol3/software/modules-sw/aspera/3.6.1/Linux/RHEL6/x86_64/connect/bin:/usr/share/Modules/bin:/opt/uge/bin/lx-amd64:/opt/uge/local:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/usr/lpp/mmfs/bin:/opt/dell/srvadmin/bin:/net/noble/vol2/home/gbonora/bin:/net/noble/vol2/home/gbonora/bin # condaon
# source activate sciHiC
# export PYTHONPATH=""

source /etc/profile.d/modules.sh
module load modules modules-init modules-gs modules-noble

module load python/2.7.3
module load pysam/0.8.4
# 20191112 julia/latest not defined & julia/1.1.1 gives errors:
# ERROR: LoadError: Failed to precompile Levenshtein [cb002d51-b850-5775-a4f8-1db47be676b8] to /net/noble/vol2/home/gbonora/.julia/compiled/v1.1/Levenshtein/NjaMZ.ji.
module load julia/0.6.0 # julia/1.1.1 # julia/latest
module load coreutils/8.24

hostname
printf "\n\nstart: %s\n\n" "$( date )"

# From README.md
# python path/to/mergeall.py
# -B [list of split.q10.sort.bams]
# -O [Path to output folder]
# -P [Prefix you want on output files]
# -C [path to file with each of the valid barcodes allowed in the experiment]
# --keep_intermediates
#
# The first argument of mergeall.py is a list of the full paths to the BAMs created as output of Part 1.
# The second argument is the full path to where you want the output to go.
# The third argument is optional, and gives each output file an experiment prefix.
# The last argument is the path to a text file with all the allowed barcodes for your experiment,
# one barcode per line (see below).


# python ./mergeall.GB.py -B "/net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouseDiff_Nmasked/readProcessing/sciATAC_mouseDiff.split.q10.sort.bam" -O "/net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouseDiff_Nmasked/readProcessing/sciATAC_mouseDiff" -P "sciATAC_mouseDiff" -C "bc_combos.txt" --override_reads_per_cell "100" --force_overwrite_all --keep_intermediates
python ./mergeall.GB.py -B "/net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouseDiff_Nmasked/readProcessing/sciATAC_mouseDiff.split.q10.sort.bam" -O "/net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouseDiff_Nmasked/readProcessing/sciATAC_mouseDiff" -P "sciATAC_mouseDiff" -C "bc_combos.txt" --override_reads_per_cell "100" --keep_intermediates
# -R "MM10_BLACKLIST"  <-- THIS IS THE DEFAULT NOW

printf "\n\nend: %s\n\n" "$( date )"
