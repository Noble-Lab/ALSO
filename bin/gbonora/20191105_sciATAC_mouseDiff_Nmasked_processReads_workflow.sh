##################################################
##################################################
##################################################
#
# 201901111 Nmasked genome
# Started
#
##################################################
##################################################
##################################################
#
# 20190216 started
#
# 20190219 PE reads
#
##################################################

source ~/.bashrc
source ~/.bash_profile

PROJDIR="2019_sciATAC_analysis"
# mkdir -p "${PROJDIR}"
WORKDIR="20191105_sciATAC_mouseDiff_Nmasked"
READDIR="data_20191105_sciATAC_mouseDiff_Nmasked"


##################################################
# Working folder
mkdir -p "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"
cd  "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"










#################################################
#################################################
#################################################
# Set up data folder

VOL8HOME=/net/noble/vol8/gbonora/sciOmics
mkdir -p "${VOL8HOME}"
cd "${VOL8HOME}"
pwd -P

# Create a folder for data
mkdir -p "${VOL8HOME}"/"${PROJDIR}"/data/"${READDIR}"

# Link to PROJDIR
cd "${HOME}"/proj/"${PROJDIR}"
ln -sf "${VOL8HOME}"/"${PROJDIR}"/data .
ls -ltrh data/










# #################################################
# #################################################
# #################################################
# # Step 0. Link to mouse raw data
#
# # Make data folder
# ls -ltrh "${HOME}"/proj/"${PROJDIR}"/data
# mkdir -p "${HOME}"/proj/"${PROJDIR}"/data/"${READDIR}"/rawReads
# cd "${HOME}"/proj/"${PROJDIR}"/data/"${READDIR}"/rawReads
#
# url=http://krishna.gs.washington.edu/content/members/cholilee/sequencing_files/180512_NS500488_0616_AHLYGNBGX5_demultiplexed/
# suffix="_001.fastq.gz"
# for sample in Undetermined_S0; do
#     echo ${sample}
#     for read in I1 R1 R2; do
#       job=wget-${sample}_${read}.job
#       echo wget $url/${sample}_${read}${suffix} > $job
#       qsub -j y -l h_rt=7:59:59 -cwd $job
#     done
# done
# qstat -u gbonora

# SRCDIR=/net/shendure/vol9/seq/NEXTSEQ/180507_NS500488_0611_AH2KHCAFXY
SRCDIR=/net/shendure/vol9/seq/NEXTSEQ/180822_NS500488_0676_AH5FCGAFXY










#################################################
#################################################
#################################################
# Get Hannah pipeline setup
#
# See '20190208_cole-trapnell-lab_ATAC-Pipeline-Dev_githubRepo.sh'










# Will have to run part one of Hannah's workflow
#
# #################################################
# #################################################
# #################################################
# # Vijay work:
# #
#
# ln -sf /net/shendure/vol8/projects/single.cell.combinatorics/nobackup/sciATAC_EB_MKII/ vijayWork
#
# ls -ltrh /net/shendure/vol8/projects/single.cell.combinatorics/nobackup/sciATAC_EB_MKII/*
#
# # Download
# WDIR="2019_sciATAC_analysis/results/gbonora/20190208_sciATAC_mouseDiff/vijayWork/"
# rsync -auvPKLn --exclude="*.counts.txt" --include "*.R" --include="*.py" --include="*.txt" --include="*.by.sample" --include="*.pdf" --include="*.csv" --include='*/knee-plots/*' --include='*/final-output/*' --prune-empty-dirs --include='*/' --exclude='*' gs:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/
# rsync -auvPKLn --include="*.counts.txt.chr12" --prune-empty-dirs --include='*/' --exclude='*' gs:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/










#################################################
#################################################
#################################################
# Set up data description files

cd "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"

# Setup DATADIR
cd "${HOME}"/proj/"${PROJDIR}"/data/"${READDIR}"
pwd -P > "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/DATADIR
cat "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/DATADIR
# /net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouseDiff_Nmasked

cd "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"
# ls -1 "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/vijayWork/sciATAC_EB_MKII.split.q10.sort.bam
# echo "sciATAC_EB_MKII" > libIDs
echo "sciATAC_mouseDiff" > libIDs
cat libIDs
# sciATAC_mouseDiff










##################################################
##################################################
##################################################
# Use Hannah's prepocessing script
# Will have to run part one of Hannah's workflow too.

mkdir -p "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/processReads
cd  "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/processReads

# ssh grid
# qlogin -l mfree=16G



##################################################
# README.md
# sci-ATAC-seq Read Processing Pipeline
# The purpose of this pipeline is to get reads from bcl (directly off the illumina machines) to a peak by cell matrix.
# Currently, this pipeline takes care of the first step: running bcl2fastq, cleaning up and correcting barcodes and mapping.
# I am still very actively debugging, so report issues if you have them.

# Upload code, inc. src, picard-tools-1.141, and Trimmomatic-0.36
# Must first copy files from ATAC-Pipeline-Dev folder
WDIR="2019_sciATAC_analysis/results/gbonora/20191105_sciATAC_mouseDiff_Nmasked/processReads"
rsync -avuPKLn --exclude=".DS_Store" ~/NobleLab/proj/"${WDIR}"/ gs:proj/"${WDIR}"/










##################################################
# ## Installation:
# # Currently, this must be run on the lab cluster because it recruits cluster modules.
# # ~~~~
# 20191112 julia/latest not defined & julia/1.1.1 gives errors:
# ERROR: LoadError: Failed to precompile Levenshtein [cb002d51-b850-5775-a4f8-1db47be676b8] to /net/noble/vol2/home/gbonora/.julia/compiled/v1.1/Levenshtein/NjaMZ.ji.
module load julia/0.6.0 # julia/1.1.1 # julia/latest
julia
using Pkg # 20191111 https://github.com/JuliaLang/julia/issues/28574
Pkg.clone("https://github.com/hpliner/Levenshtein.jl.git")
Pkg.add("GZip")
Pkg.add("ArgParse")
# # ~~~~
#
# *** OUTPUT ***
# julia> Pkg.clone("https://github.com/hpliner/Levenshtein.jl.git")
# INFO: Initializing package repository /net/noble/vol2/home/gbonora/.julia/v0.6
# INFO: Cloning METADATA from https://github.com/JuliaLang/METADATA.jl
# INFO: Cloning Levenshtein from https://github.com/hpliner/Levenshtein.jl.git
# INFO: Computing changes...
# WARNING: unknown Levenshtein commit 92586fee, metadata may be ahead of package cache
# WARNING: unknown Levenshtein commit f043381e, metadata may be ahead of package cache
# INFO: Cloning cache of Compat from https://github.com/JuliaLang/Compat.jl.git
# INFO: Installing Compat v1.5.1
# INFO: Package database updated
#
# julia> Pkg.add("GZip")
# INFO: Cloning cache of GZip from https://github.com/JuliaIO/GZip.jl.git
# INFO: Installing GZip v0.4.0
# INFO: Package database updated
# INFO: METADATA is out-of-date — you may not have the latest version of GZip
# INFO: Use `Pkg.update()` to get the latest versions of your packages
#
# julia> Pkg.add("ArgParse")
# INFO: Cloning cache of ArgParse from https://github.com/carlobaldassi/ArgParse.jl.git
# INFO: Cloning cache of TextWrap from https://github.com/carlobaldassi/TextWrap.jl.git
# INFO: Installing ArgParse v0.6.1
# INFO: Installing TextWrap v0.3.0
# INFO: Package database updated
# INFO: METADATA is out-of-date — you may not have the latest version of ArgParse
# INFO: Use `Pkg.update()` to get the latest versions of your packages
#
# julia> Pkg.update()
# INFO: Updating METADATA...
# INFO: Cloning cache of Levenshtein from https://github.com/rawrgrr/Levenshtein.jl.git
# INFO: Updating Levenshtein master...
# INFO: Computing changes...
# INFO: No packages to install, update or remove
# *** OUTPUT ***

# 20191111
Pkg.installed()
# Dict{String,Union{Nothing, VersionNumber}} with 3 entries:
#   "Levenshtein" => v"0.2.0+"
#   "GZip"        => v"0.5.0"
#   "ArgParse"    => v"0.6.2"










##################################################
# ## Part 1 - Flowcell to QC-ed Sorted BAM
# ### Basic usage:
#
# Make a script called `runcall.sh` with the following contents:
# NOTE: These defaults assume that you are mapping to hg19 and that you ran the run on a NextSeq.
# If not, use `--help` to see alternative flags.
# ~~~~
# #$ -pe serial 10
# #$ -l mfree=10G
#
# module load samtools/latest
# module load python/2.7.3
# module load julia/latest
# module load pysam/0.8.1
# module load coreutils/8.24
# python path/to/runall.py -R [Flowcell run directory] -O [Path to output folder]
#   -P [Prefix you want on output files] -p 10 --keep_intermediates
# ~~~~
#
# To execute, run `qsub runcall.sh`
#
# The first argument of runall.py is the full path to your flowcell data. The second argument is the full path to where you want the output to go. The third argument is optional, and gives each output file an experiment prefix. --keep_intermediates tells the pipeline to keep all the intermediate files which lets you restart where you were when something goes wrong. A good idea to have it if there's a chance something will get messed up.
#
# (example: `runall.py -R /net/shendure/vol9/seq/170607_NS500488_0394_AHVHL5BGX2 -O /net/trapnell/vol1/hannah/Projects/tests/ -P hifT -p 10`)
#
# For more details on further arguments, run `python runall.py --help`
#
# ### Output
# 1. A folder called `fastq` that holds:
#   a. The original fastq output from bcl2fastq along with some stats from illumina.
# 2. bcl2fastq_log.txt which holds the usual output from bcl2fastq - cluster density etc.
# 3. Prefix.log, that has times when processes started and ended.
# 4. Prefix.split.q10.sort.bam/bai - an indexed and sorted bam with all mapped reads with quality above 10.
# 5. A folder called `qc_info` that contains some QC stats and plots

cd  "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/processReads
# Copy runall.py to here

# Check trimmomatic runs -- called by runall.py:
java -Xmx1G -jar /net/noble/vol2/home/gbonora/proj/2019_sciATAC_analysis/results/gbonora/20191105_sciATAC_mouseDiff_Nmasked/processReads/Trimmomatic-0.36/trimmomatic-0.36.jar

# SRCDIR=/net/shendure/vol9/seq/NEXTSEQ/180507_NS500488_0611_AH2KHCAFXY
SRCDIR=/net/shendure/vol9/seq/NEXTSEQ/180822_NS500488_0676_AH5FCGAFXY

LIBNAME=$( cat "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/libIDs )
DATADIR=$( cat "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/DATADIR )

# REFINDEX="/net/trapnell/vol1/genomes/human/hg19/bowtie2/"
REFINDEX="/net/noble/vol2/home/gbonora/refdata/mm10.Nmasked.Cast129/bowtie2-index/mm10.Cast129_Nmasked"

mkdir ${DATADIR}/readProcessing
ln -sf ${DATADIR}/readProcessing .

# Delete previous work.
# rm -rf /net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouseDiff_Nmasked/readProcessing/*

THREADS=10
# THREADS=2

cat << EOF > runcall.sh
#!/bin/bash -x
#\$ -cwd
#\$ -j y
#$ -pe serial ${THREADS}
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
-R "${SRCDIR}" \
-O "${DATADIR}/readProcessing" \
-P "${LIBNAME}" \
-G "${REFINDEX}" \
-p ${THREADS} \
--keep_intermediates
# --force_overwrite_all

EOF

# To execute, run `qsub runcall.sh`
qsub runcall.sh










##################################################
## Part 2 - QC-ed Sorted BAM(s) to sparse matrix
### Basic usage:
#
# Make a script called `mergecall.sh` with the following contents:
# ~~~~
#$ -l mfree=100G
# module load python/2.7.3
# module load pysam/0.8.3
# module load julia/latest
# module load coreutils/8.24
# python path/to/mergeall.py -B [list of split.q10.sort.bams] \
# -O [Path to output folder] \
# -P [Prefix you want on output files] \
# -C [path to file with each of the valid barcodes allowed in the experiment] \
# --keep_intermediates
# ~~~~
#
# To execute, run `qsub mergecall.sh`
#
# The first argument of mergeall.py is a list of the full paths to the BAMs created as output of Part 1.
# The second argument is the full path to where you want the output to go.
# The third argument is optional, and gives each output file an experiment prefix.
# The last argument is the path to a text file with all the allowed barcodes for your experiment, one barcode per line (see below).
#
# For more details on further arguments, run `python mergeall.py --help`
#
#
# ### Output
# a bunch of stuff
#
# There is a helper function in src called `barcode_helper.R` that will assemble your barcode file for you. To run it do:
# `Rscript --vanilla path/to/barcode_helper.R [LIST OF BARCODE COMBOS]`
# The list of barcodes will be letter plus index number combos that were used in your experiment,
# and that correspond with the PCR and Tn5 plates.
# For example, if I used Column 5 and row A for one plate and column 7 and row B for a second plate,
# I would run `Rscript --vanilla path/to/barcode_helper.R A i5 B i7`



##################################################
# Biallelic processing:
# 1. Use mouse black list https://sites.google.com/site/anshulkundaje/projects/blacklists
# 2. Override UMI cutoff to be 100

cd  "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/processReads
# copy mergeall.py to here

# Copy barcode info from old work directory
OLDWDIR="$HOME/proj/2019_sciATAC_analysis/results/gbonora/20190208_sciATAC_mouseDiff/"
cp ${OLDWDIR}/vijayWork/barcodes.txt .
cp ${OLDWDIR}/vijayWork/bc_combos.txt .

LIBNAME=$( cat "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/libIDs )
DATADIR=$( cat "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/DATADIR )

BARCODES="bc_combos.txt"
BLACKLIST="MM10_BLACKLIST" # make sure that this file is in place /src/mm10.blacklist.bed
UMITHRESH=100

INDIR=${DATADIR}/readProcessing
# BAMFILELIST="BAMlist.txt"
# echo ${INDIR}/"${LIBNAME}".split.q10.sorted.bam > ./${BAMFILELIST}
BAMFILE=${INDIR}/"${LIBNAME}".split.q10.sort.bam

OUTDIR=${DATADIR}/readProcessing/${LIBNAME}

cat << EOF > "mergecall.sh"
#!/bin/bash
#\$ -l mfree=64G
#\$ -cwd
#\$ -j y

source $HOME/.bashrc
source $HOME/.bash_profile

# export PATH=/net/noble/vol2/home/gbonora/anaconda2/bin:$PATH # condaon
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
printf "\n\nstart: %s\n\n" "\$( date )"

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

python ./mergeall.GB.py \
-B "${BAMFILE}" \
-O "${OUTDIR}" \
-P "${LIBNAME}" \
-C "${BARCODES}" \
--override_reads_per_cell "${UMITHRESH}" \
# --force_overwrite_all \
--keep_intermediates

# -R "${BLACKLIST}"  <-- THIS IS THE DEFAULT NOW

printf "\n\nend: %s\n\n" "\$( date )"
EOF

qsub mergecall.sh

#  rm ${DATADIR}/readProcessing/sciATAC_mouseDiff/*.log;


############
# Check

cat ${DATADIR}/readProcessing/sciATAC_mouseDiff.log

ls -ltrh ${DATADIR}/readProcessing



##################################################
# 20200326
# Prepare for downstream processing

cd  "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/processReads

# ln -sf /net/noble/vol2/home/gbonora/proj/2019_sciATAC_analysis/results/gbonora/20190208_sciATAC_mouseDiff/vijayWork/sci-ATAC_out/Plate1_mappings.txt processReads/readProcessing/sciATAC_mouseDiff/.

# Copy python code from old work directory
OLDWDIR="$HOME/proj/2019_sciATAC_analysis/results/gbonora/20190208_sciATAC_mouseDiff/"
# cp -pr ${OLDWDIR}/vijayWork/barcodes.txt .
# cp -pr ${OLDWDIR}/vijayWork/bc_combos.txt .
cp -pr ${OLDWDIR}/vijayWork/sci-ATAC_out/*.py .

head ${DATADIR}/readProcessing/sciATAC_mouseDiff/sciATAC_mouseDiff.counts.txt
# chr4_115979368_115979896        CTGAAGCTGAGCGCGATGCTCTAACTCGGTACTGAC    3
# chr3_123690078_123691589        CGGCTATGTTATTGAGGCTAGTTCTAGACAGGACGT    2
# chr10_95563220_95564827 TAATGCGCGCGCGTTCATCTCTAACTCGTATAGCCT    2

python parse_sciATAC.py <( cut -f 2 ${DATADIR}/readProcessing/sciATAC_mouseDiff/sciATAC_mouseDiff.counts.txt ) > ${DATADIR}/readProcessing/sciATAC_mouseDiff/cellIDs.txt
cat ${DATADIR}/readProcessing/sciATAC_mouseDiff/cellIDs.txt | sort | uniq > ${DATADIR}/readProcessing/sciATAC_mouseDiff/cellIDs.unique.txt
cat ${DATADIR}/readProcessing/sciATAC_mouseDiff/cellIDs.unique.txt | wc -l
# 5586

# gzip -f /net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouseDiff_Nmasked/readProcessing/sciATAC_mouseDiff/sciATAC_mouseDiff.counts.txt








##################################################
##################################################
##################################################
# Get results

WDIR="2019_sciATAC_analysis/results/gbonora/20191105_sciATAC_mouseDiff_Nmasked"

# Images
rsync -auvPKLn --include='*.png' --include='*.pdf' --include='*/' --prune-empty-dirs  --exclude='*' gs:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/

# HTML and txt
rsync -auvPKLn --exclude="*counts*" --include='*.txt' --include='*.html' --include='*/' --prune-empty-dirs  --exclude='*' gs:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/

# Counts table and cell IDs
rsync -auvPKLn --include="sciATAC_mouseDiff.counts.txt.gz" --include="cellIDs.unique.txt" --include='*/' --prune-empty-dirs  --exclude='*' gs:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/

# Code
rsync -auvPKLn --include="*.py" --include='*/' --prune-empty-dirs  --exclude='*' gs:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/


# Misc
rsync -auvPKLn --include='20180516_sciATAC_WTandDel1Patski_pilot_segregationResults.txt' --include='*/' --prune-empty-dirs  --exclude='*' gs:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/

# upload code, inc. src, picard-tools-1.141, and Trimmomatic-0.36
WDIR="2019_sciATAC_analysis/results/gbonora/20191105_sciATAC_mouseDiff_Nmasked"
rsync -avuPKLn --exclude=".DS_Store" ~/NobleLab/proj/"${WDIR}"/ gs:proj/"${WDIR}"/
