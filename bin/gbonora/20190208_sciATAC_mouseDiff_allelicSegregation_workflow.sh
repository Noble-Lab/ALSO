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
WORKDIR="20190208_sciATAC_mouseDiff"
READDIR="data_20190208_sciATAC_mouseDiff"


##################################################
# Working folder
mkdir -p "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"
cd  "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"










# #################################################
# #################################################
# #################################################
# # Set up data folder
#
# VOL8HOME=/net/noble/vol8/gbonora/sciOmics
# mkdir -p "${VOL8HOME}"
# cd "${VOL8HOME}"
# pwd -P
#
# # Create a folder for data
# mkdir -p "${VOL8HOME}"/"${PROJDIR}"/data/"${READDIR}"
#
# # Link to PROJDIR
# cd "${HOME}"/proj/"${PROJDIR}"
# ln -sf "${VOL8HOME}"/"${PROJDIR}"/data .
# ls -ltrh data/










#################################################
#################################################
#################################################
# Step 0. Link to mouse raw data

# # Make data folder
# ls -ltrh "${HOME}"/proj/"${PROJDIR}"/data
# mkdir -p "${HOME}"/proj/"${PROJDIR}"/data/"${READDIR}"/rawReads
# cd "${HOME}"/proj/"${PROJDIR}"/data/"${READDIR}"/rawReads

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

# # SRCDIR=/net/shendure/vol9/seq/NEXTSEQ/180507_NS500488_0611_AH2KHCAFXY
# SRCDIR=/net/shendure/vol9/seq/NEXTSEQ/180822_NS500488_0676_AH5FCGAFXY









#################################################
#################################################
#################################################
# Get Hannah pipeline setup
#
# See '20190208_cole-trapnell-lab_ATAC-Pipeline-Dev_githubRepo.sh'










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










# #################################################
# #################################################
# #################################################
# # Set up data description files
#
# cd "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"
#
# # Setup DATADIR
# cd "${HOME}"/proj/"${PROJDIR}"/data/"${READDIR}"
# pwd -P > "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/DATADIR
# cat "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/DATADIR
# # /net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20190208_sciATAC_mouseDiff
#
# cd "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"
# ls -1 "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/vijayWork/sciATAC_EB_MKII.split.q10.sort.bam
# echo "sciATAC_EB_MKII" > libIDs
# cat libIDs
# # sciATAC_EB_MKII










#################################################
#################################################
#################################################
# Check out reads

module load samtools/latest

# NOTE: PE
samtools flagstat "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/vijayWork/sciATAC_EB_MKII.split.q10.sort.bam
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 208717262 + 0 mapped (100.00% : N/A)
# 208717262 + 0 paired in sequencing
# 104358631 + 0 read1
# 104358631 + 0 read2
# 208717262 + 0 properly paired (100.00% : N/A)
# 208717262 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

# NOTE: Sorted by position.
samtools view "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/vijayWork/sciATAC_EB_MKII.split.q10.sort.bam | head -n 4
# TCTCGCGCGAGCGCGATGCAGAGAGGTCGTACTGAC:0  99      chr1    3000321 22      50M     =       3000416 145     GTTGTTATGTCTCCTTTTTCATTTCTGATTTTGTTAATTATAGTACAGTC      AAAAAE6EEEEEEAEEEEE6EEEEAEEEEEEEEEEEAEEAEEAEEEEEEE  AS:i:-5 XS:i:-30        XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:4G45       YS:i:-5 YT:Z:CP
# TCTCGCGCGAGCGCGATGCAGAGAGGTCGTACTGAC:0  147     chr1    3000416 22      50M     =       3000321 -145    GACTTTCTCAAAGAACCAGCTACTAGTTTGGTTGATTCTTTGAATATTTC      EEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAAAAA  AS:i:-5 XS:i:-10        XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:25T24      YS:i:-5 YT:Z:CP
# TCCGGAGATCTATCGCTGTAATTGACCTTATAGCCT:0  163     chr1    3000648 32      50M     =       3000966 368     CATTGTGCCCCATATGTTTGGCTATGTTGTGGATTTATTTTCATTAAACT      AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE  AS:i:0  XS:i:-5 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YS:i:0  YT:Z:CP
# CGGCTATGAGTCCGCTGCTAGTTCTAGACAGGACGT:2  99      chr1    3000835 40      50M     =       3001007 222     TCTGAAAGGATGCATGGGAAAATTTCAATATTTTTGTATCTGTTGAGGAC      AAAAAEEEEEEEEEEEEEEEEEEEEEAEEAE/EEEEEEEEEEEEEEEEEE  AS:i:0  XS:i:-10        XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YS:i:-18        YT:Z:CP







#################################################
#################################################
#################################################
# Alellic segregation

mkdir -p "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/allelicSegregation
cd  "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/allelicSegregation

INDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/vijayWork
DATADIR=$(cat ../DATADIR)
OUTDIR=${DATADIR}/segregatedReads
mkdir -p ${OUTDIR}
ln -sf ${OUTDIR} ..


#################################################
# sort by name first

# ssh grid
# qlogin -pe serial 8 -l mfree=4G

module load samtools/1.7

samtools sort -n -m 4G --threads 8 -T sciATAC_EB_MKII.nsorting \
-o ${OUTDIR}/sciATAC_EB_MKII.split.q10.nsorted.bam  \
"${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/vijayWork/sciATAC_EB_MKII.split.q10.sort.bam

ls -ltrh ${OUTDIR}

samtools view ${OUTDIR}/sciATAC_EB_MKII.split.q10.sort.nsorted.bam | head -n 4
# AGCGATAGAACGAATTCGAAGCCTACGAATAGAGGC:6  99      chr15   58066016        42      50M     =       58066224        258     ACCCACTGTCCTTATACAAGAGAGACAGACAGACAGACATGCAAAAAGGT      AAAAAEEEEEEEEEEEEEEEEEEEEEEEAEEE/EEEEEEEEEEEEEEEEE      AS:i:0  XN:i:0  XM:i:0 XO:i:0   XG:i:0  NM:i:0  MD:Z:50 YS:i:0  YT:Z:CP
# AGCGATAGAACGAATTCGAAGCCTACGAATAGAGGC:6  147     chr15   58066224        42      50M     =       58066016        -258    AAGATAAACCTCAGGAAAATGTCAAAACCAAAAGCAGTCCAGTATTAAAG      AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEAEEEEEAAAAA      AS:i:0  XN:i:0  XM:i:0 XO:i:0   XG:i:0  NM:i:0  MD:Z:50 YS:i:0  YT:Z:CP
# AGCGATAGAACGAATTCGAAGCCTACGACAGGACGT:5  99      chr5    145103845       40      50M     =       145103909       114     CACACACACACACACATTCACACTCACACGCACAGACATGCACGGGGCCT      AAAAAEEEEEEEEEEEEAEEEAEEEEEEEEEEEEE/AEEEEEEEEEAEE/      AS:i:0  XN:i:0  XM:i:0 XO:i:0   XG:i:0  NM:i:0  MD:Z:50 YS:i:-15        YT:Z:CP
# AGCGATAGAACGAATTCGAAGCCTACGACAGGACGT:5  147     chr5    145103909       40      50M     =       145103845       -114    CGGCTGAGACTAAAACTCCTCTCATAGTGGCGGCTTACACAGGGCACGAC      ///A////<E/AAA///A//E/E/////A////E//EE/E//A/EA//AA      AS:i:-15        XN:i:0 XM:i:5   XO:i:0  XG:i:0  NM:i:5  MD:Z:0T14T2A11T15T3     YS:i:0  YT:Z:CP


#################################################
# Now segregate

SNPTHRESH=1
QTHRESH=0

old_IFS=$IFS
IFS=$'\n'
libIDs=($(cat ../libIDs)) # libIDs to array
IFS=$old_IFS

# *** NOTE: Need mm10 snps
# # snpFile=/net/noble/vol1/home/wenxiu/proj/blau/results/wenxiu/2012-08-14_spretus-genome/spretus.ATG.vcf.Fan.combined.gz
# # snpFileName=spretus.ATG.vcf.Fan.combined.gz
# # snpFile=$HOME/refdata/mm10pseudoSpretus/spretus.snps.vcf.gz
# # snpFileName=spretus.snps.vcf.gz
# snpFile=$HOME/refdata/mm10pseudoSpretus/spretus.patski.snps.vcf.gz
# snpFileName=spretus.patski.snps.vcf.gz
# snpFile=$HOME/refdata/mm10pseudoCastBWAindexV2/2016-06-27_CAST-genome-mm10/"CAST.snps.vcf.gz"
# snpFileName="CAST.snps.vcf.gz"
# Use *** CORRECT *** Sanger snps:
refSnpFile=$HOME/refdata/mm10CAST129SNPs/129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz
# refSnpFile=$HOME/refdata/mm10CAST129SNPs/129S5SvEvBrd.mgp.v5.snps.dbSNP142.PASS.vcf.gz
# refSnpFile=$HOME/refdata/mm10CAST129SNPs/129P2_OlaHsd.mgp.v5.snps.dbSNP142.vcf.gz
altSnpFile=$HOME/refdata/mm10CAST129SNPs/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz

refSnpFileName=$(basename ${refSnpFile})
altSnpFileName=$(basename ${altSnpFile})


for (( i=0; i<${#libIDs[@]}; i++ )); do
libID=${libIDs[$i]}
SUBINDIR=${OUTDIR}
SUBOUTDIR=${OUTDIR}
mkdir -p ${SUBOUTDIR}

BAMFILES=($(ls -1 ${SUBINDIR}/*.split.q10.nsorted.bam))

for (( j=0; j<${#BAMFILES[@]}; j++ )); do
BAMFILE=${BAMFILES[$j]}
BASENAME=$(basename ${BAMFILE} ".bam")

jobfile="${libID}.${BASENAME}.job"
echo ${jobfile}

cat << EOF > "${jobfile}"
#!/bin/bash
#\$ -l mfree=30G
#\$ -cwd
#\$ -j y

source $HOME/.bashrc
source $HOME/.bash_profile

export PATH=/net/noble/vol2/home/gbonora/anaconda2/bin:$PATH # condaon
# source activate sciHiC
export PYTHONPATH=""

source /etc/profile.d/modules.sh
module load modules modules-init modules-gs modules-noble
# module load bedtools/latest
module load samtools/latest

hostname
printf "\n\nstart: %s\n\n" "\$( date )"

if [[ \${TMPDIR:-} != "" ]]; then

printf "\n\nCopying files: %s\n\n" "\$( date )"
cp ${BAMFILE} \${TMPDIR:-}/${BASENAME}".bam"
cp ${refSnpFile} \${TMPDIR:-}/${refSnpFileName}
cp ${altSnpFile} \${TMPDIR:-}/${altSnpFileName}

printf "\n\nSegregating: %s\n\n" "\$( date )"
# python ./segregate-sciATAC.refAndAltSNPs.SE.SNPthresh.Qthresh.py \${TMPDIR:-}/${BASENAME}".bam" \${TMPDIR:-}/${refSnpFileName} \${TMPDIR:-}/${altSnpFileName} \${TMPDIR:-}/$BASENAME $SNPTHRESH $QTHRESH
python ./segregate-sciATAC.refAndAltSNPs.PE.SNPthresh.Qthresh.py \${TMPDIR:-}/${BASENAME}".bam" \${TMPDIR:-}/${refSnpFileName} \${TMPDIR:-}/${altSnpFileName} \${TMPDIR:-}/$BASENAME $SNPTHRESH $QTHRESH

printf "\n\nCopying back results: %s\n\n" "\$( date )"
cp -t $SUBOUTDIR/ \${TMPDIR:-}/$BASENAME.ref.bam  \${TMPDIR:-}/$BASENAME.alt.bam  \${TMPDIR:-}/$BASENAME.ambig.bam  \${TMPDIR:-}/$BASENAME.contra.bam

else

# python ./segregate-sciATAC.SE.refAndAltSNPs.SNPthresh.Qthresh.py ${BAMFILE} ${refSnpFile} ${altSnpFile} $SUBOUTDIR/$BASENAME $SNPTHRESH $QTHRESH
python ./segregate-sciATAC.PE.refAndAltSNPs.SNPthresh.Qthresh.py ${BAMFILE} ${refSnpFile} ${altSnpFile} $SUBOUTDIR/$BASENAME $SNPTHRESH $QTHRESH

fi

printf "\n\nend: %s\n\n" "\$( date )"
EOF

done
done

chmod u+x *.job

ls -1 *.job | wc -l
# 1



############
# Submit jobs
#sshgrid
while read JOBFILE; do
echo "${JOBFILE}"
qsub "${JOBFILE}"
done  < <(ls -1 *.job)

# # OR just run
JOBFILE=sciATAC_EB_MKII.split.q10.nsorted.job
./"${JOBFILE}"  2>&1 | tee "${JOBFILE}".out



############
# Check

tail -n28 sciATAC_EB_MKII.split.q10.sort.job.out  | head -n 24 > sciATAC_EB_MKII.split.q10.sort.SE.segregation.summary.txt
# 208717262 total reads.
# 29356625 reads assigned as ref based on one end.
# 32811649 reads assigned as alt based on one end.
# 1183384 reads contradictory or unresolvable.
# 145365604 reads unassigned because ambiguous.
#
# 87522221 (0.84%) SNP-associated bases in E1s.
# 10342072073 (99.16%) SNP-unassociated bases in E1s.
# 42052521 total SNP base calls assigned to 129 allele in E1s.
# 45104880 total SNP base calls assigned to Cast in E1s.
#
# 9990822 total # of Sanger 129 SNP base calls considered in E1s.
# 4698513 129 base calls at 129 SNP loci in E1s.
# 5237624 B6 base calls (inferred as Cast-eminating) at 129 SNP loci in E1s.
# 39863 non-129/B6/Cast base calls at 129 SNP loci in E1s.
# 0 Low quality base calls at 129 SNP sites in E1s.
# 14822 Ns at 129 SNP sites in E1s.
#
# 77531399 total # of Sanger Cast SNPs considered in E1s.
# 39867256 Cast base calls at Cast SNP loci in E1s.
# 37354008 B6 base calls (inferred as 129-eminating) at Cast SNP loci in E1s.
# 280811 non-129/B6/Cast bases at Cast SNP loci in E1s.
# 0 Low quality base calls at Cast SNP sites in E1s.
# 29324 Ns at Cast SNP sites in E1s.

tail -n51 sciATAC_EB_MKII.split.q10.nsorted.job.out | head -n 43 > sciATAC_EB_MKII.split.q10.nsorted.PE.segregation.summary.txt
# 104358631 total reads.
# 18804104 reads assigned as ref based on one end.
# 21101305 reads assigned as alt based on one end.
# 918951 reads contradictory or unresolvable.
# 61168931 reads unassigned because ambiguous.
#
# 43762976 (0.84%) SNP-associated bases in E1s.
# 5171031969 (99.16%) SNP-unassociated bases in E1s.
# 21030719 total SNP base calls assigned to 129 allele in E1s.
# 22550678 total SNP base calls assigned to Cast in E1s.
#
# 4995148 total # of Sanger 129 SNP base calls considered in E1s.
# 2349128 129 base calls at 129 SNP loci in E1s.
# 2618811 B6 base calls (inferred as Cast-eminating) at 129 SNP loci in E1s.
# 19782 non-129/B6/Cast base calls at 129 SNP loci in E1s.
# 0 Low quality base calls at 129 SNP sites in E1s.
# 7427 Ns at 129 SNP sites in E1s.
#
# 38767828 total # of Sanger Cast SNPs considered in E1s.
# 19931867 Cast base calls at Cast SNP loci in E1s.
# 18681591 B6 base calls (inferred as 129-eminating) at Cast SNP loci in E1s.
# 139665 non-129/B6/Cast bases at Cast SNP loci in E1s.
# 0 Low quality base calls at Cast SNP sites in E1s.
# 14705 Ns at Cast SNP sites in E1s.
#
# 43759245 (0.84%) SNP-associated bases in E2s.
# 5171040104 (99.16%) SNP-unassociated bases in E2s.
# 21021802 total SNP base calls assigned to 129 allele in E2s.
# 22554202 total SNP base calls assigned to Cast in E2s.
#
# 4995674 total # of Sanger 129 SNP base calls considered in E2s.
# 2349385 129 base calls at 129 SNP loci in E2s.
# 2618813 B6 base calls (inferred as Cast-eminating) at 129 SNP loci in E2s.
# 20081 non-129/B6/Cast base calls at 129 SNP loci in E2s.
# 0 Low quality base calls at 129 SNP sites in E2s.
# 7395 Ns at 129 SNP sites in E2s.
#
# 38763571 total # of Sanger Cast SNPs considered in E2s.
# 19935389 Cast base calls at Cast SNP loci in E2s.
# 18672417 B6 base calls (inferred as 129-eminating) at Cast SNP loci in E2s.
# 141146 non-129/B6/Cast bases at Cast SNP loci in E2s.
# 0 Low quality base calls at Cast SNP sites in E2s.
# 14619 Ns at Cast SNP sites in E2s.


############
# Sort and index output files

# ssh grid
# qlogin -pe serial 8 -l mfree=4G

module load samtools/1.7

while read FILE; do

BASEN=$(basename ${FILE} ".bam")
echo ${BASEN}

samtools sort -m 4G --threads 8 -T ${BASEN}.sorting \
-o "${OUTDIR}/${BASEN/.nsorted./.sorted.}.bam" ${FILE}
echo ""

samtools index -@ 8 "${OUTDIR}/${BASEN/.nsorted./.sorted.}.bam"
echo ""

done < <( ls -1 ${OUTDIR}/*.nsorted.*.bam  )


############
# Check output files

ls -ltrh ${OUTDIR}
# -rw-rw-r-- 1 gbonora noblelab   41M Feb 16 21:24 sciATAC_EB_MKII.split.q10.sort.contra.bam
# -rw-rw-r-- 1 gbonora noblelab  3.8G Feb 16 21:24 sciATAC_EB_MKII.split.q10.sort.ambig.bam
# -rw-rw-r-- 1 gbonora noblelab 1009M Feb 16 21:24 sciATAC_EB_MKII.split.q10.sort.alt.bam
# -rw-rw-r-- 1 gbonora noblelab  834M Feb 16 21:24 sciATAC_EB_MKII.split.q10.sort.ref.bam
#
# -rw-rw-r-- 1 gbonora noblelab   79M Feb 20 01:52 sciATAC_EB_MKII.split.q10.nsorted.contra.bam
# -rw-rw-r-- 1 gbonora noblelab  3.8G Feb 20 01:52 sciATAC_EB_MKII.split.q10.nsorted.ambig.bam
# -rw-rw-r-- 1 gbonora noblelab  1.7G Feb 20 01:52 sciATAC_EB_MKII.split.q10.nsorted.alt.bam
# -rw-rw-r-- 1 gbonora noblelab  1.4G Feb 20 01:52 sciATAC_EB_MKII.split.q10.nsorted.ref.bam
# -rw-rw-r-- 1 gbonora noblelab  6.7G Feb 19 13:46 sciATAC_EB_MKII.split.q10.nsorted.bam

samtools view ${OUTDIR}/sciATAC_EB_MKII.split.q10.sort.ref.bam
samtools view ${OUTDIR}/sciATAC_EB_MKII.split.q10.sort.alt.bam

samtools view ${OUTDIR}/sciATAC_EB_MKII.split.q10.nsorted.ref.bam
samtools view ${OUTDIR}/sciATAC_EB_MKII.split.q10.nsorted.alt.bam

samtools view ${OUTDIR}/sciATAC_EB_MKII.split.q10.sorted.ref.bam
samtools view ${OUTDIR}/sciATAC_EB_MKII.split.q10.sorted.alt.bam










##################################################
##################################################
##################################################
# Adapt Hannah's prepocessing script for allelic data
#
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
WDIR="2019_sciATAC_analysis/results/gbonora/20190208_sciATAC_mouseDiff/allelicReprocessing"
rsync -avuPKLn --exclude=".DS_Store" ~/NobleLab/proj/"${WDIR}"/ gs:proj/"${WDIR}"/



# ## Installation:
# # Currently, this must be run on the lab cluster because it recruits cluster modules.
# # ~~~~
module load julia/latest
julia
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
##################################################
##################################################
# Allelic processing RUN1: 20190222

mkdir -p "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/allelicReprocessing
cd  "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/allelicReprocessing

mkdir ${DATADIR}/allelicWork
ln -sf ${DATADIR}/allelicWork .

# sync mergall.py and src subfolder

cp ../vijayWork/barcodes.txt .
cp ../vijayWork/bc_combos.txt .
# diff bc_combos.txt /net/noble/vol2/home/gbonora/proj/2019_sciATAC_analysis/results/gbonora/20190208_sciATAC_mouseDiff/vijayWork/bc_combos.txt

DATADIR=$(cat ../DATADIR)
INDIR=${DATADIR}/segregatedReads

old_IFS=$IFS
IFS=$'\n'
libIDs=($(cat ../libIDs)) # libIDs to array
IFS=$old_IFS

BARCODES="bc_combos.txt"





##################################################
# Process reads and tally counts in allelic peaks

for (( i=0; i<${#libIDs[@]}; i++ )); do
libID=${libIDs[$i]}
for theAllele in alt ref; do

ALLELICPREFIX="${libID}.${theAllele}"

OUTDIR=${DATADIR}/allelicWork/${ALLELICPREFIX}
mkdir -p ${OUTDIR}
# ln -sf ${OUTDIR} .

# BAMFILELIST="BAMlist.${theAllele}.txt"
# echo ${INDIR}/"${libID}".split.q10.sorted.${theAllele}.bam > ./${BAMFILELIST}
BAMFILE=${INDIR}/"${libID}".split.q10.sorted.${theAllele}.bam

jobfile="mergecall.${theAllele}.job"
echo ${jobfile}

cat << EOF > "${jobfile}"
#!/bin/bash
#\$ -l mfree=50G
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
module load julia/latest
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

python ./mergeall.py \
-B "${BAMFILE}" \
-O "${OUTDIR}" \
-P "${ALLELICPREFIX}" \
-C "${BARCODES}" \
--keep_intermediates

printf "\n\nend: %s\n\n" "\$( date )"
EOF

done
done

chmod u+x mergecall.*.job



############
# Submit jobs
#sshgrid
while read JOBFILE; do
echo "${JOBFILE}"
qsub "${JOBFILE}"
done  < <(ls -1 mergecall.*.job)

# rm mergecall*.job.o*; rm sciATAC_EB_MKII.alt/*.log; rm sciATAC_EB_MKII.ref/*.log


############
# Check

cat ${DATADIR}/allelicWork/sciATAC_EB_MKII.alt/sciATAC_EB_MKII.alt.log
cat ${DATADIR}/allelicWork/sciATAC_EB_MKII.ref/sciATAC_EB_MKII.ref.log

ls -ltrh ${DATADIR}/allelicWork/sciATAC_EB_MKII.alt
ls -ltrh ${DATADIR}/allelicWork/sciATAC_EB_MKII.ref










##################################################
# Tally allelic counts in biallelic peaks

ALTINDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/vijayWork/sci-ATAC_out/macs_output/
BIALLELICMACSFILE=${ALTINDIR}/sci-ATAC_EB_MKII_macs_peaks.narrowPeak

for (( i=0; i<${#libIDs[@]}; i++ )); do
libID=${libIDs[$i]}
for theAllele in alt ref; do

ALLELICPREFIX="${libID}.${theAllele}"

OUTDIR=${DATADIR}/allelicWork/${ALLELICPREFIX}
# mkdir -p ${OUTDIR}
# ln -sf ${OUTDIR} .


jobfile="biallelicPeakCounts.${theAllele}.job"
echo ${jobfile}

cat << EOF > "${jobfile}"
#!/bin/bash
#\$ -l mfree=50G
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
module load julia/latest
module load coreutils/8.24

hostname
printf "\n\nstart: %s\n\n" "\$( date )"

python ./mergeall.GB.allelicReadsInBiallelicPeaks.py \
-M "${BIALLELICMACSFILE}" \
-O "${OUTDIR}" \
-P "${ALLELICPREFIX}" \
--keep_intermediates

printf "\n\nend: %s\n\n" "\$( date )"
EOF

done
done

chmod u+x biallelicPeakCounts.*.job



############
# Submit jobs
#sshgrid
while read JOBFILE; do
echo "${JOBFILE}"
qsub "${JOBFILE}"
done  < <(ls -1 biallelicPeakCounts.*.job)

# rm biallelicPeakCounts.*.job.o*


############
# Check

cat ${DATADIR}/allelicWork/sciATAC_EB_MKII.alt/sciATAC_EB_MKII.alt.log
cat ${DATADIR}/allelicWork/sciATAC_EB_MKII.ref/sciATAC_EB_MKII.ref.log

ls -ltrh ${DATADIR}/allelicWork/sciATAC_EB_MKII.alt
ls -ltrh ${DATADIR}/allelicWork/sciATAC_EB_MKII.ref

cat ${DATADIR}/allelicWork/sciATAC_EB_MKII.alt/sciATAC_EB_MKII.alt.counts.biallelicPeaks.txt | awk 'END{print NF, NR}'
# 3 4,793,287
cat ${DATADIR}/allelicWork/sciATAC_EB_MKII.ref/sciATAC_EB_MKII.ref.counts.biallelicPeaks.txt | awk 'END{print NF, NR}'
# 3 4,233,308

head -n 3 ${DATADIR}/allelicWork/sciATAC_EB_MKII.alt/sciATAC_EB_MKII.alt.counts.biallelicPeaks.txt
# chr7_82648214_82648994  CTGAAGCTTAAGTCCTGAAGGCGAGAGCGTACTGAC    2
# chr2_66439727_66440927  TCCGGAGATGATATTGCGTACGATCATCTAATCTTA    1
# chr4_126428979_126430257        TAATGCGCAGTTCTCATGCTCTAACTCGTATAGCCT    1




##################################################
# Merge allelic counts in biallelic peaks

for (( i=0; i<${#libIDs[@]}; i++ )); do
libID=${libIDs[$i]}

ALLELICPREFIX="${libID}.allelic"
echo ${ALLELICPREFIX}

OUTDIR=${DATADIR}/allelicWork/${ALLELICPREFIX}
mkdir -p ${OUTDIR}
# ln -sf ${OUTDIR} .

cp ../vijayWork/sci-ATAC_out/Plate1_mappings.txt ${OUTDIR}

ALTPREFIX="${libID}.alt"
ALTDIR=${DATADIR}/allelicWork/${ALTPREFIX}

REFPREFIX="${libID}.ref"
REFDIR=${DATADIR}/allelicWork/${REFPREFIX}

cat <( cat ${ALTDIR}/${ALTPREFIX}.counts.biallelicPeaks.txt | awk '{printf("%s\t%s_alt\t%d\n", $1, $2, $3)}' ) \
<( cat ${REFDIR}/${REFPREFIX}.counts.biallelicPeaks.txt | awk '{printf("%s\t%s_ref\t%d\n", $1, $2, $3)}' ) \
> ${OUTDIR}/${ALLELICPREFIX}.counts.biallelicPeaks.txt

done



############
# Check

ls -ltrh ${OUTDIR}/
cat ${OUTDIR}/${ALLELICPREFIX}.counts.biallelicPeaks.txt | awk 'END{print NF, NR}'
# 3 9,026,595










##################################################
##################################################
##################################################
# Allelic processing RUN2: 20190313
# 1. Use mouse black list https://sites.google.com/site/anshulkundaje/projects/blacklists
# 2. Override UMI cutoff to be 100

mkdir -p "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/allelicReprocessing
cd  "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/allelicReprocessing

mkdir ${DATADIR}/allelicWork.V2
ln -sf ${DATADIR}/allelicWork.V2 .

# cp ../vijayWork/barcodes.txt .
# sync mergall.py and src subfolder

DATADIR=$(cat ../DATADIR)
INDIR=${DATADIR}/segregatedReads

old_IFS=$IFS
IFS=$'\n'
libIDs=($(cat ../libIDs)) # libIDs to array
IFS=$old_IFS

BARCODES="bc_combos.txt"
BLACKLIST="MM10_BLACKLIST" # make sure that this file is in place /src/mm10.blacklist.bed
UMITHRESH=100






##################################################
# Reprocess reads and tally counts in allelic peaks

for (( i=0; i<${#libIDs[@]}; i++ )); do
libID=${libIDs[$i]}
for theAllele in alt ref; do

ALLELICPREFIX="${libID}.${theAllele}"

OUTDIR=${DATADIR}/allelicWork.V2/${ALLELICPREFIX}
mkdir -p ${OUTDIR}
# ln -sf ${OUTDIR} .

# BAMFILELIST="BAMlist.${theAllele}.txt"
# echo ${INDIR}/"${libID}".split.q10.sorted.${theAllele}.bam > ./${BAMFILELIST}
BAMFILE=${INDIR}/"${libID}".split.q10.sorted.${theAllele}.bam

jobfile="mergecall.${theAllele}.job"
echo ${jobfile}

cat << EOF > "${jobfile}"
#!/bin/bash
#\$ -l mfree=16G
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
module load julia/latest
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
-P "${ALLELICPREFIX}" \
-C "${BARCODES}" \
--override_reads_per_cell "${UMITHRESH}" \
--force_overwrite_all \
--keep_intermediates

# -R "${BLACKLIST}"  <-- THIS IS THE DEFAULT NOW

printf "\n\nend: %s\n\n" "\$( date )"
EOF

done
done

chmod u+x mergecall.*.job



############
# Submit jobs
#sshgrid
while read JOBFILE; do
echo "${JOBFILE}"
qsub "${JOBFILE}"
done  < <(ls -1 mergecall.*.job)

# rm mergecall*.job.o*; rm sciATAC_EB_MKII.alt/*.log; rm sciATAC_EB_MKII.ref/*.log


############
# Check

cat ${DATADIR}/allelicWork.V2/sciATAC_EB_MKII.alt/sciATAC_EB_MKII.alt.log
cat ${DATADIR}/allelicWork.V2/sciATAC_EB_MKII.ref/sciATAC_EB_MKII.ref.log

ls -ltrh ${DATADIR}/allelicWork.V2/sciATAC_EB_MKII.alt
ls -ltrh ${DATADIR}/allelicWork.V2/sciATAC_EB_MKII.ref










##################################################
# Tally allelic counts in biallelic peaks

ALTINDIR="${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/vijayWork/sci-ATAC_out/macs_output/
BIALLELICMACSFILE=${ALTINDIR}/sci-ATAC_EB_MKII_macs_peaks.narrowPeak

for (( i=0; i<${#libIDs[@]}; i++ )); do
libID=${libIDs[$i]}
for theAllele in alt ref; do

ALLELICPREFIX="${libID}.${theAllele}"

OUTDIR=${DATADIR}/allelicWork.V2/${ALLELICPREFIX}
# mkdir -p ${OUTDIR}
# ln -sf ${OUTDIR} .


jobfile="biallelicPeakCounts.${theAllele}.job"
echo ${jobfile}

cat << EOF > "${jobfile}"
#!/bin/bash
#\$ -l mfree=16G
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
module load julia/latest
module load coreutils/8.24

hostname
printf "\n\nstart: %s\n\n" "\$( date )"

python ./mergeall.GB.allelicReadsInBiallelicPeaks.py \
-M "${BIALLELICMACSFILE}" \
-O "${OUTDIR}" \
-P "${ALLELICPREFIX}" \
--keep_intermediates

printf "\n\nend: %s\n\n" "\$( date )"
EOF

done
done

chmod u+x biallelicPeakCounts.*.job



############
# Submit jobs
#sshgrid
while read JOBFILE; do
echo "${JOBFILE}"
qsub "${JOBFILE}"
done  < <(ls -1 biallelicPeakCounts.*.job)

# rm biallelicPeakCounts.*.job.o*


############
# Check

cat ${DATADIR}/allelicWork.V2/sciATAC_EB_MKII.alt/sciATAC_EB_MKII.alt.log
cat ${DATADIR}/allelicWork.V2/sciATAC_EB_MKII.ref/sciATAC_EB_MKII.ref.log

ls -ltrh ${DATADIR}/allelicWork.V2/sciATAC_EB_MKII.alt
ls -ltrh ${DATADIR}/allelicWork.V2/sciATAC_EB_MKII.ref

cat ${DATADIR}/allelicWork.V2/sciATAC_EB_MKII.alt/sciATAC_EB_MKII.alt.counts.biallelicPeaks.txt | awk 'END{print NF, NR}'
# 3 4,808,669
cat ${DATADIR}/allelicWork.V2/sciATAC_EB_MKII.ref/sciATAC_EB_MKII.ref.counts.biallelicPeaks.txt | awk 'END{print NF, NR}'
# 3 4,248,567

head -n 3 ${DATADIR}/allelicWork.V2/sciATAC_EB_MKII.alt/sciATAC_EB_MKII.alt.counts.biallelicPeaks.txt
# chr7_82648214_82648994  CTGAAGCTTAAGTCCTGAAGGCGAGAGCGTACTGAC    2
# chr2_66439727_66440927  TCCGGAGATGATATTGCGTACGATCATCTAATCTTA    1
# chr4_126428979_126430257        TAATGCGCAGTTCTCATGCTCTAACTCGTATAGCCT    1




##################################################
# Merge allelic counts in biallelic peaks

for (( i=0; i<${#libIDs[@]}; i++ )); do
libID=${libIDs[$i]}

ALLELICPREFIX="${libID}.allelic"
echo ${ALLELICPREFIX}

OUTDIR=${DATADIR}/allelicWork.V2/${ALLELICPREFIX}
mkdir -p ${OUTDIR}
# ln -sf ${OUTDIR} .

cp ../vijayWork/sci-ATAC_out/Plate1_mappings.txt ${OUTDIR}

ALTPREFIX="${libID}.alt"
ALTDIR=${DATADIR}/allelicWork.V2/${ALTPREFIX}

REFPREFIX="${libID}.ref"
REFDIR=${DATADIR}/allelicWork.V2/${REFPREFIX}

cat <( cat ${ALTDIR}/${ALTPREFIX}.counts.biallelicPeaks.txt | awk '{printf("%s\t%s_alt\t%d\n", $1, $2, $3)}' ) \
<( cat ${REFDIR}/${REFPREFIX}.counts.biallelicPeaks.txt | awk '{printf("%s\t%s_ref\t%d\n", $1, $2, $3)}' ) \
> ${OUTDIR}/${ALLELICPREFIX}.counts.biallelicPeaks.txt

done



############
# Check

ls -ltrh ${OUTDIR}/
cat ${OUTDIR}/${ALLELICPREFIX}.counts.biallelicPeaks.txt | awk 'END{print NF, NR}'
# 9,057,236






##################################################
##################################################
##################################################
# Get results

WDIR="2019_sciATAC_analysis/results/gbonora/20190208_sciATAC_mouseDiff"

# Images
rsync -auvPKLn --include='*.png' --include='*.pdf' --include='*/' --prune-empty-dirs  --exclude='*' gs:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/

# HTML and txt
rsync -auvPKLn --exclude="*counts*" --include='*.txt' --include='*.html' --include='*/' --prune-empty-dirs  --exclude='*' gs:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/

# Misc
rsync -auvPKLn --include='20180516_sciATAC_WTandDel1Patski_pilot_segregationResults.txt' --include='*/' --prune-empty-dirs  --exclude='*' gs:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/

# upload code, inc. R analysis
rsync -avuPKLn --include "*.txt" --include "*.py" --include "*.R" --include "*.sh" --include='*/' --prune-empty-dirs  --exclude='*' ~/NobleLab/proj/"${WDIR}"/ gs:proj/"${WDIR}"/
rsync -avuPKLn --exclude=".DS_Store" --include "*/src/*" --include='*/' --prune-empty-dirs  --exclude='*' ~/NobleLab/proj/"${WDIR}"/ gs:proj/"${WDIR}"/
rsync -avuPKLn --exclude=".DS_Store" --include "*/picard-tools-1.141/*" --include='*/' --prune-empty-dirs  --exclude='*' ~/NobleLab/proj/"${WDIR}"/ gs:proj/"${WDIR}"/
