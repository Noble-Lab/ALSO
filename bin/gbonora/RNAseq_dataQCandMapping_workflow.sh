##################################################
##################################################
##################################################
#
# 170214
# Giancarlo
# New Del1b RNA-seq libraries
#
# 170620
# New invDxz4 and WT RNA-seq libraries
#
##################################################

source ~/.bashrc
source ~/.bash_profile

PROJDIR="<projDirName>"
WORKDIR="<workDirName>"
DATADIR="<dataDirName>"
REFDIR="<refdataDirName>"


##################################################
# Working folder
mkdir -p "${PROJDIR}"/"${WORKDIR}"
cd  "${PROJDIR}"/"${WORKDIR}"






#################################################
#################################################
#################################################
# Get mm10 iGenomes data

mkdir "${REFDIR}"/iGenomes
cd "${REFDIR}"/iGenomes

wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
tar -xzvf Mus_musculus_UCSC_mm10.tar.gz

# cat fasta files
cat "${REFDIR}"/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/*.fa > ~/refdata/mm10/mm10.fa &








#################################################
#################################################
#################################################
# Tophat2 alignment
#
# http://ccb.jhu.edu/software/tophat/manual.shtml
# -T/--transcriptome-only			Only align the reads to the transcriptome and report only those mappings as genomic mappings.
#
# -x/--transcriptome-max-hits		Maximum number of mappings allowed for a read, when aligned to the transcriptome 
#									(any reads found with more then this number of mappings will be discarded).
#
# -M/--prefilter-multihits			When mapping reads on the transcriptome, some repetitive or low complexity reads that would be discarded 
#									in the context of the genome may appear to align to the transcript sequences and thus may end up reported as mapped 
#									to those genes only. This option directs TopHat to first align the reads to the whole genome in order to determine 
#									and exclude such multi-mapped reads (according to the value of the -g/--max-multihits option). 
#
# -g/--max-multihits <int>			Instructs TopHat to allow up to this many alignments to the reference for a given read, 
#									and choose the alignments based on their alignment scores if there are more than this number. 
#									The default is 20 for read mapping. Unless you use --report-secondary-alignments, 
#									TopHat will report the alignments with the best alignment score. 
#									If there are more alignments with the same score than this number, 
#									TopHat will randomly report only this many alignments. 
#									In case of using --report-secondary-alignments, TopHat will try to report alignments up to this option value, 
#									and TopHat may randomly output some of the alignments with the same score to meet this number.
#
# --report-secondary-alignments		By default TopHat reports best or primary alignments based on alignment scores (AS). 
#									Use this option if you want to output additional or secondary alignments (up to 20 alignments, 
#									this limit can be changed by using the -g/--max-multihits option above). 
#
# Will use default of reporting only primary alignments and not prefilter

WDDIR="${PROJDIR}"/"${WORKDIR}/tophatAlign"
mkdir "${WDDIR}"
cd "${WDDIR}"


DATADIR=$( cat "${PROJDIR}"/"${WORKDIR}"/DATADIR )

INPUTS=$( cat "${PROJDIR}"/"${WORKDIR}"/libIDs ) 

MMS=($( cat "${PROJDIR}"/"${WORKDIR}"/libMMs ))

INPUTDIR="${DATADIR}/fastq"

OUTPUTDIR="${DATADIR}/tophat2"
mkdir -p "${OUTPUTDIR}"

REFIDX="${REFDIR}/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"

GTF="${REFDIR}/iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"

PROCS=8

LIBTYPE="fr-unstranded"
#LIBTYPE="fr-firststrand"
#LIBTYPE="fr-secondstrand"


#################################################
# Generate jobs

TALLY=-1
for INPUT in ${INPUTS[@]}; do
TALLY=$(( TALLY + 1 ))
nMMs=${MMS[${TALLY}]}
BASEINPUT=$( basename "${INPUT}" ".fq.gz" )
echo "${INPUT}" "${BASEINPUT}" "${nMMs}"

# Triple the default to allow for SNPs (75bp reads)
# But be ware of this issue:
# Error: the read mismatches (6) and the read gap length (2) should be less than or equal to the read edit dist (2)
# Either decrease --read-mismatches or --read-gap-length, or increase --read-edit-dist
MISMATCHES="${nMMs}" # default 2
GAPLEN=2 # default 2
EDITDIST="${nMMs}" # default 2

OPTIONS="-p ${PROCS} --library-type ${LIBTYPE} \
--read-mismatches ${MISMATCHES} \
--read-gap-length ${GAPLEN} \
--read-edit-dist ${EDITDIST} \
--transcriptome-only \
--no-novel-juncs"

jobFile="${BASEINPUT}.tophat2.job"
cat << EOF > "${jobFile}"
#!/bin/bash -x
#\$ -pe serial ${PROCS} 
#\$ -l mfree=4G
#\$ -q noble-long.q 
#\$ -l h_rt=47:59:59
#\$ -cwd
#\$ -j y

source /etc/profile.d/modules.sh
module load modules modules-init modules-gs modules-noble
module load bowtie2/2.2.3
module load samtools/0.1.19
module load boost/1.52.0
module load tophat/2.0.12

which tophat

hostname

printf "\n\nstart: %s\n\n" "\$( date )"

tophat ${OPTIONS} -G "${GTF}" -o "${OUTPUTDIR}/${BASEINPUT}" "${REFIDX}" "${INPUTDIR}"/"${INPUT}"

printf "\n\nfinish: %s\n" "\$( date )"

EOF
chmod 755 "${jobFile}"
done


#################################################
# Submit jobs
#sshgrid

for INPUT in ${INPUTS[@]}; do
#for INPUT in PatskiDel1 PatskiDel2 PatskiDel5; do
echo "${INPUT}"
BASEINPUT=$( basename "${INPUT/.fq.gz/}" )

jobFile="${BASEINPUT}.tophat2.job"
qsub "${jobFile}"
done 


















#################################################
#################################################
#################################################
# Alignment summary

WDDIR="${PROJDIR}"/"${WORKDIR}/tophatAlign"
mkdir "${WDDIR}"
cd "${WDDIR}"

for INPUT in ${INPUTS[@]}; do
BASEINPUT=$( basename "${INPUT}" ".fq.gz" )
printf "\n${BASEINPUT}"
cat "${OUTPUTDIR}"/"${BASEINPUT}"/align_summary.txt 
done > align_summary.ALLsamples.txt


















#################################################
#################################################
#################################################
# Sort and index files

WDDIR="${PROJDIR}"/"${WORKDIR}/tophatAlign"
mkdir "${WDDIR}"
cd "${WDDIR}"

DATADIR=$( cat "${PROJDIR}"/"${WORKDIR}"/DATADIR )
INPUTS=$( cat "${PROJDIR}"/"${WORKDIR}"/libIDs ) 
INPUTDIR="${DATADIR}/tophat2"
OUTPUTDIR="${DATADIR}/sortedBAMs"
mkdir -p "${OUTPUTDIR}"

PROCS=8

#################################################
# Generate jobs

for INPUT in ${INPUTS[@]}; do
BASEINPUT=$( basename "${INPUT}" ".fq.gz" )
echo "${INPUT}" "${BASEINPUT}"

jobFile="${BASEINPUT}.sortBAMs.job"
cat << EOF > "${jobFile}"
#!/bin/bash -x
#\$ -pe serial ${PROCS} 
#\$ -l mfree=4G
#\$ -q noble-long.q 
#\$ -cwd
#\$ -j y


source /etc/profile.d/modules.sh
module load modules modules-init modules-gs modules-noble
module load samtools/0.1.19

hostname

printf "\n\nstart: %s\n\n" "\$( date )"

samtools sort -@ ${PROCS} -m 15G "${INPUTDIR}/${BASEINPUT}"/accepted_hits.bam "${OUTPUTDIR}"/"${BASEINPUT}".sorted

samtools index "${OUTPUTDIR}"/"${BASEINPUT}".sorted.bam

printf "\n\nfinish: %s\n" "\$( date )"

EOF
chmod 755 "${jobFile}"
done


#################################################
# Submit jobs
#sshgrid
for INPUT in ${INPUTS[@]}; do
echo "${INPUT}"
BASEINPUT=$( basename "${INPUT/.fq.gz/}" )

jobFile="${BASEINPUT}.sortBAMs.job"
qsub "${jobFile}"
done 


#################################################
# Check
ls -lrtrh "${OUTPUTDIR}"























#################################################
#################################################
#################################################
# Regions of interest

WDDIR="${PROJDIR}"/"${WORKDIR}/tophatAlign"
mkdir "${WDDIR}"
cd "${WDDIR}"

DATADIR=$( cat "${PROJDIR}"/"${WORKDIR}"/DATADIR )

INPUTS=$( cat "${PROJDIR}"/"${WORKDIR}"/libIDs ) 

INPUTDIR="${DATADIR}/sortedBAMs"

OUTPUTDIR="${DATADIR}/ROIBAMs"
mkdir -p "${OUTPUTDIR}"

ls -lrtrh "${INPUTDIR}"


# Load tools
module load samtools/1.3
#module load bedtools/2.24.0




#//# Del1: chrX:72882857-73010093 in MM9 = chrX:75637518-75764754 in MM10
#//# Inv1: chrX:72919385-72920291 in MM9 = chrX:75674046-75674952 in MM10
#//# Del2: chrX:72882839-72919375in MM9 = chrX:75637500-756740364 in MM10
#//# Del5: chrX:72966422-73010099 in MM9 = chrX:75721083-75764760 in MM10 

# http://www.artificialworlds.net/blog/2012/10/17/bash-associative-array-examples/
declare -A ROI
ROI[Del1]="chrX:75637518-75764754"
ROI[Del5]="chrX:75721083-75764760"
ROI[Del2]="chrX:75637500-756740364"
ROI[Xist]="chrX:103460373-103483233"
ROI[Firre]="chrX:50563120-50635321"
for K in "${!ROI[@]}"; do echo $K ${ROI[$K]}; done

offset=100000

for INPUT in ${INPUTS[@]}; do
INPUT="${INPUT/.fq.gz}"

# Loop through all keys in an associative array
for K in "${!ROI[@]}"; do 
echo $K ${ROI[$K]}

ADJROI=$( echo ${ROI[$K]} | cut -d':' -f1 )
ADJROI=$( printf "%s:%s" $ADJROI $( echo ${ROI[$K]} | cut -d':' -f2 | sed 's/-/\t/' | awk -v offset=$offset '{printf("%d-%d", $1-offset, $2+offset)}' ) )
echo $ADJROI

# Extract reads of interest
samtools view -bh "${INPUTDIR}"/"${INPUT}".sorted.bam $ADJROI > "${OUTPUTDIR}"/"${INPUT}".$K.sorted.bam

# ... & index
samtools index "${OUTPUTDIR}"/"${INPUT}".$K.sorted.bam "${OUTPUTDIR}"/"${INPUT}".$K.sorted.bai 

# Generate coverage track from chrX reads
# genomeCoverageBed -bg -trackline -trackopts name=${INPUT} -ibam ${INPUT}.$K.sorted.ba -g ${HOME}/refdata/mm10/chromInfo.txt | gzip > ${INPUT}.$K.sorted.bedgraph.gz &

done # ROIs

done # samples




ls -lrtrh "${OUTPUTDIR}"