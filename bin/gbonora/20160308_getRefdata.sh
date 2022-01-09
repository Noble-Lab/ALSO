##################################################################################################
##################################################################################################
##################################################################################################
# 160308
# Download refdata


source ~/.bashrc
source ~/.bash_profile


##################################################
##################################################
##################################################
# Move to new refdata folder
mkdir /net/noble/vol4/noble/user/gbonora/refdata
ln -s /net/noble/vol4/noble/user/gbonora/refdata $HOME




















##################################################################################################
##################################################################################################
##################################################################################################
# Mouse reference assemblies


##################################################
##################################################
##################################################
# Get mm9 iGenomes data

mkdir ~/refdata/iGenomes
cd ~/refdata/iGenomes

wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm9/Mus_musculus_UCSC_mm9.tar.gz
tar -xzvf Mus_musculus_UCSC_mm9.tar.gz





##################################################
##################################################
##################################################
# Get UCSC mm9 FASTA

mkdir -p ~/refdata/mm9/FASTA
cd ~/refdata/mm9/FASTA

wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/*'





##################################################
##################################################
##################################################
# Get  mm9 Bowtie2 index

mkdir -p ~/refdata/mm9/Bowtie2Index
cd ~/refdata/mm9/Bowtie2Index

wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm9.zip
unzip mm9.zip





##################################################
##################################################
##################################################
# Get UCSC mm9 chromSizes

mkdir -p ~/refdata/mm9
cd ~/refdata/mm9

wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/chromInfo.txt.gz'
gunzip chromInfo.txt.gz





##################################################
##################################################
##################################################
# Get UCSC mm10 chromSizes

mkdir -p ~/refdata/mm10
cd ~/refdata/mm10

wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/chromInfo.txt.gz'
gunzip chromInfo.txt.gz


###################
# 20180809
# 2bit
wget http://hgdownload.cse.ucsc.edu/gbdb/mm10/mm10.2bit .




















##################################################################################################
##################################################################################################
##################################################################################################
# 160331
# Download mappability tracks

mkdir -p ~/refdata/mm9/mappability
cd ~/refdata/mm9/mappability

# http://hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeMapability/md5sum.txt
# wgEncodeCrgMapabilityAlign100mer.bigWig 0fc0a6ff28f4e2ee780bf6702ed7aea7
# wgEncodeCrgMapabilityAlign36mer.bigWig cea40c2bca1572ebfab8edc8936ea9ab
# wgEncodeCrgMapabilityAlign40mer.bigWig 13cc3bf69325fb580586e79e75d82654
# wgEncodeCrgMapabilityAlign50mer.bigWig def4bb1db43b8d93f16cf91f6987ecad
# wgEncodeCrgMapabilityAlign75mer.bigWig 39b8af3492edd29d935f0b8d5b1dbdda

for Kmer in 36 40 50 75 100; do
# Kmer=36
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign${Kmer}mer.bigWig > /dev/null &
done

# md5
for Kmer in 36 40 50 75 100; do
printf "${Kmer}\t"
md5sum wgEncodeCrgMapabilityAlign${Kmer}mer.bigWig
done
#36      cea40c2bca1572ebfab8edc8936ea9ab  wgEncodeCrgMapabilityAlign36mer.bigWig
#40      13cc3bf69325fb580586e79e75d82654  wgEncodeCrgMapabilityAlign40mer.bigWig
#50      def4bb1db43b8d93f16cf91f6987ecad  wgEncodeCrgMapabilityAlign50mer.bigWig
#75      39b8af3492edd29d935f0b8d5b1dbdda  wgEncodeCrgMapabilityAlign75mer.bigWig
#100     0fc0a6ff28f4e2ee780bf6702ed7aea7  wgEncodeCrgMapabilityAlign100mer.bigWig


# Check it out
theChr="chrX"
theStart="100655712"
theEnd="100678572"
bigWigToBedGraph -chrom="${theChr}" -start="${theStart}" -end="${theEnd}" \
wgEncodeCrgMapabilityAlign100mer.bigWig \
wgEncodeCrgMapabilityAlign100mer."${theChr}"."${theStart}"."${theEnd}".bedGraph



##################################################
# 160405
# Mappable regions

# To bedGraph
for Kmer in 36 40 50 75 100; do
printf "${Kmer}\t"
bigWigToBedGraph wgEncodeCrgMapabilityAlign${Kmer}mer.bigWig wgEncodeCrgMapabilityAlign${Kmer}mer.bedGraph &
done

gzip *.bedGraph &

# Genome size
TOTbps=$( cat ../chromInfo.txt | grep -v "_random" | awk '{totbp+=$2} END{print totbp}' )
echo ${TOTbps}

# Mappable genome
printf "windowSize\tcoveredBps\ttotalBps\tcoverage%%\n" > mm9MappableGenome.txt
for Kmer in 36 40 50 75 100; do
# Kmer=36
printf "${Kmer}\n"
zcat wgEncodeCrgMapabilityAlign${Kmer}mer.bedGraph.gz | \
#awk '$4>0.5{print $0, $3-$2}'
#awk '$4==1{tallyA++}$4>0.5{tallyB++}END{print NR, tallyA, tallyA/NR, tallyB, tallyB/NR}' # NOTE: score > 0.5 ~ score == 1
awk -v Kmer="${Kmer}" -v totbp="${TOTbps}" '$4>0.5{mapbp+=($3-$2)} END{printf("%d\t%d\t%d\t%0.2f\n",Kmer,mapbp,totbp,mapbp/totbp*100)}' >> mm9MappableGenome.txt
done




















##################################################################################################
##################################################################################################
##################################################################################################
# 160415
# Rebuild BWA indices for new version bwa/0.7.3
# https://sourceforge.net/p/samtools/mailman/message/29193818/

WENXIUDIR="/net/noble/vol1/home/wenxiu"
OLDBWAIndexDir=${WENXIUDIR}/proj/blau/results/wenxiu/2012-08-15_bwa-index

ls -ltrh *.fa
# -rw-r--r-- 1 wenxiu noblelab 2.6G Mar 26  2013 spretus.fa
# -rw-r--r-- 1 wenxiu noblelab 2.6G Apr 29  2013 spretus.Fan.combined.fa
# -rw-r--r-- 1 wenxiu noblelab 2.6G Jul  8  2013 black6.fa

mkdir ~/refdata/patskiBWAindex
cd  ~/refdata/patskiBWAindex

# Just link files.
# cp -t . "${OLDBWAIndexDir}"/black6.fa  "${OLDBWAIndexDir}"/spretus.fa "${OLDBWAIndexDir}"/spretus.Fan.combined.fa
ln -sf "${OLDBWAIndexDir}"/black6.fa .
ln -sf "${OLDBWAIndexDir}"/spretus.fa .
ln -sf "${OLDBWAIndexDir}"/spretus.Fan.combined.fa .

module load bwa/0.7.3

while read FASTA; do
bwa index "${FASTA}" > "${FASTA/.fa/.bwaIndex.out}" 2>&1 &
done < <( ls -1 *.fa )

# *** NOTE: THREE INDIVIDUAL GENOMES INDEXED HERE! ***
# ***       NOT CONCATENATED.                      ***




















##################################################################################################
##################################################################################################
##################################################################################################
# 160415
# Bl6 genome from sanger
# ftp://ftp-mouse.sanger.ac.uk/current_denovo/README
#######################
# FILE DESCRIPTION    #
#######################
# <strain>.fa.gz : Assembled chromosomes+unplaced scaffolds (>2kbp)
# <strain>.fa.masked.gz : Repeatmasked chromsomes+unplaced scaffolds (>2kbp)
# <strain>.fa.out.gz : Repeatmasker output file of annotated repeats

# *** BUT NOT BL6. ACTUALLY, C57BL_6J IS BL6! ***

mkdir ~/refdata/C57BL_6NJ
cd ~/refdata/C57BL_6NJ

wget ftp://ftp-mouse.sanger.ac.uk/current_denovo/C57BL_6NJ.chromosomes.unplaced.gt2k.fa.gz
md5sum C57BL_6NJ.chromosomes.unplaced.gt2k.fa.gz
# 52c1d2b1d19327417e4de08e8c49d35f  C57BL_6NJ.chromosomes.unplaced.gt2k.fa.gz

# ftp://ftp-mouse.sanger.ac.uk/current_denovo/C57BL_6NJ.chromosomes.unplaced.gt2k.fa.gz.md5
# 52c1d2b1d19327417e4de08e8c49d35f  C57BL_6NJ/C57BL_6NJ.chromosomes.unplaced.gt2k.fa.gz

########
# 160421
# New indices for Sanger mm10 for BWA mem
gunzip C57BL_6NJ.chromosomes.unplaced.gt2k.fa.gz &
mv C57BL_6NJ.chromosomes.unplaced.gt2k.fa C57BL_6NJ.fa
module load bwa/0.7.3
bwa index "C57BL_6NJ.fa" > "C57BL_6NJ.bwaIndex.out" 2>&1 &





##################################################
##################################################
##################################################
# 160415
# Spretus genome from sanger
# ftp://ftp-mouse.sanger.ac.uk/current_denovo/README
#######################
# FILE DESCRIPTION    #
#######################
# <strain>.fa.gz : Assembled chromosomes+unplaced scaffolds (>2kbp)
# <strain>.fa.masked.gz : Repeatmasked chromsomes+unplaced scaffolds (>2kbp)
# <strain>.fa.out.gz : Repeatmasker output file of annotated repeats

mkdir ~/refdata/SPRET_EiJ
cd ~/refdata/SPRET_EiJ

wget ftp://ftp-mouse.sanger.ac.uk/current_denovo/SPRET_EiJ.chromosomes.unplaced.gt2k.fa.gz
md5sum SPRET_EiJ.chromosomes.unplaced.gt2k.fa.gz
# 7a907776b07085b0bd4fcce6af61fa53  SPRET_EiJ.chromosomes.unplaced.gt2k.fa.gz

# ftp://ftp-mouse.sanger.ac.uk/current_denovo/SPRET_EiJ.chromosomes.unplaced.gt2k.fa.gz.md5sum
# 7a907776b07085b0bd4fcce6af61fa53  SPRET_EiJ.chromosomes.unplaced.gt2k.fa.gz

########
# 160421
# New indices for BWA mem
gunzip SPRET_EiJ.chromosomes.unplaced.gt2k.fa.gz &
mv SPRET_EiJ.chromosomes.unplaced.gt2k.fa SPRET_EiJ.fa

module load bwa/0.7.3
bwa index "SPRET_EiJ.fa" > "SPRET_EiJ.bwaIndex.out" 2>&1 &





##################################################
##################################################
##################################################
# 160415
# Make a concatenated BL6 & Spretus genome from sanger

mkdir ~/refdata/C57BL_6NJ_SPRET_EiJ
cd ~/refdata/C57BL_6NJ_SPRET_EiJ

zcat ~/refdata/C57BL_6NJ/C57BL_6NJ.chromosomes.unplaced.gt2k.fa.gz | grep "^>" | more

cat <( zcat ~/refdata/C57BL_6NJ/C57BL_6NJ.chromosomes.unplaced.gt2k.fa.gz | sed 's/^\(>.*\)/\1_black6/' ) \
 <( zcat ~/refdata/SPRET_EiJ/SPRET_EiJ.chromosomes.unplaced.gt2k.fa.gz | sed 's/^\(>.*\)/\1_spretus/' ) \
 > ~/refdata/C57BL_6NJ_SPRET_EiJ/C57BL_6NJ_SPRET_EiJ.fa &

grep "^>"  ~/refdata/C57BL_6NJ_SPRET_EiJ/C57BL_6NJ_SPRET_EiJ.fa | more

#########################
# New indices for BWA mem
module load bwa/0.7.3
bwa index "C57BL_6NJ_SPRET_EiJ.fa" > "C57BL_6NJ_SPRET_EiJ.bwaIndex.out" 2>&1 &





##################################################
##################################################
##################################################
# 160503
# True Bl6 genome from sanger. See '# Reference genome and Ensembl version' in
# ftp://ftp-mouse.sanger.ac.uk/current_snps/README
# Also see:
# ftp://ftp-mouse.sanger.ac.uk/README
# # Reference genome and Ensembl version
# All SNP and indel calls are relative to the reference mouse genome C57BL/6J (GRCm38).
# The reference genome used for the alignment can be found here:
# ftp://ftp-mouse.sanger.ac.uk/ref/
#
# C57BL_6J vs C57BL_6NJ
# http://www.research.uci.edu/facilities-services/ular/ular-documents/jax_2015_know_your_b6_mouse.pdf

mkdir ~/refdata/C57BL_6J
cd ~/refdata/C57BL_6J

wget ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa

# 160531
# Add 'chr' to chrom IDs and remove annotation
cat GRCm38_68.fa | sed 's/^>\(.*\) dna:.*/>chr\1/' | sed 's/^>chrMT/>chrM/' | sed 's/^>chrJH/>JH/' | sed 's/^>chrGL/>GL/' > C57BL_6J.fa &

# check
grep "^>"  C57BL_6J.fa | more

###########
# BWA index
module load bwa/0.7.3
bwa index "C57BL_6J.fa" > "C57BL_6J.fa.bwaIndex.out" 2>&1 &





##################################################
##################################################
##################################################
# 160524
# Make a concatenated proper BL6 & Spretus genome from sanger
# I.e. Use C57BL_6J, not C57BL_6NJ

mkdir ~/refdata/C57BL_6J_SPRET_EiJ
cd ~/refdata/C57BL_6J_SPRET_EiJ

# Note above changes:
# cd ~/refdata/SPRET_EiJ/
# mv SPRET_EiJ.chromosomes.unplaced.gt2k.fa SPRET_EiJ.fa
# gzip SPRET_EiJ.fa

# zcat ~/refdata/C57BL_6J/GRCm38_68.fa.gz | grep "^>" | more
# zcat ~/refdata/C57BL_6J/GRCm38_68.fa.gz | sed 's/^\(>.*\) dna:.*/\1_black6/' |  grep "^>" | more
# zcat ~/refdata/SPRET_EiJ/SPRET_EiJ.fa.gz  | grep "^>" | more

# 160531
# Had to rerun this b/c failed to add 'chr' prefix to GRCm38_68/C57BL_6J chromosomes IDs
#cat <( cat ~/refdata/C57BL_6J/GRCm38_68.fa| sed 's/^\(>.*\) dna:.*/\1_black6/' )
cat <( cat ~/refdata/C57BL_6J/C57BL_6J.fa| sed 's/^\(>.*\)/\1_black6/' ) \
 <( cat ~/refdata/SPRET_EiJ/SPRET_EiJ.fa | sed 's/^\(>.*\)/\1_spretus/' ) \
 > ~/refdata/C57BL_6J_SPRET_EiJ/C57BL_6J_SPRET_EiJ.fa &

# check
grep "^>"  ~/refdata/C57BL_6J_SPRET_EiJ/C57BL_6J_SPRET_EiJ.fa | more

# BWA index
module load bwa/0.7.3
bwa index "C57BL_6J_SPRET_EiJ.fa" > "C57BL_6J_SPRET_EiJ.fa.bwaIndex.out" 2>&1 &




















##################################################################################################
##################################################################################################
##################################################################################################
# 160503
# Build mm10 pseudo spretus genome

# Based on "/net/noble/vol1/home/wenxiu/proj/blau/results/wenxiu/2012-08-14_spretus-genome"
# Specifically code in "/net/noble/vol1/home/wenxiu/proj/blau/results/wenxiu/2012-08-14_spretus-genome/runall.sh"
mkdir -p ~/proj/2010blau/results/wenxiu
cd ~/proj/2010blau/results/wenxiu
svn list -v file:///net/noble/vol1/svn/grch-projects/2010blau/results/wenxiu/2012-08-14_spretus-genome
svn co file:///net/noble/vol1/svn/grch-projects/2010blau/results/wenxiu/2012-08-14_spretus-genome 2012-08-14_spretus-genome
less ~/proj/2010blau/results/wenxiu/2012-08-14_spretus-genome/runall.sh

# Actually use "2010blau/results/wenxiu/2014-11-05_spretus-genome"
#
# On 3 May, 2016, at 6:21 PM, Wenxiu Ma wrote:
#
# I reconstructed the mm10 pseudo-genome in
# results/wenxiu/2014-11-05_spretus-genome/ and built the bwa index in
# 2014-11-06_bwa-index-mm10/. I used mgp.v4.* version from Sanger so I
# guess that version of scripts would work better for the latest Sanger
# mm10 SNPs.
#
mkdir -p ~/proj/2010blau/results/wenxiu
cd ~/proj/2010blau/results/wenxiu
svn list -v file:///net/noble/vol1/svn/grch-projects/2010blau/results/wenxiu/2014-11-05_spretus-genome
svn co file:///net/noble/vol1/svn/grch-projects/2010blau/results/wenxiu/2014-11-05_spretus-genome 2014-11-05_spretus-genome
less ~/proj/2010blau/results/wenxiu/2014-11-05_spretus-genome/runall.sh


# And "2010blau/results/wenxiu/2014-11-06_bwa-index-mm10/" to cat genomes and build index
mkdir -p ~/proj/2010blau/results/wenxiu
cd ~/proj/2010blau/results/wenxiu
svn list -v file:///net/noble/vol1/svn/grch-projects/2010blau/results/wenxiu/2014-11-06_bwa-index-mm10
svn co file:///net/noble/vol1/svn/grch-projects/2010blau/results/wenxiu/2014-11-06_bwa-index-mm10 2014-11-06_bwa-index-mm10
less ~/proj/2010blau/results/wenxiu/2014-11-06_bwa-index-mm10/runall.sh





##################################################
# New mm10 SNPs

mkdir $HOME/refdata/mm10pseudoSpretus
cd $HOME/refdata/mm10pseudoSpretus

cp -p ~/proj/2010blau/results/wenxiu/2014-11-05_spretus-genome/*.sh .
cp -p ~/proj/2010blau/results/wenxiu/2014-11-05_spretus-genome/*.py .

vi runall.sh # Add steps

# update to use new 'ftp://ftp-mouse.sanger.ac.uk/current_snps/' files
./runall.sh 1 > runall.step1.out 2>&1 &

# Compare new vcf to old one
zcat /net/noble/vol1/home/wenxiu/proj/blau/results/wenxiu/2012-08-14_spretus-genome/spretus.ATG.vcf.Fan.combined.gz | grep -v "^#" | wc -l
# 35,544,005
zcat /net/noble/vol1/home/wenxiu/proj/blau/results/wenxiu/2012-08-14_spretus-genome/spretus.ATG.vcf.Fan.combined.gz | grep "^X" | wc -l
# 1,634,281
zcat spretus.snps.vcf.gz | grep -v "^#" | wc -l
# 36,574,600
zcat spretus.snps.vcf.gz | grep "^X" | wc -l
# 1,367,287

# *** NOTE: No longer use Fan chrX SNPs ***


#####
# 160517
# BUG
#
# Traceback (most recent call last):
#   File "./summarize-spretus-snps.py", line 86, in <module>
#     snpChrFile.write("%s\t%d\t%.2f%%\t%d\t%.2f\n" % (k, chrDict[k], chrDict[k]*100.0/numSpretusSNP, chrsizeDict[k], chrsizeDict[k]*1.0/chrDict[k]))
# KeyError: 'MT'

cd $HOME/refdata/mm10pseudoSpretus

zcat spretus.snps.vcf.gz | grep -v "^##" | cut -f 1 | sort | uniq -c
# ...
# 1870163 7
# 1892910 8
# 1817612 9
#       1 #CHROM
#     380 MT
# 1367287 X
#    1615 Y

cat mm10.sizes
# ...
# chrUn_GL456396  21240
# chrUn_GL456368  20208
# chrM    16299
# chr4_JH584292_random    14945
# chr4_JH584295_random    1976

# Had to add this code to "summarize-spretus-snps.py":
#  # GB 160507
#  if words[0] == "MT":
#    words[0] = "M"

./runall.sh 1 > runall.step1.take2.out 2>&1 &


#####
# 160518

sshgrid
source ~/.bashrc
source ~/.bash_profile

cd $HOME/refdata/mm10pseudoSpretus

# From Wenxiu 160523:
# merge-snps-into-genome.py substitute snps to mm10/chr*.fa.gz and save
# the results to spretus/chr*.snp.fa.gz
# merge-indels-into-genome.py add indels to spretus/chr*.snp.fa.gz and
# save results to spretus/chr*.fa.gz

./runall.sh 2 > runall.step2.out 2>&1 &

#./runall.sh 3 > runall.step3.out 2>&1 &
#./runall.sh 4 > runall.step4.out 2>&1 &
#
# ./merge-snps-indels.py spretus.snps.vcf.gz spretus.indels.vcf.gz spretus.snps-indels.2.vcf.gz
# read 36574600 snps.
# Error! duplicate chr/pos variant 1/7076472
#
# See email exchange with Wenxiu about this.
#
# From Wenxiu 160523:
# The script merge-snps-indels.py is not useful. At the beginning I
# probably wanted to use that script to combined two vcfs and then make
# the pseudo-genome. But since the two vcfs are inconsistent so in the
# end I wrote two scripts merge-snps-into-genome.py and
# merge-indels-into-genome.py. Running them one after the other will
# make sure that indels overwrite snps.

# From Wenxiu:
# If you need snps substitution only, use spretus/chr*.snp.fa.gz



##################################################
# 160524
# BWA mem index

mkdir $HOME/refdata/mm10pseudoSpretusBWAindex
cd $HOME/refdata/mm10pseudoSpretusBWAindex

cp -p ~/proj/2010blau/results/wenxiu/2014-11-06_bwa-index-mm10/runall.sh .

vi runall.sh # Make edits

./runall.sh > runall.out 2>&1 &




















##################################################################################################
##################################################################################################
##################################################################################################
# Mouse transcriptomes


##################################################
##################################################
##################################################
# 160517
# Get mm9 transcriptome

mkdir ~/refdata/mm9/transcriptome
cd ~/refdata/mm9/transcriptome

# BUT THIS FILE JUST GIVES CO-ORDINATES. i.e. similar to a GTF
wget hgdownload.cse.ucsc.edu/goldenPath/mm9/database/transcriptome.txt.gz
# gunzip transcriptome.txt.gz

# RNA transcript RNA https://www.biostars.org/p/43594/
# ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/
# ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/
# ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.2/
# UCSC mm9 is based on 37.1

wget -O NCBI.build.37.1.rna.fa.gz ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/RNA/rna.fa.gz

zcat NCBI.build.37.1.rna.fa.gz | more



#####
# 160519
# Use iGenomes GTF file

GTF="${HOME}/refdata/iGenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"
DNAFA="${HOME}/refdata/iGenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa"
RNAFA="${HOME}/refdata/iGenomes/Mus_musculus/UCSC/mm9/Sequence/TranscriptomeFasta/transcriptome.fa"
mkdir $(dirname "${RNAFA}")

source /etc/profile.d/modules.sh
module load modules modules-init modules-gs modules-noble
module load bowtie2/2.2.3
module load samtools/0.1.19
module load boost/1.52.0
module load tophat/2.0.12

gtf_to_fasta "${GTF}" "${DNAFA}" "${RNAFA}"
ls -ltrh "${RNAFA}"
gzip "${RNAFA}"

ln -sf "${RNAFA}".gz ./iGenomes.mm9.rna.fa.gz





##################################################
##################################################
##################################################
# 160606
# Get mm10 iGenomes data

mkdir ~/refdata/iGenomes
cd ~/refdata/iGenomes

wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
tar -xzvf Mus_musculus_UCSC_mm10.tar.gz

# 180428
# cat fasta files
cat ~/refdata/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/*.fa > ~/refdata/mm10/mm10.fa &





##################################################
##################################################
##################################################
# 160606
# Get mm10 iGenomes ensemble data

cd ~/refdata/iGenomes

wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/Ensembl/GRCm38/Mus_musculus_Ensembl_GRCm38.tar.gz
tar -xzvf Mus_musculus_Ensembl_GRCm38.tar.gz




















##################################################################################################
##################################################################################################
##################################################################################################
# 160608
# Build mm10 color space index

cd /net/noble/vol2/home/gbonora/refdata/iGenomes/Mus_musculus/UCSC/mm10/Sequence

module load bowtie/1.0.0

mkdir BowtieColorIndex

bowtie-build --color BowtieIndex/genome.fa BowtieColorIndex/genomeC > genomeCbuild.out 2>&1 &

# To prevent Tophat from having to reconsistute BowtieColorIndex/genomeC.fa FASTA file from Bowtie index
cp -p BowtieIndex/genome.fa  BowtieColorIndex/genomeC.fa




















##################################################################################################
##################################################################################################
##################################################################################################
# 170412
# Refflat vs genes.gtf

cd ~/refdata/iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/

zcat  refFlat.txt.gz | cut -f1 | sort | uniq | wc -l
# 24421

cat genes.gtf | cut -f9 | sed 's/.*gene_id "\(.*\)"; gene_name ".*/\1/' | sort | uniq | wc -l
# 24421

zcat  refFlat.txt.gz | cut -f1 | sort | uniq > refFlat.txt.Genes.txt &
cat genes.gtf | cut -f9 | sed 's/.*gene_id "\(.*\)"; gene_name ".*/\1/' | sort | uniq > genes.gtf.Genes.txt &
diff refFlat.txt.Genes.txt genes.gtf.Genes.txt
# nada




















##################################################################################################
##################################################################################################
##################################################################################################
# 180209
# Get cast index
#
# On Jan 24, 2018, at 10:26 AM, Wenxiu Ma wrote:
#
# Hi Giancarlo,
#
# Thank you for willing to help. I have listed the relevant file paths below.
#
# We have two paired-end ATAC-seq samples, one for the control and one for CRISPR.
# http://giancarlo:patski@biocluster.ucr.edu/~wenxiu/proj/cdistechelab/X-inactivation/data/BL6xCast/2017-05-26_ES-ATAC-seq/
#
# The cells are male ES cells (BL6 X CAST). I downloaded the
# CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz file from Sanger and built a
# mm10/CAST pseudo-genome.
# http://giancarlo:patski@biocluster.ucr.edu/~wenxiu/proj/cdistechelab/X-inactivation//results/wenxiu/2016-06-27_CAST-genome-mm10/
#
# I used my allelic mapping pipeline and the results folder can be found at,
# http://giancarlo:patski@biocluster.ucr.edu/~wenxiu/proj/cdistechelab/X-inactivation/results/wenxiu/2017-06-18_BL6xCAST-bwa/
#
# The problem is that we observed mapping ratio between BL6 and CAST to
# be around 1.8 and 2.2 for CRISPR and control respectively. I could not
# figure out why the ratio is so biased (previously I got about 1.3
# ratio for RNA-seq samples in the same cells). Xinxian mentioned that
# you have a mapping pipeline that might help with the allelic biases.
# So it would be great if you could help run your allelic mapping
# pipeline on the same data.
#
# Thanks a lot!
#
# Wenxiu

mkdir $HOME/refdata/mm10pseudoCastBWAindex
cd $HOME/refdata/mm10pseudoCastBWAindex

# https://stackoverflow.com/questions/5043239/how-do-i-mirror-a-directory-with-wget-without-creating-parent-directories
wget -r --no-parent --cut-dirs=7 http://giancarlo:patski@biocluster.ucr.edu/~wenxiu/proj/cdistechelab/X-inactivation/results/wenxiu/2016-06-27_CAST-genome-mm10/ .

mv biocluster.ucr.edu/ 2016-06-27_CAST-genome-mm10/

ls -ltrh 2016-06-27_CAST-genome-mm10/
ls -ltrh 2016-06-27_CAST-genome-mm10/CAST

# *** THIS DID NOT WORK. INCORRECT BWA VERSION ***
# cp -p ~/proj/2010blau/results/wenxiu/2014-11-06_bwa-index-mm10/runall.sh .
# mv runall.sh runall.V1.sh

cp -p $HOME/refdata/mm10pseudoSpretusBWAindex/runall.sh .
vi runall.sh # Make edits... cast for spretus
./runall.sh > runall.out 2>&1 &




#################################################
#################################################
#################################################
# 180217
# Get updated cast index
#
# On Feb 15, 2018, at 4:16 PM, Wenxiu Ma wrote:
#
# Hi Giancarlo,
#
# I have updated the CAST pseudo-genome at
# http://biocluster.ucr.edu/~wenxiu/proj/cdistechelab/X-inactivation//results/wenxiu/2016-06-27_CAST-genome-mm10/CAST/.
# And the CAST.snps.vcf.gz SNP file is available at its parent folder. I
# will try my segregation script to check if it helps with the mapping
# biases.
# If possible, could you run your mapping pipeline again using the new
# CAST.snps.vcf.gz file?
#
# Thanks a lot!
#
# Wenxiu

mkdir $HOME/refdata/mm10pseudoCastBWAindexV2
cd $HOME/refdata/mm10pseudoCastBWAindexV2

# https://stackoverflow.com/questions/5043239/how-do-i-mirror-a-directory-with-wget-without-creating-parent-directories
wget -r --no-parent --cut-dirs=7 http://giancarlo:patski@biocluster.ucr.edu/~wenxiu/proj/cdistechelab/X-inactivation/results/wenxiu/2016-06-27_CAST-genome-mm10/ .

mv biocluster.ucr.edu/ 2016-06-27_CAST-genome-mm10/

ls -ltrh 2016-06-27_CAST-genome-mm10/
ls -ltrh 2016-06-27_CAST-genome-mm10/CAST

# *** THIS DID NOT WORK. INCORRECT BWA VERSION ***
# cp -p ~/proj/2010blau/results/wenxiu/2014-11-06_bwa-index-mm10/runall.sh .
# mv runall.sh runall.V1.sh

cp -p $HOME/refdata/mm10pseudoCastBWAindex/runall.sh .
vi runall.sh # Make edits... cast for spretus
./runall.sh > runall.out 2>&1 &




















##################################################################################################
##################################################################################################
##################################################################################################
# 180303
# Generate new spretus assembly for 4DN allelic test
# Just the Sanger SNPS: no Fan SNPs, no validation

mkdir ~/refdata/mm10pseudoSpretus.dbSNP142
cd ~/refdata/mm10pseudoSpretus.dbSNP142

mm10dir=/net/noble/vol1/data/ucsc/goldenPath/mm10/chromosomes/

#################################################
# Download variant files

# Indels
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/SPRET_EiJ.mgp.v5.indels.dbSNP142.normed.vcf.gz

# SNPs
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz


################
# Compare to previous spretus vcf (which was manually derived from merged vcf file (mgp.v5.merged.snps_all.dbSNP142.vcf.gz)
# Genome-wide

zcat ~/refdata/mm10pseudoSpretus/spretus.snps.vcf.gz | grep -v "^#" | grep -c PASS
# 36,574,600
zcat SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz | grep -v "^#" | grep -c PASS
# 41,667,479

# Indels
zcat SPRET_EiJ.mgp.v5.indels.dbSNP142.normed.vcf.gz | grep -v "^#" | grep -c PASS
# 5253579


################
# compare chrX

zcat $HOME/refdata/mm10pseudoSpretus/spretus.snps.vcf.gz | grep "^X" | wc -l
# 1,367,287
zcat SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz | grep -v "^#" | grep PASS | awk '$1=="X"{tally++} END{print tally}'
# 1,906,894

# Wow! Many more!
# I followed Wenxiu's approach of extracting spretus SNPs from the merged file.
# But, it seems that extracting SNPs from the merged file had issues!

# *** Compare to '/Users/ajs/Documents/Gian's Stuff/NobleLab/labDocWD/20170315_vSNPdensity/vSNPdensity.xlsx' ***


################
# Why different?

zcat ~/refdata/mm10pseudoSpretus/spretus.snps.vcf.gz | grep -v "^#" | grep PASS | head
# 1	3000185	rs585444580	G	T	999	PASS	AC=2;AN=2	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:114:34:0:285,114,0:255,102,0:2:57:34:0,0,32,2:0:-0.693132:.:1
# 1	3000234	rs579469519	G	A	999	PASS	AC=2;AN=2	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:47:0:288,154,0:255,141,0:2:54:47:0,0,41,6:0:-0.693147:.:1

zcat SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz | grep -v "^#" | grep PASS | head
# 1	3000023	.	C	A	94	PASS	DP4=0,0,9,0;DP=9;CSQ=A||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:33:9:0:140,33,0:121,27,0:2:23:9:0,0,9,0:0:-0.662043:.:1
# 1	3000185	rs585444580	G	T	228	PASS	DP4=0,0,32,2;DP=34;CSQ=T||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:114:34:0:285,114,0:255,102,0:2:57:34:0,0,32,2:0:-0.693132:.:1
# 1	3000234	rs579469519	G	A	228	PASS	DP4=0,0,41,6;DP=47;CSQ=A||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:47:0:288,154,0:255,141,0:2:54:47:0,0,41,6:0:-0.693147:.:1
# But qualityr (94) lower than other PASS records (228)
# And has no SNP ID.
# https://www.biostars.org/p/91893/

zcat SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz | grep -v "^#" | grep PASS | awk '$6==228{tally++} END{print tally}'
# 38,705,223

zcat SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz | grep -v "^#" | grep PASS | awk '$3~/rs/{tally++} END{print tally}'
# 39,956,443

zcat ~/refdata/mm10pseudoSpretus/spretus.snps.vcf.gz | grep -v "^#" | grep PASS | awk '$3~/rs/{tally++} END{print tally}'
# 35,355,946
# Not all previous spretus vcf file SNP have IDs either.

# Note in merged file the SNP does not have PASS Set!
zcat ~/refdata/mm10pseudoSpretus/mgp.v5.merged.snps_all.dbSNP142.vcf.gz | grep -v "^#" | grep PASS | awk '$3~/rs/{tally++} END{print tally}'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	129P2_OlaHsd	129S1_SvImJ	129S5SvEvBrd	AKR_J
# A_J	BALB_cJ	BTBR_T+_Itpr3tf_J	BUB_BnJ	C3H_HeH	C3H_HeJ	C57BL_10J	C57BL_6NJ	C57BR_cdJ	C57L_J	C58_J
# CAST_EiJ	CBA_J	DBA_1J	DBA_2J	FVB_NJ	I_LnJ	KK_HiJ	LEWES_EiJ	LP_J	MOLF_EiJ	NOD_ShiLtJ	NZB_B1NJ
# 	NZO_HlLtJ	NZW_LacJ	PWK_PhJ	RF_J	SEA_GnJ	SPRET_EiJ	ST_bJ	WSB_EiJ	ZALENDE_EiJ
# 1	3000023	.	C	A	153	MinDP;MinAB;Qual;Het	DP=170;DP4=2,0,168,0;CSQ=A||||intergenic_variant||||||||
# 	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:20:8:0:133,20,0:109,11,0:2:29:7:1,0,7,0:0:-0.636426:.:1	1/1:22:6
# :0.166667:152,22,0:137,18,0:2:36:6:0,0,6,0:0:-0.616816:.:1	1/1:11:4:0:70,11,0:51,4,0:2:24:3:1,0,3,0:0:-0.511536:.:0
# 1/1:15:4:0.25:80,15,0:68,12,0:2:25:4:0,0,4,0:0:-0.556411:.:0	1/1:11:3:0:72,11,0:63,9,0:2:34:3:0,0,3,0:0:-0.511536:.:0
# 0/1:3:1:0:40,3,3:37,3,0:2:40:1:0,0,1,0:0:-0.379885:.:0	1/1:36:10:0:194,36,0:174,30,0:2:38:10:0,0,10,0:0:-0.670168:.:1	1/1:22:6
# :0.333333:111,22,0:96,18,0:2:20:6:0,0,6,0:0:-0.616816:.:1	0/1:3:1:0:35,3,3:32,3,0:2:40:1:0,0,1,0:0:-0.379885:.:0	1/1:19:5
# :0:135,19,0:121,15,0:2:34:5:0,0,5,0:0:-0.590765:.:1	1/1:15:4:0:106,15,0:94,12,0:2:29:4:0,0,4,0:0:-0.556411:.:0	1/1:11:3
# :0:83,11,0:74,9,0:2:28:3:0,0,3,0:0:-0.511536:.:0	1/1:33:9:0:199,33,0:180,27,0:2:42:9:0,0,9,0:0:-0.662043:.:1	1/1:6:2:
# 0:56,6,0:50,6,0:2:27:2:0,0,2,0:0:-0.453602:.:0	1/1:15:4:0.5:75,15,0:63,12,0:2:25:4:0,0,4,0:0:-0.556411:.:0	1/1:15:4:0:79,15
# ,0:67,12,0:2:24:4:0,0,4,0:0:-0.556411:.:0	1/1:22:6:0:133,22,0:118,18,0:2:31:6:0,0,6,0:0:-0.616816:.:1	1/1:15:4:0:108,1
# 5,0:96,12,0:2:35:4:0,0,4,0:0:-0.556411:.:0	1/1:15:4:0.25:88,15,0:76,12,0:2:31:4:0,0,4,0:0:-0.556411:.:0	0/0:.:1:0:.,.,.:
# .,.,.:2:40:1:0,0,1,0:0:-0.379885:.:0	1/1:11:3:0:79,11,0:70,9,0:2:31:3:0,0,3,0:0:-0.511536:.:0	1/1:6:2:0.5:48,6,0:42,6,
# 0:2:20:2:0,0,2,0:0:-0.453602:.:0	1/1:15:4:0.5:70,15,0:58,12,0:2:22:4:0,0,4,0:0:-0.556411:.:0	1/1:6:2:0:69,6,0:63,6,0:
# 2:40:2:0,0,2,0:0:-0.453602:.:0	1/1:19:5:0.2:94,19,0:80,15,0:2:24:5:0,0,5,0:0:-0.590765:.:1	1/1:15:4:0:87,15,0:75,12,0:2:30:
# 4:0,0,4,0:0:-0.556411:.:0	1/1:26:7:0.142857:150,26,0:134,21,0:2:31:7:0,0,7,0:0:-0.636426:.:1	1/1:19:5:0:132,19,0:118,
# 15,0:2:48:5:0,0,5,0:0:-0.590765:.:1	1/1:30:8:0:171,30,0:153,24,0:2:37:8:0,0,8,0:0:-0.651104:.:1	1/1:11:3:0.333333:72,11,
# 0:63,9,0:2:33:3:0,0,3,0:0:-0.511536:.:0	1/1:19:5:0:137,19,0:123,15,0:2:38:5:0,0,5,0:0:-0.590765:.:1	1/1:22:6:0.166667:137,22
# ,0:122,18,0:2:42:6:0,0,6,0:0:-0.616816:.:1	1/1:33:9:0:140,33,0:121,27,0:2:23:9:0,0,9,0:0:-0.662043:.:1	1/1:33:9:0.11111
# 1:170,33,0:151,27,0:2:29:9:0,0,9,0:0:-0.662043:.:1	1/1:26:7:0:140,26,0:124,21,0:2:31:7:0,0,7,0:0:-0.636426:.:1	1/1:6:2:
# 0:48,6,0:42,6,0:2:21:2:0,0,2,0:0:-0.453602:.:0


################
# Different question: Any randomn chrs in vcf?

#zcat SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz | grep -v "^#" | cut -f 1 | sort | uniq
zcat SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz | grep -v "^#" | awk 'NR==1{prevchr=""}{if ($1!=prevchr) {print $1}; prevchr=$1}'
# 1
# 2
# 3
# 4
# 5
# 6
# 7
# 8
# 9
# 10
# 11
# 12
# 13
# 14
# 15
# 16
# 17
# 18
# 19
# X
# Y


# which chroms?
zcat SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz  | grep -v "^#" | grep PASS | cut -f 1 | sort | uniq -c | awk '{printf("chr%s\t%d\n", $2, $1)}' | sort -k1V
# chr1    3279072
# chr2    2844052
# chr3    2699289
# chr4    2466685
# chr5    2461122
# chr6    2432543
# chr7    2111798
# chr8    2129225
# chr9    1999558
# chr10   2218545
# chr11   1943814
# chr12   1860064
# chr13   1927204
# chr14   1959713
# chr15   1756001
# chr16   1650396
# chr17   1533446
# chr18   1489615
# chr19   976513
# chrX    1906894
# chrY    21930




#################################################
# Generate new pseudo spretus

# cp -p ~/proj/2010blau/results/wenxiu/2014-11-05_spretus-genome/*.sh .
# cp -p ~/proj/2010blau/results/wenxiu/2014-11-05_spretus-genome/*.py .
cp -p  $HOME/refdata/mm10pseudoSpretus/*.sh .
cp -p  $HOME/refdata/mm10pseudoSpretus/*.py .

vi runall.sh # No need to extract SNPs from merged file i.e. no step 1

vi merge-snps-into-genome.py # Add code to filter SNPs that was in extract-spretus-snps.py

./runall.sh

############
# Check chr1

more spretus.chr1.job.o50644479
# n014.grid.gs.washington.edu
# Read chromosome 1 of length 195471971.
# Changed 3279072 bases.

zcat SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz | grep -v "^#" | grep PASS | awk '$1==1{tally++} END{print tally}'
# 3279072
# Yay!


############
# Check chrX, Y and M

more spretus.chrX.job.o50644510
# n014.grid.gs.washington.edu
# Read chromosome X of length 171031299.
# Changed 1906894 bases.

more spretus.chrY.job.o50644511
# n004.grid.gs.washington.edu
# Read chromosome Y of length 91744698.
# Changed 21930 bases.

more spretus.chrM.job.o50644513
# n014.grid.gs.washington.edu
# Read chromosome M of length 16299.
# Changed 0 bases.



##################################################
# 160524
# BWA mem index
# mm10dir=/net/noble/vol1/data/ucsc/goldenPath/mm10/chromosomes/

mkdir $HOME/refdata/mm10pseudoSpretus.dbSNP142.bwa.7.17-r1188.index
cd $HOME/refdata/mm10pseudoSpretus.dbSNP142.bwa.0.7.17-r1188.index

# cp -p ~/proj/2010blau/results/wenxiu/2014-11-06_bwa-index-mm10/runall.sh .
cp -p $HOME/refdata/mm10pseudoSpretusBWAindex/runall.sh .

vi runall.sh # Make edits

./runall.sh > runall.out 2>&1 &




##################################################
# 180309
# Generate PASS-filtered VCF?

zcat $HOME/refdata/mm10pseudoSpretus/spretus.patski.snps.vcf.gz | grep -v PASS

grep "SNPs from" /net/noble/vol2/home/gbonora/proj/2018_4DNomics_analysis/results/gbonora/20180307_4DNomicsAllelicHiCpipelineTest_bwamem5M_Patski/segregation/insituDNaseHiC.WG.patski.WT.rep3.H3W5NAFXX.L1_1.segregation.job.o50897960
# Read 45340081 SNPs from /tmp/50897960.1.noble-long.q/SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz.

zcat $HOME/refdata/mm10pseudoSpretus.dbSNP142/SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz | grep -vc "^#"
# 45,340,081
# CRAP!

zcat $HOME/refdata/mm10pseudoSpretus.dbSNP142/SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz | awk '$7=="PASS"{tally++}END{print tally}'
# 41,667,479

# *** RATHER THAN FILTER VCF, EDIT 'segregate-mapped-reads-CIGARsed.py' TO CHECK FOR THIS! ***
#     # GB 180309
#     # Filter loaded SNPs
#     if words[6] != 'PASS':
#         continue




















##################################################################################################
##################################################################################################
##################################################################################################
# 180308
# Generate new hg38 assembly for 4DN allelic test
# hg38 allelic references
# https://docs.google.com/document/d/11OjVk26qlZguoGXD4NZWS61ezFEyJ99-L66LsVGYPoQ/edit#
# https://docs.google.com/document/d/1770GqtfY7YO1naGG2NqZAraEUJpyNw-iRMMJECxVnJ8/edit
# ftp://platgene_ro:''@ussd-ftp.illumina.com/2017-1.0/hg38/hybrid/hg38.hybrid.vcf.gz
# https://github.com/Illumina/PlatinumGenomes/blob/master/files/2017-1.0.files
#
# Useful links:
# https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it
# https://www.biostars.org/p/119731/
#

mkdir ~/refdata/hg38.phased.NA12878
cd ~/refdata/hg38.phased.NA12878

hg38dir=/net/noble/vol1/data/ucsc/goldenPath/hg38/chroms/


#################################################
# Get VCF file

mkdir ~/refdata/hg38.phased.NA12878/vcf
cd ~/refdata/hg38.phased.NA12878/vcf

wget ftp://platgene_ro@ussd-ftp.illumina.com/2017-1.0/hg38/hybrid/hg38.hybrid.vcf.gz
# Login incorrect.

# https://github.com/Illumina/PlatinumGenomes
wget ftp://platgene_ro:''@ussd-ftp.illumina.com/2017-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz
# That worked!

# Try:
wget ftp://platgene_ro:''@ussd-ftp.illumina.com/2017-1.0/hg38/hybrid/hg38.hybrid.vcf.gz
# Worked too!

# So do files from https://github.com/Illumina/PlatinumGenomes/blob/master/files/2017-1.0.files
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/hybrid/README.md
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/hybrid/hg38.hybrid.vcf.gz # <--------- *** USE THIS ONE ***
wget https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz

#unzip *.gz

cat hg38.hybrid.vcf | grep "^chr1" | head -n 10
# chr1    727477  .       G       A       0       PASS    SOURCE=PG       GT      1|1
# chr1    727717  .       G       C       0       PASS    SOURCE=PG       GT      1|1
# chr1    740738  .       C       T       0       PASS    SOURCE=PG       GT      1|0
# chr1    744224  .       C       A       0       PASS    SOURCE=PG       GT      1|1
# chr1    766566  .       A       G       0       PASS    SOURCE=PG       GT      1|1
# chr1    768987  .       T       C       0       PASS    SOURCE=PG       GT      1|1
# chr1    770502  .       G       A       0       PASS    SOURCE=PG       GT      1|0
# chr1    770988  .       A       G       0       PASS    SOURCE=PG       GT      1|1
# chr1    771398  .       G       A       0       PASS    SOURCE=PG       GT      1|1
# chr1    772711  .       A       G       0       PASS    SOURCE=PG       GT      1|1

for LOC in $( cat hg38.hybrid.vcf | grep "^chr1" | head -n 10 | cut -f 2); do
#LOC=5275059
cat ${hg38dir}/chr1.fa | awk -v loc=${LOC} '$0!~/^>/{blocklen=length($0); totlen+=blocklen; blocks++;
if(totlen>=loc){
blockloc=loc-(blocklen*(blocks-1));
print((blocklen*(blocks-1))+blockloc, substr($0, blockloc, 1), blocklen, totlen);
exit }
}'
done
# 727477 G 50 727500
# 727717 G 50 727750
# 740738 C 50 740750
# 744224 C 50 744250
# 766566 a 50 766600
# 768987 t 50 769000
# 770502 g 50 770550
# 770988 A 50 771000
# 771398 g 50 771400
# 772711 A 50 772750

# All pass filter!
cat hg38.hybrid.vcf | grep -v "^#" | cut -f 7 | grep -v "PASS" | wc -l
# 0

# filter for SNPs only
cat hg38.hybrid.vcf | awk '$0~/^#/{print $0} length($4)==1 && length($5)==1{print $0}' > hg38.hybrid.snps.vcf



#################################################
# Generate consensus genomes for each haplotype
#
# See '20160202_getApps.sh' to download and install bcftools
#

cd ~/refdata/hg38.phased.NA12878

# Concatenate the hg38 genome
hg38dir=/net/noble/vol1/data/ucsc/goldenPath/hg38/chroms/
cat $hg38dir/chr[0-9]{,?}.fa $hg38dir/chrX.fa $hg38dir/chrY.fa $hg38dir/chrM.fa | gzip > hg38.fa.gz

# bgzip compress vcf
bgzip -c vcf/hg38.hybrid.snps.vcf > vcf/hg38.hybrid.snps.vcf.bgz
tabix -p vcf vcf/hg38.hybrid.snps.vcf.bgz

# incorporate na12878 SNPs/INDELS into hg38 and create hg38 -> na12878 chain files
# bcftools consensus -f hg38.fa.gz -H 1 vcf/hg38.hybrid.vcf.bgz -i  | gzip > hg38_na12878_hap1.fa.gz &
bcftools consensus -f hg38.fa.gz -H 1 vcf/hg38.hybrid.snps.vcf.bgz | gzip > hg38_na12878_hap1.fa.gz &
bcftools consensus -f hg38.fa.gz -H 2 vcf/hg38.hybrid.snps.vcf.bgz | gzip > hg38_na12878_hap2.fa.gz &


# check:
cat vcf/hg38.hybrid.snps.vcf | grep "^chr1" | head -n 10
# chr1    727477  .       G       A       0       PASS    SOURCE=PG       GT      1|1
# chr1    727717  .       G       C       0       PASS    SOURCE=PG       GT      1|1
# chr1    740738  .       C       T       0       PASS    SOURCE=PG       GT      1|0
# chr1    744224  .       C       A       0       PASS    SOURCE=PG       GT      1|1
# chr1    766566  .       A       G       0       PASS    SOURCE=PG       GT      1|1
# chr1    768987  .       T       C       0       PASS    SOURCE=PG       GT      1|1
# chr1    770502  .       G       A       0       PASS    SOURCE=PG       GT      1|0
# chr1    770988  .       A       G       0       PASS    SOURCE=PG       GT      1|1
# chr1    771398  .       G       A       0       PASS    SOURCE=PG       GT      1|1
# chr1    772711  .       A       G       0       PASS    SOURCE=PG       GT      1|1


# Check hap1:
for LOC in $( cat vcf/hg38.hybrid.snps.vcf | grep "^chr1" | head -n 10 | cut -f 2); do
#LOC=792149
zcat hg38_na12878_hap1.fa.gz  | awk -v loc=${LOC} '$0!~/^>/{blocklen=length($0); totlen+=blocklen; blocks++;
if(totlen>=loc){
blockloc=loc-(blocklen*(blocks-1));
print((blocklen*(blocks-1))+blockloc, substr($0, blockloc, 1), blocklen, totlen);
exit }
}'
done
# 727477 A 60 727500
# 727717 C 60 727740
# 740738 T 60 740760
# 744224 A 60 744240
# 766566 G 60 766620
# 768987 C 60 769020
# 770502 A 60 770520
# 770988 G 60 771000
# 771398 A 60 771420
# 772711 G 60 772740


# Check hap2:
for LOC in $( cat vcf/hg38.hybrid.snps.vcf | grep "^chr1" | head -n 10 | cut -f 2); do
#LOC=792149
zcat hg38_na12878_hap2.fa.gz | awk -v loc=${LOC} '$0!~/^>/{blocklen=length($0); totlen+=blocklen; blocks++;
if(totlen>=loc){
blockloc=loc-(blocklen*(blocks-1));
print((blocklen*(blocks-1))+blockloc, substr($0, blockloc, 1), blocklen, totlen);
exit }
}'
done
# 727477 A 60 727500
# 727717 C 60 727740
# 740738 C 60 740760
# 744224 A 60 744240
# 766566 G 60 766620
# 768987 C 60 769020
# 770502 g 60 770520
# 770988 G 60 771000
# 771398 A 60 771420
# 772711 G 60 772740

##########
# Examples:

# chr1 792149     .       A       G       0       PASS    SOURCE=PG       GT      0|1
# hap1 A
# hap2 G

# chr1 727477     .       G       A       0       PASS    SOURCE=PG       GT      1|1
# hap1 A
# hap2 A
# *** THESE DON'T SEEM USEFUL FOR SEGREGATION? ***

# chr1    740738  .       C       T       0       PASS    SOURCE=PG       GT      1|0
# hap1 T
# hap2 C

# *** NO "0|0" SNPs ***



#################################################
# Generate 'phased' vcf file for each haplotype -- necessary for segregation
#
# THIS IS NOT NECESSARY AS READS AND VCF BASES ARE BOTH IN UPPER CASE
# # *** -------------- ALSO ENSURE ALL BASES ARE IN UPPER CASE ------------------ ***
# # *** Segregation script does not consider case of SNPs and reads are uppercase ***
#
# 180311
# *** And must remove "chr" from chrom field to be compatible with segregation script! ***
#

cd ~/refdata/hg38.phased.NA12878/vcf

# Check out mouse vcf
zcat $HOME/refdata/mm10pseudoSpretus.dbSNP142/SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz | grep -v "^#" | head
# #CHROM  POS     ID              REF     ALT     QUAL    FILTER  INFO    FORMAT  SPRET_EiJ
# 1    3000023    .               C       A       94      PASS    DP4=0,0,9,0;DP=9;CSQ=A||||intergenic_variant||||||||    GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI  1/1:33:9:0:140,33,0:121,27,0:2:23:9:0,0,9,0:0:-0.662043:.:1
# 1    3000126    rs580370473     G       T       31.2311 Het     DP4=4,0,7,0;DP=11;CSQ=T||||intergenic_variant||||||||   GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI  0/1:20:11:0.0909091:72,0,20:65,0,21:2:45:7:4,0,7,0:0:-0.636426:.:0
# 1    3000185    rs585444580     G       T       228     PASS    DP4=0,0,32,2;DP=34;CSQ=T||||intergenic_variant||||||||  GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI  1/1:114:34:0:285,114,0:255,102,0:2:57:34:0,0,32,2:0:-0.693132:.:1
# 1    3000234    rs579469519     G       A       228     PASS    DP4=0,0,41,6;DP=47;CSQ=A||||intergenic_variant||||||||  GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI  1/1:127:47:0:288,154,0:255,141,0:2:54:47:0,0,41,6:0:-0.693147:.:1

cat hg38.hybrid.snps.vcf | awk 'BEGIN{FS="\t"; OFS="\t"}; \
$0~/^#/{print $0}; \
$10=="0|1"{sub("chr", "", $1); print $1, $2, $3, $4, $5, $6, $7, $8, $9, "0|1"}; \
$10=="1|1"{sub("chr", "", $1); print $1, $2, $3, $5, $5, $6, $7, $8, $9, "0|1"}; \
$10=="1|0"{sub("chr", "", $1); print $1, $2, $3, $5, $4, $6, $7, $8, $9, "0|1"}' > hg38.hybrid.snps.phased.vcf

# But 1|1 genotype not informative so leave out
cat hg38.hybrid.snps.vcf | awk 'BEGIN{FS="\t"; OFS="\t"}; \
$0~/^#/{print $0}; \
$10=="0|1"{sub("chr", "", $1); print $1, $2, $3, $4, $5, $6, $7, $8, $9, "0|1"}; \
$10=="1|0"{sub("chr", "", $1); print $1, $2, $3, $5, $4, $6, $7, $8, $9, "0|1"}' > hg38.hybrid.snps.filtered.phased.vcf

cat hg38.hybrid.vcf | grep -v "^#" | wc -l
# 4,251,128

cat hg38.hybrid.snps.vcf | grep -v "^#" | wc -l
# 3,619,127

cat hg38.hybrid.snps.phased.vcf | grep -v "^#" | wc -l
# 3,619,171

cat hg38.hybrid.snps.filtered.phased.vcf | grep -v "^#"  | wc -l
# 2,188,782

cat hg38.hybrid.snps.filtered.phased.vcf  | gzip > hg38.hybrid.snps.filtered.phased.vcf.gz

# Check:
zcat hg38.hybrid.snps.filtered.phased.vcf.gz  | grep -v "^#" | head


#################################################
# Compare to snps to indels
# Remember all pass filter.

# filter for informative snps only
cat hg38.hybrid.vcf | grep -v "^#" | awk '$10!="1|1" && (length($4)==1 && length($5)==1){tally++} END{print tally}'
# 2,188,782

# filter for informative indels only
cat hg38.hybrid.vcf | grep -v "^#" | awk '$10!="1|1" && (length($4)>1 || length($5)>1){tally++} END{print tally}'
# 632,001


#################################################
# which chroms?

# snps
cat hg38.hybrid.vcf | grep -v "^#" | awk '$10!="1|1" && (length($4)==1 && length($5)==1){print $0}' | cut -f 1 | sort | uniq -c
#  158602 chr1
#  111886 chr10
#  104830 chr11
#   99880 chr12
#   77334 chr13
#   68823 chr14
#   59809 chr15
#   67190 chr16
#   54595 chr17
#   61852 chr18
#   51285 chr19
#  172660 chr2
#   50087 chr20
#   30687 chr21
#   31573 chr22
#  147607 chr3
#  151899 chr4
#  140412 chr5
#  148742 chr6
#  121783 chr7
#  116040 chr8
#   95069 chr9
#   66137 chrX

# indels
cat hg38.hybrid.vcf | grep -v "^#" | awk '$10!="1|1" && (length($4)>1 || length($5)>1){print $0}' | cut -f 1 | sort | uniq -c
#   29738 chr1
#   19823 chr10
#   18275 chr11
#   19272 chr12
#   14301 chr13
#   12815 chr14
#   10626 chr15
#   10573 chr16
#   10754 chr17
#   10973 chr18
#   10680 chr19
#   31432 chr2
#    8873 chr20
#    5402 chr21
#    5625 chr22
#   27024 chr3
#   26829 chr4
#   25304 chr5
#   26069 chr6
#   22168 chr7
#   19343 chr8
#   16224 chr9
#   13384 chrX




#################################################
#################################################
# Build bwa indices
for genome in hap1 hap2; do
    ~/bin/bwa index hg38_na12878_$genome.fa.gz
done




















##################################################################################################
##################################################################################################
##################################################################################################
# 180314
# Compare to snps to indels in mm10 and hg38
# Remember all pass filter in hg39 but not in mm10 vcf

mkdir ~/refdata/mm10hg38vcfComparison
cd ~/refdata/mm10hg38vcfComparison

mm10vcfDir="$HOME/refdata/mm10pseudoSpretus.dbSNP142"
hg39vcfDir="$HOME/refdata/hg38.phased.NA12878/vcf"

for assembly in mm10 hg38; do
$HOME/bin/UCSC.utilities/fetchChromSizes $assembly | grep -v "_" | sort -k1V | grep -v "chrM" > $assembly.sizes
done


#####
# Totals

mm10snps=$(zcat ${mm10vcfDir}/SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz | grep -v "^#" | grep -c PASS)
# 41,667,479
mm10indels=$(zcat ${mm10vcfDir}/SPRET_EiJ.mgp.v5.indels.dbSNP142.normed.vcf.gz | grep -v "^#" | grep -c PASS)
# 5,253,579

# filter for informative snps only
hg39snps=$(cat ${hg39vcfDir}/hg38.hybrid.vcf | grep -v "^#" | awk '$10!="1|1" && (length($4)==1 && length($5)==1){tally++} END{print tally}')
# 2,188,782

# filter for informative indels only
hg39indels=$(cat ${hg39vcfDir}/hg38.hybrid.vcf | grep -v "^#" | awk '$10!="1|1" && (length($4)>1 || length($5)>1){tally++} END{print tally}')
# 632,001


#####
# Per chrom

zcat ${mm10vcfDir}/SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz  | grep -v "^#" | grep PASS | cut -f 1 | sort | uniq -c | awk '{printf("chr%s\t%d\n", $2, $1)}' | sort -k1V > mm10.snpsPerChr.txt

zcat ${mm10vcfDir}/SPRET_EiJ.mgp.v5.indels.dbSNP142.normed.vcf.gz | grep -v "^#" | grep PASS | cut -f 1 | sort | uniq -c | awk '{printf("chr%s\t%d\n", $2, $1)}' | sort -k1V > mm10.indelsPerChr.txt

cat ${hg39vcfDir}/hg38.hybrid.vcf | grep -v "^#" | awk '$10!="1|1" && (length($4)==1 && length($5)==1){print $0}' | cut -f 1 | sort | uniq -c | awk 'BEGIN{OFS="\t"} {print $2, $1}' | sort -k1V > hg38.snpsPerChr.txt
printf "chrY\t0\n" >> hg38.snpsPerChr.txt

cat ${hg39vcfDir}/hg38.hybrid.vcf | grep -v "^#" | awk '$10!="1|1" && (length($4)>1 || length($5)>1){print $0}' | cut -f 1 | sort | uniq -c | awk 'BEGIN{OFS="\t"} {print $2, $1}' | sort -k1V > hg38.indelsPerChr.txt
printf "chrY\t0\n" >> hg38.indelsPerChr.txt

more *.txt


#####
# Aggregate & calc densities

for assembly in mm10 hg38; do

echo $assembly

paste <( cat $assembly.sizes ) <( cut -f 2 $assembly.snpsPerChr.txt ) <( cut -f 2 $assembly.indelsPerChr.txt ) > $assembly.snpVindels.tmp

printf "chr\tlength\tSNPs\tbasesPerSNP\tIndels\tbasesPerIndel\tSNP2IndelRatio\n" > $assembly.snpVindels.tsv

cat $assembly.snpVindels.tmp | awk 'BEGIN{FS="\t"; OFS="\t"}
{
if($3!=0){snpDensity=$2/$3; snpDensity=sprintf("%0.2f", snpDensity)} else {snpDensity="nan"};
if($4!=0){indelDensity=$2/$4; indelDensity=sprintf("%0.2f", indelDensity)} else {indelDensity="nan"};
if($4!=0){snp2indelRatio=$3/$4; snp2indelRatio=sprintf("%0.2f", snp2indelRatio)} else {snp2indelRatio="nan"};
printf("%s\t%d\t%d\t%0.2f\t%d\t%0.2f\t%s\n", $1, $2, $3, snpDensity, $4, indelDensity, snp2indelRatio);
totBases+=$2; totSNPs+=$3; totIndels+=$4
}
END{
snpDensity=totBases/totSNPs; indelDensity=totBases/totIndels; snp2indelRatio=totSNPs/totIndels;
printf("%s\t%d\t%d\t%0.2f\t%d\t%0.2f\t%0.2f\n", "genome", totBases, totSNPs, snpDensity, totIndels, indelDensity, snp2indelRatio);
}' >> $assembly.snpVindels.tsv

rm $assembly.snpVindels.tmp

# HTML
/net/noble/vol1/home/noble/bin/tsv2html.py "$assembly.snpVindels.tsv" > "$assembly.snpVindels.html"

done

more *.tsv










#################################################
#################################################
#################################################
# 180314
# Compare to snps to indels in different mm10 vcfs

mkdir ~/refdata/mm10vcfComparison
cd ~/refdata/mm10vcfComparison

for assembly in mm10; do
$HOME/bin/UCSC.utilities/fetchChromSizes $assembly | grep -v "_" | sort -k1V | grep -v "chrM" > $assembly.sizes
done


#####
#####
# Per chrom SNPs

zcat $HOME/refdata/mm10pseudoSpretus.dbSNP142/SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf.gz  | grep -v "^#" | grep PASS | cut -f 1 | sort | uniq -c | awk '{printf("chr%s\t%d\n", $2, $1)}' | sort -k1V > spretus.dbSNP142.snpsPerChr.txt

zcat $HOME/refdata/mm10pseudoSpretus/spretus.snps.vcf.gz | grep -v "^#" | grep PASS | cut -f 1 | sort | uniq -c | awk '{printf("chr%s\t%d\n", $2, $1)}' | sort -k1V > spretus.extracted.snpsPerChr.txt
grep -v chrMT spretus.extracted.snpsPerChr.txt > spretus.extracted.snpsPerChr.txt.tmp
mv spretus.extracted.snpsPerChr.txt.tmp spretus.extracted.snpsPerChr.txt

zcat $HOME/refdata/mm10pseudoSpretus/spretus.patski.snps.vcf.gz | grep -v "^#" | grep PASS | cut -f 1 | sort | uniq -c | awk '{printf("chr%s\t%d\n", $2, $1)}' | sort -k1V > patski.validated.snpsPerChr.txt
printf "chrY\t0\n" >> patski.validated.snpsPerChr.txt

zcat $HOME/refdata/mm10pseudoSpretus/spretus.brain.snps.vcf.gz | grep -v "^#" | grep PASS | cut -f 1 | sort | uniq -c | awk '{printf("chr%s\t%d\n", $2, $1)}' | sort -k1V > brain.validated.snpsPerChr.txt
printf "chrY\t0\n" >> brain.validated.snpsPerChr.txt

more *.snpsPerChr.txt


#####
# Per chrom SNP densities

for varPerchrFile in \
spretus.dbSNP142 \
spretus.extracted \
patski.validated \
brain.validated \
; do

echo $varPerchrFile
paste <( cat mm10.sizes ) <( cut -f 2 $varPerchrFile.snpsPerChr.txt ) > $varPerchrFile.snpDensities.tmp

printf "chr\tlength\t%sSNPs\tbasesPerSNP\n" > $varPerchrFile.snpDensities.tsv
cat  $varPerchrFile.snpDensities.tmp | awk 'BEGIN{FS="\t"; OFS="\t"}
{
if($3!=0){snpDensity=$2/$3; snpDensity=sprintf("%0.2f", snpDensity)} else {snpDensity="nan"};
printf("%s\t%d\t%d\t%0.2f\n", $1, $2, $3, snpDensity);
totBases+=$2; totSNPs+=$3; totIndels+=$4
}
END{
snpDensity=totBases/totSNPs;
printf("%s\t%d\t%d\t%0.2f\n", "genome", totBases, totSNPs, snpDensity);
}' >> $varPerchrFile.snpDensities.tsv

rm $varPerchrFile.snpDensities.tmp

# HTML
/net/noble/vol1/home/noble/bin/tsv2html.py "$varPerchrFile.snpDensities.tsv" > "$varPerchrFile.snpDensities.html"

done


#####
# Aggregate SNP densities

cat spretus.dbSNP142.snpDensities.tsv > mm10.vcfs.snpDensities.tsv

for varPerchrFile in \
spretus.extracted \
patski.validated \
brain.validated \
; do

echo $varPerchrFile

paste <( cat mm10.vcfs.snpDensities.tsv ) <( cat $varPerchrFile.snpDensities.tsv | awk 'BEGIN{OFS="\t"} {print $3, $4}' ) > mm10.vcfs.snpDensities.tsv.tmp
mv mm10.vcfs.snpDensities.tsv.tmp mm10.vcfs.snpDensities.tsv

done

printf "\t\tspretus.dbSNP142\t\tspretus.extracted\t\tpatski.validated\t\tbrain.validated\t\n" > header.snpDensities.tmp
cat header.snpDensities.tmp mm10.vcfs.snpDensities.tsv > mm10.vcfs.snpDensities.tsv.tmp
mv mm10.vcfs.snpDensities.tsv.tmp mm10.vcfs.snpDensities.tsv
rm header.snpDensities.tmp

# HTML
/net/noble/vol1/home/noble/bin/tsv2html.py "mm10.vcfs.snpDensities.tsv" > "mm10.vcfs.snpDensities.html"

more *.snpDensities.tsv


#####
#####
# Per chrom Indels

zcat $HOME/refdata/mm10pseudoSpretus.dbSNP142/SPRET_EiJ.mgp.v5.indels.dbSNP142.normed.vcf.gz | grep -v "^#" | grep PASS | cut -f 1 | sort | uniq -c | awk '{printf("chr%s\t%d\n", $2, $1)}' | sort -k1V > spretus.dbSNP142.indelsPerChr.txt

zcat $HOME/refdata/mm10pseudoSpretus/spretus.indels.vcf.gz | grep -v "^#" | grep PASS | cut -f 1 | sort | uniq -c | awk '{printf("chr%s\t%d\n", $2, $1)}' | sort -k1V > spretus.extracted.indelsPerChr.txt

more *.indelsPerChr.txt


#####
# Per chrom indel densities

for varPerchrFile in \
spretus.dbSNP142 \
spretus.extracted \
; do

echo $varPerchrFile
paste <( cat mm10.sizes ) <( cut -f 2 $varPerchrFile.snpsPerChr.txt ) > $varPerchrFile.indelDensities.tmp

printf "chr\tlength\t%sIndels\tbasesPerIndel\n" > $varPerchrFile.indelDensities.tsv
cat  $varPerchrFile.indelDensities.tmp | awk 'BEGIN{FS="\t"; OFS="\t"}
{
if($3!=0){snpDensity=$2/$3; snpDensity=sprintf("%0.2f", snpDensity)} else {snpDensity="nan"};
printf("%s\t%d\t%d\t%0.2f\n", $1, $2, $3, snpDensity);
totBases+=$2; totSNPs+=$3; totIndels+=$4
}
END{
snpDensity=totBases/totSNPs;
printf("%s\t%d\t%d\t%0.2f\n", "genome", totBases, totSNPs, snpDensity);
}' >> $varPerchrFile.indelDensities.tsv

rm $varPerchrFile.indelDensities.tmp

# HTML
/net/noble/vol1/home/noble/bin/tsv2html.py "$varPerchrFile.indelDensities.tsv" > "$varPerchrFile.indelDensities.html"

done


#####
# Aggregate indel densities

cat spretus.dbSNP142.indelDensities.tsv > mm10.vcfs.indelDensities.tsv

for varPerchrFile in \
spretus.extracted \
; do

echo $varPerchrFile

paste <( cat mm10.vcfs.indelDensities.tsv ) <( cat $varPerchrFile.indelDensities.tsv | awk 'BEGIN{OFS="\t"} {print $3, $4}' ) > mm10.vcfs.indelDensities.tsv.tmp
mv mm10.vcfs.indelDensities.tsv.tmp mm10.vcfs.indelDensities.tsv

done

printf "\t\tspretus.dbSNP142\t\tspretus.extracted\t\n" > header.indelDensities.tmp
cat header.indelDensities.tmp mm10.vcfs.indelDensities.tsv > mm10.vcfs.indelDensities.tsv.tmp
mv mm10.vcfs.indelDensities.tsv.tmp mm10.vcfs.indelDensities.tsv
rm header.indelDensities.tmp

# HTML
/net/noble/vol1/home/noble/bin/tsv2html.py "mm10.vcfs.indelDensities.tsv" > "mm10.vcfs.indelDensities.html"

more *.indelDensities.tsv




















##################################################################################################
##################################################################################################
##################################################################################################
# 180312
# Mappability data: only available for mm9
#
# https://davetang.org/muse/2013/07/08/encode-mappability/
# http://genome-preview.soe.ucsc.edu/cgi-bin/hgTrackUi?hgsid=2810193_hgP8XmpWWWLB04uqhv8EZhmkjpmn&c=chr12&g=wgEncodeMapability
# http://hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeMapability/

# Mappability file
mkdir $HOME/refdata/mm9mappability
cd $HOME/refdata/mm9mappability
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig

# Convert bigwig to bedgraph
 ~/bin/UCSC.utilities.20160601/bigWigToBedGraph $HOME/refdata/mm9Mapability/wgEncodeCrgMapabilityAlign75mer.bigWig $HOME/refdata/mm9Mapability/wgEncodeCrgMapabilityAlign75mer.bedgraph

# mm9 to mm10 chaing for liftover
mkdir $HOME/refdata/mm9Liftover
cd $HOME/refdata/mm9Liftover
wget hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz .

# liftover mm9 mappabilities to mm10
mkdir $HOME/refdata/mm10mappability
cd $HOME/refdata/mm10mappability
 ~/bin/UCSC.utilities.20160601/liftOver $HOME/refdata/mm9Mapability/wgEncodeCrgMapabilityAlign75mer.bedgraph $HOME/refdata/mm9Liftover/mm9ToMm10.over.chain.gz wgEncodeCrgMapabilityAlign75mer.mm10.bedgraph wgEncodeCrgMapabilityAlign75mer.mm9tomm10unlifted.bedgraph

# sort mappability
mkdir /scratch/sorttemp
cat wgEncodeCrgMapabilityAlign75mer.mm10.bedgraph | sort -T /scratch/sorttemp -k1V -k2 > wgEncodeCrgMapabilityAlign75mer.mm10.sorted.bedgraph

# # clean chromInfo
# grep -v "_" ~/refdata/m10/chromInfo.txt | sort -k1V -k2 > ~/refdata/m10/chromInfo.cleaned.txt

# Average mappability w/in 100mers
~/bin/pythonScripts/binbedgraph.py wgEncodeCrgMapabilityAlign75mer.mm10.sorted.bedgraph ~/refdata/mm10/chromInfo.txt 100 1 True wgEncodeCrgMapabilityAlign75mer.mm10.100bp.bedgraph

# check
grep -v nan wgEncodeCrgMapabilityAlign75mer.mm10.100bp.bedgraph | head
# chr1    3012300 3012400 1.00000000
# chr1    3014700 3014800 0.11866989 <--- *
# chr1    3014800 3014900 0.01029242
# chr1    3014900 3015000 0.02346738
# chr1    3015000 3015100 0.02492471
# chr1    3015100 3015200 0.00524817
# chr1    3015200 3015300 0.20324028
# chr1    3015300 3015400 0.17148243
# chr1    3015400 3015500 0.15775181
# chr1    3015500 3015600 0.07673367

cat wgEncodeCrgMapabilityAlign75mer.mm10.sorted.bedgraph | awk -v start=3014700 '$1=="chr1" && (($2>=(start+1) && $3<=(start+100)) || ($2<=start && $3>=(start+100))){print $0; tally++; value+=$4}; $1=="chr2"{print value, tally, value/tally; exit}'
cat wgEncodeCrgMapabilityAlign75mer.mm10.sorted.bedgraph | awk -v start=3014700 '$1=="chr1" && (($2+(($3-$2)/2))>=(start+1) && ($2+(($3-$2)/2))<=(start+100)){print $0; tally++; value+=$4}; $1=="chr2"{print value, tally, value/tally; exit}'
# chr1    3014743 3014744 0.333333
# chr1    3014744 3014745 0.125
# chr1    3014745 3014747 0.2
# chr1    3014747 3014748 0.1
# chr1    3014748 3014749 0.125
# chr1    3014749 3014750 0.1
# chr1    3014750 3014751 0.05
# chr1    3014751 3014753 0.0769231
# chr1    3014753 3014756 0.125
# chr1    3014756 3014758 0.1
# chr1    3014758 3014760 0.0769231
# chr1    3014760 3014763 0.0625
# chr1    3014763 3014766 0.0769231
# chr1    3014766 3014767 0.125
# chr1    3014767 3014771 0.1
# chr1    3014771 3014772 0.125
# chr1    3014772 3014775 0.166667
# chr1    3014775 3014790 0.1
# chr1    3014790 3014793 0.166667
# chr1    3014793 3014801 0.0384615
# 2.3734 20 0.11867




















##################################################################################################
##################################################################################################
##################################################################################################
# 180416
# hg19 annotations

mkdir ~/refdata/hg19
cd  ~/refdata/hg19


###################
# 2bit
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/hg19.2bit .


###################
# UCSC known genes
# Via UCSC browser
# http://genome.ucsc.edu/cgi-bin/hgTables?command=start
# Use 'selected fields from primary and related tables' in 'output format' dropdown.
# Download to PC

wc -l hg19knownGenes.tsv
# 82961 hg19knownGenes.tsv

# Format as BED
cat hg19knownGenes.tsv | awk 'BEGIN{OFS="\t"}NR>1{print $1, $3, $4, $10, "nan", $2, $5, $6, $7, $8, $9}' > hg19knownGenes.bed

# chr 19 genes
grep "^chr19[[:blank:]]" hg19knownGenes.bed > hg19knownGenes.chr19.bed
# 4711


##################
# UCSC refseq genes
# Via UCSC browser
# http://genome.ucsc.edu/cgi-bin/hgTables?command=start
# Use 'selected fields from primary and related tables' in 'output format' dropdown.
# Download to PC

wc -l hg19refSeq.tsv
# 61178 hg19refSeq.tsv

# Format as BED
cat hg19refSeq.tsv | awk 'BEGIN{OFS="\t"}NR>1{print $1, $3, $4, $11, $10, $2, $5, $6, $7, $8, $9}' > hg19refSeq.bed

# chr 19 genes
grep "^chr19[[:blank:]]" hg19refSeq.bed > hg19refSeq.chr19.bed
# 3681


###################
# DEC up/down genes

# Make sure these have been given the tr '\r' '\n' treatment
DECupGenes="$HOME/proj/2018_CUTnRUN_analysis/results/gbonora/20180410_CUTnRUN_analysis/externalData/DEscRNAseqData/DECupGenes.txt"
DECdnGenes="$HOME/proj/2018_CUTnRUN_analysis/results/gbonora/20180410_CUTnRUN_analysis/externalData/DEscRNAseqData/DECdnGenes.txt"


cat /dev/null > hg19refSeq.DECupGenes.bed
while read pat; do
echo $pat
grep "[[:blank:]]${pat}[[:blank:]]" hg19refSeq.bed >> hg19refSeq.DECupGenes.bed
done < ${DECupGenes}

cat /dev/null > hg19refSeq.DECdnGenes.bed
while read pat; do
echo $pat
grep "[[:blank:]]${pat}[[:blank:]]" hg19refSeq.bed >> hg19refSeq.DECdnGenes.bed
#>> hg19refSeq.DECdnGenes.bed
done < ${DECdnGenes}


###################
# Enhancer Atlas H1 enhancers
# http://enhanceratlas.org/data/enhseq/H1.fasta
wget http://enhanceratlas.org/data/enhseq/H1.fasta
mv  H1.fasta  H1.enhancerSeqs.fasta

# >chr1:629800-630330_
cat  H1.enhancerSeqs.fasta | grep "^>" | sed 's/>//' | sed 's/_.*$//' | tr ':' '\t' | tr '-' '\t' | gzip > H1.enhancerSeqs.bed.gz
zcat H1.enhancerSeqs.bed.gz | wc -l
# 58821

zcat H1.enhancerSeqs.bed.gz | grep "^chr19" > H1.enhancerSeqs.chr19.bed
# 3495


###################
# Epigenome roadmap chromHMM H1 (E003) enhancers
# http://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
#
# 18 state model:
#
# STATE NO.	MNEMONIC	DESCRIPTION	COLOR NAME	COLOR CODE
# 1	TssA	Active TSS	Red	255,0,0
# 2	TssFlnk	Flanking TSS	Orange Red	255,69,0
# 3	TssFlnkU	Flanking TSS Upstream	Orange Red	255,69,0
# 4	TssFlnkD	Flanking TSS Downstream	Orange Red	255,69,0
# 5	Tx	Strong transcription	Green	0,128,0
# 6	TxWk	Weak transcription	DarkGreen	0,100,0
# 7	EnhG1	Genic enhancer1	GreenYellow	194,225,5
# 8	EnhG2	Genic enhancer2	GreenYellow	194,225,5
# 9	EnhA1	Active Enhancer 1	Orange	255,195,77
# 10	EnhA2	Active Enhancer 2	Orange	255,195,77
# 11	EnhWk	Weak Enhancer	Yellow	255,255,0
# 12	ZNF/Rpts	ZNF genes & repeats	Medium Aquamarine	102,205,170
# 13	Het	Heterochromatin	PaleTurquoise	138,145,208
# 14	TssBiv	Bivalent/Poised TSS	IndianRed	205,92,92
# 15	EnhBiv	Bivalent Enhancer	DarkKhaki	189,183,107
# 16	ReprPC	Repressed PolyComb	Silver	128,128,128
# 17	ReprPCWk	Weak Repressed PolyComb	Gainsboro	192,192,192
# 18	Quies	Quiescent/Low	White	255,255,255

# Downloaded to PC and uploaded to server
SRCDIR="$HOME/proj/2018_CUTnRUN_analysis/results/gbonora/20180410_CUTnRUN_analysis/externalData/18state.chromHMM.all.dense.browserFiles/E003_18_core_K27ac_dense.bed.gz"

# chr19
zcat "${SRCDIR}" | grep -e "^chr19[[:blank:]]" > H1_18chromHMM.chr19.bed
cat H1_18chromHMM.chr19.bed | grep -e "_TssA" > H1_TssA.chr19.bed
cat H1_18chromHMM.chr19.bed | grep -e "_TssBiv" > H1_TssBiv.chr19.bed
cat H1_18chromHMM.chr19.bed | grep -e "_Tx[[:blank:]]" > H1_Tx.chr19.bed
cat H1_18chromHMM.chr19.bed | grep -e "_TxWk" > H1_TxWk.chr19.bed
cat H1_18chromHMM.chr19.bed | grep -e "_EnhA" > H1_EnhA.chr19.bed
cat H1_18chromHMM.chr19.bed | grep -e "_EnhG" > H1_EnhG.chr19.bed
cat H1_18chromHMM.chr19.bed | grep -e "_EnhWk" > H1_EnhWk.chr19.bed
cat H1_18chromHMM.chr19.bed | grep -e "_EnhBiv" > H1_EnhBiv.chr19.bed
cat H1_18chromHMM.chr19.bed | grep -e "_ReprPC" > H1_ReprPC.chr19.bed
cat H1_18chromHMM.chr19.bed | grep -e "_Het" > H1_Het.chr19.bed
cat H1_18chromHMM.chr19.bed | grep -e "_ZNF" > H1_ZNF.chr19.bed

wc -l H1_*.chr19.bed
#   22479 H1_18chromHMM.chr19.bed
#    1169 H1_EnhA.chr19.bed
#     623 H1_EnhBiv.chr19.bed
#     503 H1_EnhG.chr19.bed
#    2688 H1_EnhWk.chr19.bed
#    2200 H1_Het.chr19.bed
#    1681 H1_ReprPC.chr19.bed
#    1235 H1_TssA.chr19.bed
#     599 H1_TssBiv.chr19.bed
#    1546 H1_Tx.chr19.bed
#    3420 H1_TxWk.chr19.bed
#    2140 H1_ZNF.chr19.bed


# Autosomes
zcat "${SRCDIR}" | grep -e "^chr[0-9]*[[:blank:]]" > H1_18chromHMM.chrA.bed
cat H1_18chromHMM.chrA.bed | grep -e "_TssA" > H1_TssA.chrA.bed
cat H1_18chromHMM.chrA.bed | grep -e "_TssBiv" > H1_TssBiv.chrA.bed
cat H1_18chromHMM.chrA.bed | grep -e "_Tx[[:blank:]]" > H1_Tx.chrA.bed
cat H1_18chromHMM.chrA.bed | grep -e "_TxWk" > H1_TxWk.chrA.bed
cat H1_18chromHMM.chrA.bed | grep -e "_EnhA" > H1_EnhA.chrA.bed
cat H1_18chromHMM.chrA.bed | grep -e "_EnhG" > H1_EnhG.chrA.bed
cat H1_18chromHMM.chrA.bed | grep -e "_EnhWk" > H1_EnhWk.chrA.bed
cat H1_18chromHMM.chrA.bed | grep -e "_EnhBiv" > H1_EnhBiv.chrA.bed
cat H1_18chromHMM.chrA.bed | grep -e "_ReprPC" > H1_ReprPC.chrA.bed
cat H1_18chromHMM.chrA.bed | grep -e "_Het" > H1_Het.chrA.bed
cat H1_18chromHMM.chrA.bed | grep -e "_ZNF" > H1_ZNF.chrA.bed

wc -l H1_*.chrA.bed
#   548223 H1_18chromHMM.chrA.bed
#    48589 H1_EnhA.chrA.bed
#    14593 H1_EnhBiv.chrA.bed
#     6725 H1_EnhG.chrA.bed
#   103587 H1_EnhWk.chrA.bed
#    47385 H1_Het.chrA.bed
#    40162 H1_ReprPC.chrA.bed
#    14958 H1_TssA.chrA.bed
#    10837 H1_TssBiv.chrA.bed
#    25725 H1_Tx.chrA.bed
#    90849 H1_TxWk.chrA.bed
#    28901 H1_ZNF.chrA.bed


###################
# Then rsync to server
rsync -avuPLn  --exclude ".DS_Store"  ~/refdata/hg19 gs:refdata




















##################################################################################################
##################################################################################################
##################################################################################################
# 180515
# mm10 gene sets

mkdir ~/refdata/mm10
cd  ~/refdata/mm10


##################
# UCSC known genes genes
# Via UCSC browser
# http://genome.ucsc.edu/cgi-bin/hgTables?command=start
# Use 'selected fields from primary and related tables' in 'output format' dropdown.
# Download to PC

wc -l mm10.UCSC.knownGenes.tsv
# 63760

# Format as BED
cat mm10.UCSC.knownGenes.tsv | awk 'BEGIN{OFS="\t"}NR>1{print $2, $4, $5, $11, $1, $3, $6, $7, $8, $9, $10}' > mm10.UCSC.knownGenes.bed

wc -l mm10.UCSC.knownGenes.bed
# 63759 mm10.UCSC.knownGenes.bed

# chrX genes
grep "^chrX[[:blank:]]" mm10.UCSC.knownGenes.bed > mm10.UCSC.knownGenes.chrX.bed
# 2430


##################
# 20190430
# Diplosomal genes
# (For DelFirreXa, exclude chromosomes 2,3,4,5,7,8,9,12, and 14.)
grep "^chr1[[:blank:]]" mm10.UCSC.knownGenes.bed > mm10.UCSC.knownGenes.diplosomes.bed
grep "^chr6[[:blank:]]" mm10.UCSC.knownGenes.bed >> mm10.UCSC.knownGenes.diplosomes.bed
grep "^chr10[[:blank:]]" mm10.UCSC.knownGenes.bed >> mm10.UCSC.knownGenes.diplosomes.bed
grep "^chr11[[:blank:]]" mm10.UCSC.knownGenes.bed >> mm10.UCSC.knownGenes.diplosomes.bed
grep "^chr13[[:blank:]]" mm10.UCSC.knownGenes.bed >> mm10.UCSC.knownGenes.diplosomes.bed
grep "^chr15[[:blank:]]" mm10.UCSC.knownGenes.bed >> mm10.UCSC.knownGenes.diplosomes.bed
grep "^chr16[[:blank:]]" mm10.UCSC.knownGenes.bed >> mm10.UCSC.knownGenes.diplosomes.bed
grep "^chr17[[:blank:]]" mm10.UCSC.knownGenes.bed >> mm10.UCSC.knownGenes.diplosomes.bed
grep "^chr18[[:blank:]]" mm10.UCSC.knownGenes.bed >> mm10.UCSC.knownGenes.diplosomes.bed
grep "^chr19[[:blank:]]" mm10.UCSC.knownGenes.bed >> mm10.UCSC.knownGenes.diplosomes.bed
wc -l mm10.UCSC.knownGenes.diplosomes.bed
# 28003 mm10.UCSC.knownGenes.diplosomes.bed


##################
# Escape genes
ESCAPEES=(
Xist
Zfp300
Pbdc1
Ftx
Eif2s3x
Kdm6a
Kdm5c
5530601H04Rik
Xkrx
Atp6ap1
Rp2h
Ddx3x
1810030O07Rik
Mid1
Shroom4
Gm14827
Firre
Slc16a2
Car5b
Magt1
Tmem29
Efnb1
Pim2
Hmgn5
Gnl3l
Yipf6
Asmt
Sh3bgrl
Rlim)

echo "${#ESCAPEES[@]}"
# 29

cat /dev/null > mm10.UCSC.knownGenes.chrX.escapees.bed
for theEscapee in "${ESCAPEES[@]}"; do
echo "${theEscapee}"
awk -v escapee="${theEscapee}" 'BEGIN{OFS="\t"} $4==escapee{print $0}' mm10.UCSC.knownGenes.chrX.bed >> mm10.UCSC.knownGenes.chrX.escapees.bed
done

wc -l mm10.UCSC.knownGenes.chrX.escapees.bed
# 88


########################
# Conserved escape genes
# 2018108

CONSESCAPEES=(
1810030O07Rik
Ddx3x
Kdm6a
Pbdc1
Eif2s3x
Xist
Mid1
5530601H04Rik
Kdm5c
)

echo "${#CONSESCAPEES[@]}"
# 9

cat /dev/null > mm10.UCSC.knownGenes.chrX.conservedEscapees.bed
for theEscapee in "${CONSESCAPEES[@]}"; do
echo "${theEscapee}"
awk -v escapee="${theEscapee}" 'BEGIN{OFS="\t"} $4==escapee{print $0}' mm10.UCSC.knownGenes.chrX.bed >> mm10.UCSC.knownGenes.chrX.conservedEscapees.bed
done

wc -l mm10.UCSC.knownGenes.chrX.conservedEscapees.bed
# 42




##################
# UCSC known genes promoters - 2kb
# Via UCSC browser
# http://genome.ucsc.edu/cgi-bin/hgTables?command=start
# Use 'selected fields from primary and related tables' in 'output format' dropdown.
# Download to PC

wc -l mm10.UCSC.knownGenePromoters2kb.tsv
# 63760

# Format as BED
cat mm10.UCSC.knownGenePromoters2kb.tsv | awk 'BEGIN{OFS="\t"}NR>1{print $0}' > mm10.UCSC.knownGenePromoters2kb.bed

# chr1 genes
grep "^chr1[[:blank:]]" mm10.UCSC.knownGenePromoters2kb.bed > mm10.UCSC.knownGenePromoters2kb.chr1.bed
# 2430

# chrX genes
grep "^chrX[[:blank:]]" mm10.UCSC.knownGenePromoters2kb.bed > mm10.UCSC.knownGenePromoters2kb.chrX.bed
# 2430



##################
# UCSC CGIs
# Via UCSC browser
# http://genome.ucsc.edu/cgi-bin/hgTables?command=start
# Use 'selected fields from primary and related tables' in 'output format' dropdown.
# Download to PC

wc -l mm10.UCSC.CGIs.tsv
# 16024 - 1

cat mm10.UCSC.CGIs.tsv

# Format as BED
cat mm10.UCSC.CGIs.tsv | awk 'BEGIN{OFS="\t"}NR>1{print $1, $2, $3}' > mm10.UCSC.CGIs.bed

# chrX genes
grep "^chrX[[:blank:]]" mm10.UCSC.CGIs.bed > mm10.UCSC.CGIs.chrX.bed
# 458


###################
# Then rsync to server
rsync -avuPLn  --exclude ".DS_Store"  ~/refdata/mm10/ gs:refdata/mm10/




















##################################################################################################
##################################################################################################
##################################################################################################
# Firre project genes


#################################################
#################################################
#################################################
# 20190426
# DelFirre DE genes

cd  ~/refdata/mm10

hexdump -c "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreUpGenes.txt"  | head
# 0000000       2   8   1   0   4   0   8   I   1   1   R   i   k  \r  \n
# 0000010   4   9   3   0   5   5   8   J   1   8   R   i   k  \r  \n   6
# 0000020   4   3   0   7   0   6   D   2   2   R   i   k  \r  \n   A   7
# 0000030   3   0   0   0   8   H   2   3   R   i   k  \r  \n   A   s   p
# 0000040   m  \r  \n   B   C   0   5   5   3   2   4  \r  \n   B   a   r
# 0000050   d   1  \r  \n   C   d   k   1   5  \r  \n   C   e   n   p   f
# 0000060  \r  \n   C   e   n   p   l  \r  \n   D   t   l  \r  \n   E   x
# 0000070   o   1  \r  \n   F   m   o   1  \r  \n   G   m   1   9   7   0
# 0000080   5  \r  \n   H   j   u   r   p  \r  \n   I   f   i   2   0   3
# 0000090  \r  \n   I   f   i   2   0   5  \r  \n   I   l   d   r   2  \r

# Clean data
cat "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreUpGenes.txt"  | tr -d '\r' > "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreUpGenes.txt.tmp"
mv "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreUpGenes.txt.tmp" "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreUpGenes.txt"

cat "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreDownGenes.txt" | tr -d '\r' > "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreDownGenes.txt.tmp"
mv "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreDownGenes.txt.tmp" "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreDownGenes.txt"

DFUPGENES=( $( cat "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreUpGenes.txt" ) )
echo "${#DFUPGENES[@]}"
# 384

DFDOWNGENES=( $( cat "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreDownGenes.txt" ) )
echo "${#DFDOWNGENES[@]}"
# 276

cat /dev/null > mm10.UCSC.knownGenes.DelFirreDEup.bed
cat /dev/null > mm10.UCSC.knownGenes.chr1.DelFirreDEup.bed
cat /dev/null > mm10.UCSC.knownGenes.chrX.DelFirreDEup.bed
cat /dev/null > mm10.UCSC.knownGenes.diplosomes.DelFirreDEup.bed
for theGene in "${DFUPGENES[@]}"; do
echo ">${theGene}<"
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.bed >> mm10.UCSC.knownGenes.DelFirreDEup.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.chr1.bed >> mm10.UCSC.knownGenes.chr1.DelFirreDEup.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.chrX.bed >> mm10.UCSC.knownGenes.chrX.DelFirreDEup.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.diplosomes.bed >> mm10.UCSC.knownGenes.diplosomes.DelFirreDEup.bed
done

wc -l mm10.UCSC.knownGenes.*DelFirreDEup.bed
#    944 mm10.UCSC.knownGenes.DelFirreDEup.bed
#    119 mm10.UCSC.knownGenes.chr1.DelFirreDEup.bed
#     45 mm10.UCSC.knownGenes.chrX.DelFirreDEup.bed
#    899 mm10.UCSC.knownGenes.diplosomes.DelFirreDEup.bed

cat /dev/null > mm10.UCSC.knownGenes.DelFirreDEdown.bed
cat /dev/null > mm10.UCSC.knownGenes.chr1.DelFirreDEdown.bed
cat /dev/null > mm10.UCSC.knownGenes.chrX.DelFirreDEdown.bed
cat /dev/null > mm10.UCSC.knownGenes.diplosomes.DelFirreDEdown.bed
for theGene in "${DFDOWNGENES[@]}"; do
echo ">${theGene}<"
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.bed >> mm10.UCSC.knownGenes.DelFirreDEdown.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.chr1.bed >> mm10.UCSC.knownGenes.chr1.DelFirreDEdown.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.chrX.bed >> mm10.UCSC.knownGenes.chrX.DelFirreDEdown.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.diplosomes.bed >> mm10.UCSC.knownGenes.diplosomes.DelFirreDEdown.bed
done

wc -l mm10.UCSC.knownGenes.*DelFirreDEdown.bed
#    721 mm10.UCSC.knownGenes.DelFirreDEdown.bed
#     80 mm10.UCSC.knownGenes.chr1.DelFirreDEdown.bed
#     41 mm10.UCSC.knownGenes.chrX.DelFirreDEdown.bed
#    673 mm10.UCSC.knownGenes.diplosomes.DelFirreDEdown.bed


cat mm10.UCSC.knownGenes.DelFirreDEup.bed mm10.UCSC.knownGenes.DelFirreDEdown.bed > mm10.UCSC.knownGenes.DelFirreDE.bed
cat mm10.UCSC.knownGenes.chr1.DelFirreDEup.bed mm10.UCSC.knownGenes.chr1.DelFirreDEdown.bed > mm10.UCSC.knownGenes.chr1.DelFirreDE.bed
cat mm10.UCSC.knownGenes.chrX.DelFirreDEup.bed mm10.UCSC.knownGenes.chrX.DelFirreDEdown.bed > mm10.UCSC.knownGenes.chrX.DelFirreDE.bed
cat mm10.UCSC.knownGenes.diplosomes.DelFirreDEup.bed mm10.UCSC.knownGenes.diplosomes.DelFirreDEdown.bed > mm10.UCSC.knownGenes.diplosomes.DelFirreDE.bed
wc -l mm10.UCSC.knownGenes.*DelFirreDE.bed
#   1665 mm10.UCSC.knownGenes.DelFirreDE.bed
#    199 mm10.UCSC.knownGenes.chr1.DelFirreDE.bed
#     86 mm10.UCSC.knownGenes.chrX.DelFirreDE.bed
#   1572 mm10.UCSC.knownGenes.diplosomes.DelFirreDE.bed



#################################################
# 20190430
# DelFirre non DE genes

# Clean data
cat "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenes.txt"  | tr -d '\r' > "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenes.txt.tmp"
mv "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenes.txt.tmp" "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenes.txt"

DFNONDEGENES=( $( cat "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenes.txt" ) )
echo "${#DFNONDEGENES[@]}"
# 853

cat /dev/null > mm10.UCSC.knownGenes.DelFirreNotDE.bed
cat /dev/null > mm10.UCSC.knownGenes.chr1.DelFirreNotDE.bed
cat /dev/null > mm10.UCSC.knownGenes.chrX.DelFirreNotDE.bed
cat /dev/null > mm10.UCSC.knownGenes.diplosomes.DelFirreNotDE.bed
for theGene in "${DFNONDEGENES[@]}"; do
echo ">${theGene}<"
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.bed >> mm10.UCSC.knownGenes.DelFirreNotDE.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.chr1.bed >> mm10.UCSC.knownGenes.chr1.DelFirreNotDE.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.chrX.bed >> mm10.UCSC.knownGenes.chrX.DelFirreNotDE.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.diplosomes.bed >> mm10.UCSC.knownGenes.diplosomes.DelFirreNotDE.bed
done

wc -l mm10.UCSC.knownGenes.*DelFirreNotDE.bed
#    2327 mm10.UCSC.knownGenes.DelFirreNotDE.bed
#     367 mm10.UCSC.knownGenes.chr1.DelFirreNotDE.bed
#       0 mm10.UCSC.knownGenes.chrX.DelFirreNotDE.bed
#    2327 mm10.UCSC.knownGenes.diplosomes.DelFirreNotDE.bed

rm mm10.UCSC.knownGenes.chrX.DelFirreNotDE.bed



#################################################
# 20190503
# DelFirre non DE genes: high and low expression

###########
# Highly expressed non DE genes

# Clean data
cat "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenesHi.txt"  | tr -d '\r' > "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenesHi.txt.tmp"
mv "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenesHi.txt.tmp" "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenesHi.txt"

DFNONDEGENESHI=( $( cat "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenesHi.txt" ) )
echo "${#DFNONDEGENESHI[@]}"
# 200

cat /dev/null > mm10.UCSC.knownGenes.DelFirreNotDEhi.bed
cat /dev/null > mm10.UCSC.knownGenes.chr1.DelFirreNotDEhi.bed
cat /dev/null > mm10.UCSC.knownGenes.chrX.DelFirreNotDEhi.bed
cat /dev/null > mm10.UCSC.knownGenes.diplosomes.DelFirreNotDEhi.bed
for theGene in "${DFNONDEGENESHI[@]}"; do
echo ">${theGene}<"
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.bed >> mm10.UCSC.knownGenes.DelFirreNotDEhi.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.chr1.bed >> mm10.UCSC.knownGenes.chr1.DelFirreNotDEhi.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.chrX.bed >> mm10.UCSC.knownGenes.chrX.DelFirreNotDEhi.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.diplosomes.bed >> mm10.UCSC.knownGenes.diplosomes.DelFirreNotDEhi.bed
done

wc -l mm10.UCSC.knownGenes.*DelFirreNotDEhi.bed
#    554 mm10.UCSC.knownGenes.DelFirreNotDEhi.bed
#     79 mm10.UCSC.knownGenes.chr1.DelFirreNotDEhi.bed
#      0 mm10.UCSC.knownGenes.chrX.DelFirreNotDEhi.bed
#    554 mm10.UCSC.knownGenes.diplosomes.DelFirreNotDEhi.bed

rm mm10.UCSC.knownGenes.chrX.DelFirreNotDEhi.bed


###########
# Lowly expressed non DE genes

# Clean data
cat "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenesLo.txt"  | tr -d '\r' > "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenesLo.txt.tmp"
mv "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenesLo.txt.tmp" "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenesLo.txt"
   
DFNONDEGENESLO=( $( cat "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenesLo.txt" ) )
echo "${#DFNONDEGENESLO[@]}"
# 200

cat /dev/null > mm10.UCSC.knownGenes.DelFirreNotDElo.bed
cat /dev/null > mm10.UCSC.knownGenes.chr1.DelFirreNotDElo.bed
cat /dev/null > mm10.UCSC.knownGenes.chrX.DelFirreNotDElo.bed
cat /dev/null > mm10.UCSC.knownGenes.diplosomes.DelFirreNotDElo.bed
for theGene in "${DFNONDEGENESLO[@]}"; do
echo ">${theGene}<"
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.bed >> mm10.UCSC.knownGenes.DelFirreNotDElo.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.chr1.bed >> mm10.UCSC.knownGenes.chr1.DelFirreNotDElo.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.chrX.bed >> mm10.UCSC.knownGenes.chrX.DelFirreNotDElo.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.diplosomes.bed >> mm10.UCSC.knownGenes.diplosomes.DelFirreNotDElo.bed
done

wc -l mm10.UCSC.knownGenes.*DelFirreNotDElo.bed
#    537 mm10.UCSC.knownGenes.DelFirreNotDElo.bed
#    112 mm10.UCSC.knownGenes.chr1.DelFirreNotDElo.bed
#      0 mm10.UCSC.knownGenes.chrX.DelFirreNotDElo.bed
#    537 mm10.UCSC.knownGenes.diplosomes.DelFirreNotDElo.bed

rm mm10.UCSC.knownGenes.chrX.DelFirreNotDElo.bed










#################################################
#################################################
#################################################
# *** 20190503 ***
# And filter out duplicate known genes
# Take canonical gene as that one which is the longest

cd  ~/refdata/mm10

# Some known genes are highly over-represented!
cat mm10.UCSC.knownGenes.bed | cut -f 4 | sort | uniq -c | sort -k1,1n  | tail -n 10
#      33 Gm20736
#      33 Pcdh15
#      37 AB335791
#      38 TRNA_Ala
#      46 AB349653
#      57 TRNA_Cys
#      60 AB353040
#      63 AB349811
#     105 AB351889
#     189 JA187517


cat mm10.UCSC.knownGenes.tsv | sort -k11,11 | awk 'BEGIN{OFS="\t"; lastSymbol="X-X-X"}
NR>1{ \
geneLength=$5-$4;
geneInfo=sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", $2, $4, $5, $11, $1, $3, $6, $7, $8, $9, $10);
# print lastSymbol, geneInfo;
if ($11==lastSymbol) {
if (geneLength>longGeneLength) { longGeneLength=geneLength; longGeneInfo=geneInfo }
}
else { 
if (lastSymbol != "X-X-X") { print longGeneInfo };
lastSymbol=$11; lastGeneLength=geneLength; longGeneInfo=geneInfo;
}
}
END {
print longGeneInfo 
}' > mm10.UCSC.knownGenes.unique.bed

wc -l mm10.UCSC.knownGenes.unique.bed
# 32695 mm10.UCSC.knownGenes.unique.bed

cat mm10.UCSC.knownGenes.unique.bed | cut -f 4 | sort | uniq -c | sort -k1,1n  | tail -n 10
#       1 mKIAA0208
#       1 mKIAA0217
#       1 mKIAA0383
#       1 mKIAA1554
#       1 mKIAA1630
#       1 mannose
#       1 mdab1
#       1 mm10.kgXref.geneSymbol
#       1 rjs
#       1 unknown

# chr1 genes
grep "^chr1[[:blank:]]" mm10.UCSC.knownGenes.unique.bed > mm10.UCSC.knownGenes.unique.chr1.bed
# 1131 mm10.UCSC.knownGenes.unique.chrX.bed

# chrX genes
grep "^chrX[[:blank:]]" mm10.UCSC.knownGenes.unique.bed > mm10.UCSC.knownGenes.unique.chrX.bed
# 1682 mm10.UCSC.knownGenes.unique.chr1.bed


##################
# Diplosomal genes
# (For DelFirreXa, exclude chromosomes 2,3,4,5,7,8,9,12, and 14.)
grep "^chr1[[:blank:]]" mm10.UCSC.knownGenes.unique.bed > mm10.UCSC.knownGenes.unique.diplosomes.bed
grep "^chr6[[:blank:]]" mm10.UCSC.knownGenes.unique.bed >> mm10.UCSC.knownGenes.unique.diplosomes.bed
grep "^chr10[[:blank:]]" mm10.UCSC.knownGenes.unique.bed >> mm10.UCSC.knownGenes.unique.diplosomes.bed
grep "^chr11[[:blank:]]" mm10.UCSC.knownGenes.unique.bed >> mm10.UCSC.knownGenes.unique.diplosomes.bed
grep "^chr13[[:blank:]]" mm10.UCSC.knownGenes.unique.bed >> mm10.UCSC.knownGenes.unique.diplosomes.bed
grep "^chr15[[:blank:]]" mm10.UCSC.knownGenes.unique.bed >> mm10.UCSC.knownGenes.unique.diplosomes.bed
grep "^chr16[[:blank:]]" mm10.UCSC.knownGenes.unique.bed >> mm10.UCSC.knownGenes.unique.diplosomes.bed
grep "^chr17[[:blank:]]" mm10.UCSC.knownGenes.unique.bed >> mm10.UCSC.knownGenes.unique.diplosomes.bed
grep "^chr18[[:blank:]]" mm10.UCSC.knownGenes.unique.bed >> mm10.UCSC.knownGenes.unique.diplosomes.bed
grep "^chr19[[:blank:]]" mm10.UCSC.knownGenes.unique.bed >> mm10.UCSC.knownGenes.unique.diplosomes.bed
wc -l mm10.UCSC.knownGenes.unique.diplosomes.bed
# 14433 mm10.UCSC.knownGenes.unique.diplosomes.bed


##################
# Escape genes
ESCAPEES=(
Xist
Zfp300
Pbdc1
Ftx
Eif2s3x
Kdm6a
Kdm5c
5530601H04Rik
Xkrx
Atp6ap1
Rp2h
Ddx3x
1810030O07Rik
Mid1
Shroom4
Gm14827
Firre
Slc16a2
Car5b
Magt1
Tmem29
Efnb1
Pim2
Hmgn5
Gnl3l
Yipf6
Asmt
Sh3bgrl
Rlim)

echo "${#ESCAPEES[@]}"
# 29

cat /dev/null > mm10.UCSC.knownGenes.unique.chrX.escapees.bed
for theEscapee in "${ESCAPEES[@]}"; do
echo "${theEscapee}"
awk -v escapee="${theEscapee}" 'BEGIN{OFS="\t"} $4==escapee{print $0}' mm10.UCSC.knownGenes.unique.chrX.bed >> mm10.UCSC.knownGenes.unique.chrX.escapees.bed
done

wc -l mm10.UCSC.knownGenes.unique.chrX.escapees.bed
# 28


########################
# Conserved escape genes

CONSESCAPEES=(
1810030O07Rik
Ddx3x
Kdm6a
Pbdc1
Eif2s3x
Xist
Mid1
5530601H04Rik
Kdm5c
)

echo "${#CONSESCAPEES[@]}"
# 9

cat /dev/null > mm10.UCSC.knownGenes.unique.chrX.conservedEscapees.bed
for theEscapee in "${CONSESCAPEES[@]}"; do
echo "${theEscapee}"
awk -v escapee="${theEscapee}" 'BEGIN{OFS="\t"} $4==escapee{print $0}' mm10.UCSC.knownGenes.unique.chrX.bed >> mm10.UCSC.knownGenes.unique.chrX.conservedEscapees.bed
done

wc -l mm10.UCSC.knownGenes.unique.chrX.conservedEscapees.bed
# 9





#################################################
# DelFirre DE genes

DFUPGENES=( $( cat "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreUpGenes.txt" ) )
echo "${#DFUPGENES[@]}"
# 384

DFDOWNGENES=( $( cat "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreDownGenes.txt" ) )
echo "${#DFDOWNGENES[@]}"
# 276

cat /dev/null > mm10.UCSC.knownGenes.unique.DelFirreDEup.bed
cat /dev/null > mm10.UCSC.knownGenes.unique.chr1.DelFirreDEup.bed
cat /dev/null > mm10.UCSC.knownGenes.unique.chrX.DelFirreDEup.bed
cat /dev/null > mm10.UCSC.knownGenes.unique.diplosomes.DelFirreDEup.bed
for theGene in "${DFUPGENES[@]}"; do
echo ">${theGene}<"
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.bed >> mm10.UCSC.knownGenes.unique.DelFirreDEup.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.chr1.bed >> mm10.UCSC.knownGenes.unique.chr1.DelFirreDEup.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.chrX.bed >> mm10.UCSC.knownGenes.unique.chrX.DelFirreDEup.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.diplosomes.bed >> mm10.UCSC.knownGenes.unique.diplosomes.DelFirreDEup.bed
done

wc -l mm10.UCSC.knownGenes.unique.*DelFirreDEup.bed
#    377 mm10.UCSC.knownGenes.unique.DelFirreDEup.bed
#     44 mm10.UCSC.knownGenes.unique.chr1.DelFirreDEup.bed
#     17 mm10.UCSC.knownGenes.unique.chrX.DelFirreDEup.bed
#    360 mm10.UCSC.knownGenes.unique.diplosomes.DelFirreDEup.bed



cat /dev/null > mm10.UCSC.knownGenes.unique.DelFirreDEdown.bed
cat /dev/null > mm10.UCSC.knownGenes.unique.chr1.DelFirreDEdown.bed
cat /dev/null > mm10.UCSC.knownGenes.unique.chrX.DelFirreDEdown.bed
cat /dev/null > mm10.UCSC.knownGenes.unique.diplosomes.DelFirreDEdown.bed
for theGene in "${DFDOWNGENES[@]}"; do
echo ">${theGene}<"
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.bed >> mm10.UCSC.knownGenes.unique.DelFirreDEdown.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.chr1.bed >> mm10.UCSC.knownGenes.unique.chr1.DelFirreDEdown.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.chrX.bed >> mm10.UCSC.knownGenes.unique.chrX.DelFirreDEdown.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.diplosomes.bed >> mm10.UCSC.knownGenes.unique.diplosomes.DelFirreDEdown.bed
done

wc -l mm10.UCSC.knownGenes.unique.*DelFirreDEdown.bed
#    273 mm10.UCSC.knownGenes.unique.DelFirreDEdown.bed
#     31 mm10.UCSC.knownGenes.unique.chr1.DelFirreDEdown.bed
#     15 mm10.UCSC.knownGenes.unique.chrX.DelFirreDEdown.bed
#    255 mm10.UCSC.knownGenes.unique.diplosomes.DelFirreDEdown.bed



cat mm10.UCSC.knownGenes.unique.DelFirreDEup.bed mm10.UCSC.knownGenes.unique.DelFirreDEdown.bed > mm10.UCSC.knownGenes.unique.DelFirreDE.bed
cat mm10.UCSC.knownGenes.unique.chr1.DelFirreDEup.bed mm10.UCSC.knownGenes.unique.chr1.DelFirreDEdown.bed > mm10.UCSC.knownGenes.unique.chr1.DelFirreDE.bed
cat mm10.UCSC.knownGenes.unique.chrX.DelFirreDEup.bed mm10.UCSC.knownGenes.unique.chrX.DelFirreDEdown.bed > mm10.UCSC.knownGenes.unique.chrX.DelFirreDE.bed
cat mm10.UCSC.knownGenes.unique.diplosomes.DelFirreDEup.bed mm10.UCSC.knownGenes.unique.diplosomes.DelFirreDEdown.bed > mm10.UCSC.knownGenes.unique.diplosomes.DelFirreDE.bed
wc -l mm10.UCSC.knownGenes.unique.*DelFirreDE.bed
#    650 mm10.UCSC.knownGenes.unique.DelFirreDE.bed
#     75 mm10.UCSC.knownGenes.unique.chr1.DelFirreDE.bed
#     32 mm10.UCSC.knownGenes.unique.chrX.DelFirreDE.bed
#    615 mm10.UCSC.knownGenes.unique.diplosomes.DelFirreDE.bed



#################################################
# 20190430
# DelFirre non DE genes

DFNONDEGENES=( $( cat "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenes.txt" ) )
echo "${#DFNONDEGENES[@]}"
# 853

cat /dev/null > mm10.UCSC.knownGenes.unique.DelFirreNotDE.bed
cat /dev/null > mm10.UCSC.knownGenes.unique.chr1.DelFirreNotDE.bed
cat /dev/null > mm10.UCSC.knownGenes.unique.chrX.DelFirreNotDE.bed
cat /dev/null > mm10.UCSC.knownGenes.unique.diplosomes.DelFirreNotDE.bed
for theGene in "${DFNONDEGENES[@]}"; do
echo ">${theGene}<"
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.bed >> mm10.UCSC.knownGenes.unique.DelFirreNotDE.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.chr1.bed >> mm10.UCSC.knownGenes.unique.chr1.DelFirreNotDE.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.chrX.bed >> mm10.UCSC.knownGenes.unique.chrX.DelFirreNotDE.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.diplosomes.bed >> mm10.UCSC.knownGenes.unique.diplosomes.DelFirreNotDE.bed
done

wc -l mm10.UCSC.knownGenes.unique.*DelFirreNotDE.bed
#    839 mm10.UCSC.knownGenes.unique.DelFirreNotDE.bed
#    122 mm10.UCSC.knownGenes.unique.chr1.DelFirreNotDE.bed
#      0 mm10.UCSC.knownGenes.unique.chrX.DelFirreNotDE.bed
#    839 mm10.UCSC.knownGenes.unique.diplosomes.DelFirreNotDE.bed

rm mm10.UCSC.knownGenes.unique.chrX.DelFirreNotDE.bed



#################################################
# 20190503
# DelFirre non DE genes: high and low expression

###########
# Highly expressed non DE genes

DFNONDEGENESHI=( $( cat "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenesHi.txt" ) )
echo "${#DFNONDEGENESHI[@]}"
# 200

cat /dev/null > mm10.UCSC.knownGenes.unique.DelFirreNotDEhi.bed
cat /dev/null > mm10.UCSC.knownGenes.unique.chr1.DelFirreNotDEhi.bed
cat /dev/null > mm10.UCSC.knownGenes.unique.chrX.DelFirreNotDEhi.bed
cat /dev/null > mm10.UCSC.knownGenes.unique.diplosomes.DelFirreNotDEhi.bed
for theGene in "${DFNONDEGENESHI[@]}"; do
echo ">${theGene}<"
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.bed >> mm10.UCSC.knownGenes.unique.DelFirreNotDEhi.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.chr1.bed >> mm10.UCSC.knownGenes.unique.chr1.DelFirreNotDEhi.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.chrX.bed >> mm10.UCSC.knownGenes.unique.chrX.DelFirreNotDEhi.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.diplosomes.bed >> mm10.UCSC.knownGenes.unique.diplosomes.DelFirreNotDEhi.bed
done

wc -l mm10.UCSC.knownGenes.unique.*DelFirreNotDEhi.bed
#    198 mm10.UCSC.knownGenes.unique.DelFirreNotDEhi.bed
#     31 mm10.UCSC.knownGenes.unique.chr1.DelFirreNotDEhi.bed
#      0 mm10.UCSC.knownGenes.unique.chrX.DelFirreNotDEhi.bed
#    198 mm10.UCSC.knownGenes.unique.diplosomes.DelFirreNotDEhi.bed

rm mm10.UCSC.knownGenes.unique.chrX.DelFirreNotDEhi.bed


###########
# Lowly expressed non DE genes

DFNONDEGENESLO=( $( cat "$HOME/proj/2019_RNAseq_analysis/results/gbonora/20190423_RNAseq_DelFirreAndRescue/notes/20190426 metaplots for DE genes/delFirreNonDEgenesLo.txt" ) )
echo "${#DFNONDEGENESLO[@]}"
# 200

cat /dev/null > mm10.UCSC.knownGenes.unique.DelFirreNotDElo.bed
cat /dev/null > mm10.UCSC.knownGenes.unique.chr1.DelFirreNotDElo.bed
cat /dev/null > mm10.UCSC.knownGenes.unique.chrX.DelFirreNotDElo.bed
cat /dev/null > mm10.UCSC.knownGenes.unique.diplosomes.DelFirreNotDElo.bed
for theGene in "${DFNONDEGENESLO[@]}"; do
echo ">${theGene}<"
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.bed >> mm10.UCSC.knownGenes.unique.DelFirreNotDElo.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.chr1.bed >> mm10.UCSC.knownGenes.unique.chr1.DelFirreNotDElo.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.chrX.bed >> mm10.UCSC.knownGenes.unique.chrX.DelFirreNotDElo.bed
awk -v theGene="${theGene}" 'BEGIN{OFS="\t"} $4==theGene{print $0}' mm10.UCSC.knownGenes.unique.diplosomes.bed >> mm10.UCSC.knownGenes.unique.diplosomes.DelFirreNotDElo.bed
done

wc -l mm10.UCSC.knownGenes.unique.*DelFirreNotDElo.bed
#    197 mm10.UCSC.knownGenes.unique.DelFirreNotDElo.bed
#     30 mm10.UCSC.knownGenes.unique.chr1.DelFirreNotDElo.bed
#      0 mm10.UCSC.knownGenes.unique.chrX.DelFirreNotDElo.bed
#    197 mm10.UCSC.knownGenes.unique.diplosomes.DelFirreNotDElo.bed

rm mm10.UCSC.knownGenes.unique.chrX.DelFirreNotDElo.bed



#################################################
# Then rsync to server
rsync -avuPLn  --exclude ".DS_Store"  ~/refdata/mm10/ gs:refdata/mm10/




















##################################################################################################
##################################################################################################
##################################################################################################
# Cast/129 SNPs


#################################################
#################################################
#################################################
# 180908
# Get CAST and 129/sv SNPS
# https://data.4dnucleome.org/biosources/4DNSRMG5APUM/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC148623/pdf/273685.pdf
# F1 ES cell lines 2–2 (female) and 2–3 (male) were prepared by
# standard methods from delayed blastocysts of a cross of a
# Mus musculus casteneus/ei male with a Mus musculus domesticus 129/Sv/Jae female (10).
#
# ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/

mkdir ~/refdata/mm10CAST129SNPs
cd ~/refdata/mm10CAST129SNPs

##########
# 129

wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf.gz
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/129S5SvEvBrd.mgp.v5.snps.dbSNP142.vcf.gz
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/129P2_OlaHsd.mgp.v5.snps.dbSNP142.vcf.gz

# *** USE THIS ONE ****
zcat 129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf.gz | grep -v "^#" | grep -vc PASS
# >0
zcat 129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf.gz | grep PASS | gzip -f > 129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz &
zcat 129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz | grep -vc "^#"
# 5,173,412

zcat 129S5SvEvBrd.mgp.v5.snps.dbSNP142.vcf.gz | grep -v "^#" | grep -vc PASS
# >0
zcat 129S5SvEvBrd.mgp.v5.snps.dbSNP142.vcf.gz | grep PASS | gzip -f > 129S5SvEvBrd.mgp.v5.snps.dbSNP142.PASS.vcf.gz &
zcat 129S5SvEvBrd.mgp.v5.snps.dbSNP142.PASS.vcf.gz | grep -vc "^#"
# 4,697,165

zcat 129P2_OlaHsd.mgp.v5.snps.dbSNP142.vcf.gz | grep -v "^#" | grep -vc PASS
# >0
zcat 129P2_OlaHsd.mgp.v5.snps.dbSNP142.vcf.gz | grep PASS | gzip -f > 129P2_OlaHsd.mgp.v5.snps.dbSNP142.PASS.vcf.gz &
zcat 129P2_OlaHsd.mgp.v5.snps.dbSNP142.PASS.vcf.gz | grep -vc "^#"
# 5,333,316


##########
# CAST

wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz

zcat CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz | grep -v "^#" | grep -vc PASS
# >0
zcat CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz  | grep PASS | gzip -f > CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz  &
zcat CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz | grep -vc "^#"
# 20,668,274

# Compare to Wenxiu's CAST SNPs
zcat $HOME/refdata/mm10pseudoCastBWAindexV2/2016-06-27_CAST-genome-mm10/"CAST.snps.vcf.gz" | grep -vc "^#"
# 17,494,409

zcat $HOME/refdata/mm10pseudoCastBWAindexV2/2016-06-27_CAST-genome-mm10/"CAST.snps.vcf.gz" | grep -v "^#" | head -n 15 | tail -n 10
# 1	3001490	.	C	A	999	PASS	AC=2;AN=2	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:59:0:290,189,0:255,175,0:2:44:59:0,0,25,34:0:-0.693147:.:1
# 1	3001712	.	C	G	999	PASS	AC=2;AN=2	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:41:0:287,136,0:255,123,0:2:57:41:0,0,18,23:0:-0.693145:.:1
# 1	3001745	rs579493097	A	G	999	PASS	AC=2;AN=2	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:117:35:0:285,117,0:255,105,0:2:58:35:0,0,15,20:0:-0.693136:.:1
# 1	3003268	rs30748911	A	G	999	PASS	AC=2;AN=2	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:79:23:0:264,79,0:237,69,0:2:60:23:0,0,22,1:0:-0.692717:.:1
# 1	3003414	rs31953890	A	G	999	PASS	AC=2;AN=2	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:51:0:288,161,0:255,148,0:2:58:51:0,0,34,17:0:-0.693147:.:1
# 1	3003449	rs32186899	T	C	999	PASS	AC=2;AN=2	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:52:0:289,171,0:255,157,0:2:59:52:0,0,27,25:0:-0.693147:.:1
# 1	3003464	rs31079645	G	A	999	PASS	AC=2;AN=2	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:53:0:289,171,0:255,157,0:2:59:53:0,0,28,25:0:-0.693147:.:1
# 1	3003508	rs32044173	C	T	999	PASS	AC=2;AN=2	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:68:0:291,220,0:255,205,0:2:59:68:0,0,39,29:0:-0.693147:.:1
# 1	3003561	rs31477589	T	C	999	PASS	AC=2;AN=2	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:62:0:290,201,0:255,187,0:2:60:62:0,0,34,28:0:-0.693147:.:1
# 1	3004324	rs30948750	T	G	999	PASS	AC=2;AN=2	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:45:0:288,148,0:255,135,0:2:60:45:0,0,18,27:0:-0.693147:.:1

zcat CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz | head -n 16 | tail -n 10
# 1	3001490	.	C	A	228	PASS	DP4=0,0,25,34;DP=59;CSQ=A||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:59:0:290,189,0:255,175,0:2:44:59:0,0,25,34:0:-0.693147:.:1
# 1	3001712	.	C	G	228	PASS	DP4=0,0,18,23;DP=41;CSQ=G||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:41:0:287,136,0:255,123,0:2:57:41:0,0,18,23:0:-0.693145:.:1
# 1	3001745	rs579493097	A	G	228	PASS	DP4=0,0,15,20;DP=35;CSQ=G||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:117:35:0:285,117,0:255,105,0:2:58:35:0,0,15,20:0:-0.693136:.:1
# *** 1	3002064	rs578738275	T	C	173	PASS	DP4=0,0,0,20;DP=20;CSQ=C||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:69:20:0:226,69,0:200,60,0:2:50:20:0,0,0,20:0:-0.692067:.:1
# *** 1	3002104	rs584286359	G	A	171	PASS	DP4=0,0,0,19;DP=19;CSQ=A||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:66:19:0:223,66,0:198,57,0:2:38:19:0,0,0,19:0:-0.69168:.:1
# 1	3003268	rs30748911	A	G	210	PASS	DP4=0,0,22,1;DP=23;CSQ=G||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:79:23:0:264,79,0:237,69,0:2:60:23:0,0,22,1:0:-0.692717:.:1
# 1	3003414	rs31953890	A	G	228	PASS	DP4=0,0,34,17;DP=51;CSQ=G||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:51:0:288,161,0:255,148,0:2:58:51:0,0,34,17:0:-0.693147:.:1
# 1	3003449	rs32186899	T	C	228	PASS	DP4=0,0,27,25;DP=52;CSQ=C||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:52:0:289,171,0:255,157,0:2:59:52:0,0,27,25:0:-0.693147:.:1
# 1	3003464	rs31079645	G	A	228	PASS	DP4=0,0,28,25;DP=53;CSQ=A||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:53:0:289,171,0:255,157,0:2:59:53:0,0,28,25:0:-0.693147:.:1
# 1	3003508	rs32044173	C	T	228	PASS	DP4=0,0,39,29;DP=68;CSQ=T||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:127:68:0:291,220,0:255,205,0:2:59:68:0,0,39,29:0:-0.693147:.:1

# It seems based on first two SNPs with scores < 200 were filtered from Wenxiu's SNP listed though the SNP score has a different scale.
# Notably the DP4 shows that only one strand has SNP coverage
# http://seqanswers.com/forums/archive/index.php/t-10181.html


diff <( zcat $HOME/refdata/mm10CAST129SNPs/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz | grep -v "^#" | awk '{printf("%s_%s\n", $1, $2)}' ) \
 <( zcat $HOME/refdata/mm10pseudoCastBWAindexV2/2016-06-27_CAST-genome-mm10/"CAST.snps.vcf.gz" | grep -v "^#"  | awk '{printf("%s_%s\n", $1, $2)}' ) \
| grep "^<" | sed 's#< ##' > CASTsnpsMissingInWenxiuList.txt &

wc -l CASTsnpsMissingInWenxiuList.txt
# 3,174,015 CASTsnpsMissingInWenxiuList.txt
# ~ 3,173,866 = 20,668,275 - 17,494,409

head CASTsnpsMissingInWenxiuList.txt

awk 'FNR==NR{snpHash[$1]=$1}\
FNR<NR{snpKey=sprintf("%s_%s", $1, $2); if (snpKey in snpHash) {print $0}}' \
<( cat CASTsnpsMissingInWenxiuList.txt | head -n 1000 ) \
<( zcat $HOME/refdata/mm10CAST129SNPs/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz | grep -v "^#" | head -n 2000 ) \
| head -n 10
# 1	3002064	rs578738275	T	C	173	PASS	DP4=0,0,0,20;DP=20;CSQ=C||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:69:20:0:226,69,0:200,60,0:2:50:20:0,0,0,20:0:-0.692067:.:1
# 1	3002104	rs584286359	G	A	171	PASS	DP4=0,0,0,19;DP=19;CSQ=A||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:66:19:0:223,66,0:198,57,0:2:38:19:0,0,0,19:0:-0.69168:.:1
# 1	3004940	rs578487906	C	T	228	PASS	DP4=0,0,25,13;DP=38;CSQ=T||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:126:38:0:286,126,0:255,114,0:2:60:38:0,0,25,13:0:-0.693143:.:1
# 1	3008029	rs582760841	T	C	228	PASS	DP4=0,0,5,15;DP=20;CSQ=C||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:69:20:0:281,69,0:255,60,0:2:60:20:0,0,5,15:0:-0.692067:.:1
# 1	3008068	rs32171215	C	T	228	PASS	DP4=0,0,8,18;DP=26;CSQ=T||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:89:26:0:283,89,0:255,78,0:2:60:26:0,0,8,18:0:-0.692976:.:1
# 1	3010142	rs264965476	T	A	204	PASS	DP4=0,0,3,6;DP=9;CSQ=A||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:33:9:0:250,33,0:231,27,0:2:60:9:0,0,3,6:0:-0.662043:.:1
# 1	3010167	rs32640266	G	T	228	PASS	DP4=0,0,4,7;DP=11;CSQ=T||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:40:11:0:275,40,0:255,33,0:2:60:11:0,0,4,7:0:-0.676189:.:1
# 1	3012078	rs578840132	C	A	185	PASS	DP4=0,0,21,0;DP=21;CSQ=A||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:73:21:0:238,73,0:212,63,0:2:60:21:0,0,21,0:0:-0.692352:.:1
# 1	3012125	rs45877293	T	A	191	PASS	DP4=0,0,26,0;DP=26;CSQ=A||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:89:26:0:246,89,0:218,78,0:2:60:26:0,0,26,0:0:-0.692976:.:1
# 1	3013521	rs252435275	G	A	228	PASS	DP4=0,0,12,22;DP=34;CSQ=A||||intergenic_variant||||||||	GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI	1/1:114:34:0:285,114,0:255,102,0:2:60:34:0,0,12,22:0:-0.693132:.:1

# But SNPs with scores also missing from Wenxiu's list???










#################################################
#################################################
#################################################
# 20190329
# It seems I may have used the incorrect 129 SNPS.
# 
# https://data.4dnucleome.org/biosources/4DNSRMG5APUM/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC148623/pdf/273685.pdf
# F1 ES cell lines 2–2 (female) and 2–3 (male) were prepared by
# standard methods from delayed blastocysts of a cross of a 
# Mus musculus casteneus/ei male with a Mus musculus domesticus 129/Sv/Jae female (10).  
#                                                               ^^^^^^^^^^^^^^^^^^^^^^
# 
# But Sanger MGP used '129S1/SvImJ [JAX #2448]':
# https://media.nature.com/original/nature-assets/nature/journal/v477/n7364/extref/nature10413-s1.pdf
#
# There is a difference:
# http://www.informatics.jax.org/inbred_strains/mouse/docs/129.shtml
# http://www.informatics.jax.org/mgihome/nomen/strain_129.shtml
# 
# 129S1/SvImJ (129/SvImj) 	129S1 	Aw/Aw 	 white (or light)-bellied agouti
# ...
# 129S4/SvJae (129/SvJae) 	129S4 	Aw/Aw 	 white (or light)-bellied agouti
#
# AND:
# # http://www.informatics.jax.org/mgihome/nomen/strain_129.shtml
# Note that 129/Sv (JR000094), which exists as frozen embryos and was discovered to be heterozygous at many loci, 
# is not included on this list as it is clearly not an inbred strain.
#
# Change/Heard allelic atac paper used 129S1/SvImJ:
# https://www.nature.com/articles/ng.3769
#
# Did Gilbert inadvertantly use 129/svImj too?
# https://www.biorxiv.org/content/biorxiv/early/2017/11/21/221762.full.pdf
# 'Reads  of  quality scores  above  30  were independently  mapped  to  129/Sv and castaneus reference  genomes'
#
# Try to get 129/Sv/Jae from dbSNP archive,
# But archived:
# https://ncbiinsights.ncbi.nlm.nih.gov/2017/05/09/phasing-out-support-for-non-human-genome-organism-data-in-dbsnp-and-dbvar/
# Now at:
# https://www.ebi.ac.uk/eva/?Home
# https://www.ebi.ac.uk/ena/data/view/PRJEB11471

cd ~/refdata/mm10CAST129SNPs

# Archived NCBI dbSNP vcf
# https://ncbiinsights.ncbi.nlm.nih.gov/2017/05/09/phasing-out-support-for-non-human-genome-organism-data-in-dbsnp-and-dbvar/
# https://www.ncbi.nlm.nih.gov/mailman/pipermail/dbsnp-announce/2018q2/000186.html
wget ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz

# No apparent strain information:
zcat 00-All.vcf.gz

# ftp://ftp-mouse.sanger.ac.uk/current_snps/README
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz

# No apparent svJae SNPs
zcat mgp.v5.merged.snps_all.dbSNP142.vcf.gz | grep SvJae | head

# https://www.ebi.ac.uk/ena/data/view/PRJEB11471
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ126/ERZ126133/mgp.v5.merged.snps_all.dbSNP142.noALT_X.done.vcf.gz

# What's in this?
zcat mgp.v5.merged.snps_all.dbSNP142.noALT_X.done.vcf.gz | more








#################################################
#################################################
#################################################
# 20181115
# Snps & indels density w/in windows for 129 and Cast

mkdir ~/refdata/mm10_binnedVariantDensity_129andCast
cd ~/refdata/mm10_binnedVariantDensity_129andCast

for assembly in mm10; do
$HOME/bin/UCSC.utilities/fetchChromSizes $assembly | grep -v "_" | sort -k1V | grep -v "chrM" > $assembly.sizes
done

BINSIZE=500000

SNPDIR=$HOME/refdata/mm10CAST129SNPs
unset SNPFILES
declare -A SNPFILES
SNPFILES["129"]=${SNPDIR}/129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz
SNPFILES["Cast"]=${SNPDIR}/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz
echo ${SNPFILES[@]}

#####
#####
# Binned SNP density

module load gawk/4.1.4

for SNPLINE in 129 Cast; do
# SNPLINE=129
# SNPLINE=Cast

outFile=binnedSNPDensity_${SNPLINE}_${BINSIZE}.bedGraph

awk -v binSize=${BINSIZE} '
FNR==NR{ chroms[NR]=$1; chrLength[$1]=$2 }
FNR==1 && NR>1 {
  nChroms = NR - 1;
  for (i=1; i<=nChroms; i++) {
    snpChrom = chroms[i];
    for (j=0; j<=chrLength[snpChrom]; j+=binSize) {
      snpBin = sprintf("%d", j/binSize);
      # print i, snpChrom, chrLength[snpChrom], j, snpBin;
      binnedSNPcounts[snpChrom][snpBin]=0;
    }
  }
}
FNR<NR{
  snpChrom=sprintf("chr%s", $1);
  snpBin=sprintf("%d", $2/binSize);
  # print $1, snpChrom, $2, snpBin;
  binnedSNPcounts[snpChrom][snpBin]++;
}
END{
  # print nChroms;
  for (i=1; i<=nChroms; i++) {
    snpChrom=chroms[i];
    # print i, snpChrom, chrLength[snpChrom];
    for (j=0; j<=chrLength[snpChrom]; j+=binSize) {
      snpBin = sprintf("%d", j/binSize);
      # print $1, snpChrom, $2, snpBin;
      printf("%s\t%d\t%d\t%d\n", snpChrom, snpBin*binSize, (snpBin+1)*binSize, binnedSNPcounts[snpChrom][snpBin]);
    }
  }
}' \
<( cat $assembly.sizes ) <( zcat ${SNPFILES[${SNPLINE}]} | grep -v "^#" ) > ${outFile}

done

ls -ltrh

#####
# Aggregate SNP densities

paste binnedSNPDensity_129_500000.bedGraph <( cat binnedSNPDensity_Cast_500000.bedGraph | cut -f 4) > binnedSNPDensity_129andCast_500000.tsv

more $HOME/refdata/mm10_binnedVariantDensity_129andCast/binnedSNPDensity_129andCast_500000.tsv





# #####
# #####
# # Per chrom Indels
#
# zcat $HOME/refdata/mm10pseudoSpretus.dbSNP142/SPRET_EiJ.mgp.v5.indels.dbSNP142.normed.vcf.gz | grep -v "^#" | grep PASS | cut -f 1 | sort | uniq -c | awk '{printf("chr%s\t%d\n", $2, $1)}' | sort -k1V > spretus.dbSNP142.indelsPerChr.txt
#
# zcat $HOME/refdata/mm10pseudoSpretus/spretus.indels.vcf.gz | grep -v "^#" | grep PASS | cut -f 1 | sort | uniq -c | awk '{printf("chr%s\t%d\n", $2, $1)}' | sort -k1V > spretus.extracted.indelsPerChr.txt
#
# more *.indelsPerChr.txt










#################################################
#################################################
#################################################
# 20190402
# Cast 129 SNP overlap analysis
#
# Where SNPs overlap and alt is the same for both Cast and 129,
# reads from both Cast and 129 alleles will be assigned as alt by samtools and allelic ratios will be 1.
#
# Where SNPs overlap and alt is different for Cast and 129,
# only reads from corresponding allele will be considered and allelic ratios will still be 1.

SNPDIR=$HOME/refdata/mm10CAST129SNPs

awk 'BEGIN{OFS="\t"} \
NR==FNR{ snphash[sprintf("%s_%s", $1, $2)]=sprintf("%s_%s\t%s%s", $1, $2, $4, $5) } \
NR>FNR{ key=sprintf("%s_%s", $1, $2); if (key in snphash) {printf("%s\t%s_%s\t%s%s\n", snphash[key], $1, $2, $4, $5)} }' \
<( zcat $SNPDIR/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz | grep -v "^#" | grep "PASS" ) \
<( zcat $SNPDIR/129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz | grep -v "^#" | grep "PASS" ) \
> $SNPDIR/Cast129SNPoverlap.txt

# head  $SNPDIR/Cast129SNPoverlap.txt

cat $SNPDIR/Cast129SNPoverlap.txt | awk '$2!=$4{ print $0 }' > $SNPDIR/Cast129SNPoverlapDiffSNP.txt

nCast=$( zcat $SNPDIR/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz | grep -v "^#" | grep -c "PASS" )
n129=$( zcat $SNPDIR/129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz | grep -v "^#" | grep -c "PASS" )
nOverlap=$( cat $SNPDIR/Cast129SNPoverlap.txt | wc -l )
nDiffOverlap=$( cat $SNPDIR/Cast129SNPoverlapDiffSNP.txt | wc -l )

echo ${nCast} ${n129} ${nOverlap} ${nDiffOverlap}
# 20668274 5173412 2452473 27254

echo ${nCast} ${n129} ${nOverlap} ${nDiffOverlap} | \
awk 'BEGIN{ printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", \
"nCast", "n129", "nOverlap", "%Cast", "%129", "nOverlapDiffSNP", "%CastOverlapDiffSNP", "%129OverlapDiffSNP" ) }; \
{ printf("%d\t%d\t%d\t%0.2f\t%0.2f\t%d\t%0.2f\t%0.2f\n", \
$1, $2, $3, $3/$1*100, $3/$2*100, $4, $4/$1*100, $4/$2*100 ) }' > $SNPDIR/Cast129SNPoverlapSummary.txt

cat $SNPDIR/Cast129SNPoverlapSummary.txt
# nCast   n129    nOverlap        %CastOverlap   %129Overlap    nOverlapDiffSNP    %CastOverlapDiffSNP   %129OverlapDiffSNP
# 20,668,274        5,173,412 2,452,473 11.87   47.41   27,254   0.13    0.53










#################################################
#################################################
#################################################
# Eliminate all overlapping Cast 129 SNPs

# Cast
awk 'BEGIN{OFS="\t"} \
NR==FNR{ snphash[$1]=$1 } \
NR>FNR{ key=sprintf("%s_%s", $1, $2); if (!(key in snphash)) {print $0} }' \
<( cat $SNPDIR/Cast129SNPoverlap.txt ) \
<( zcat $SNPDIR/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz ) \
| gzip -f > $SNPDIR/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.No129.vcf.gz &

zcat $SNPDIR/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.No129.vcf.gz | grep -v "^#" | grep -c "PASS"
#  20,668,274 - 2,452,473 = 18,215,801


# 129
awk 'BEGIN{OFS="\t"} \
NR==FNR{ snphash[$1]=$1 } \
NR>FNR{ key=sprintf("%s_%s", $1, $2); if (!(key in snphash)) {print $0} }' \
<( cat $SNPDIR/Cast129SNPoverlap.txt ) \
<( zcat $SNPDIR/129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz ) \
| gzip -f > $SNPDIR/129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.NoCast.vcf.gz &

zcat $SNPDIR/129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.NoCast.vcf.gz | grep -v "^#" | grep -c "PASS"
# 5,173,412 - 2,452,473 = 2,720,939

# 20,936,740 total








#################################################
#################################################
#################################################
# 20190508
# vcf from Nicolas Servant
# 
# On May 6, 2019, at 4:08 AM, Servant Nicolas <Nicolas.Servant@curie.fr> wrote:
# I uploaded our vcf file for Cast/129 here ;
# http://zerkalo.curie.fr/partage/temp/snps_CAST_129S1.vcf
# 
# On May 7, 2019, at 1:15 AM, Servant Nicolas <Nicolas.Servant@curie.fr> wrote:
# This is an old server, but still secured ;)
# I added a login/pass
# guest/g0est1!
# If you still have the message, you should be able to force firefox to access the server.
#
# Note: Servant SNP file is 129/Cast

cd ~/refdata/mm10CAST129SNPs
SNPDIR=$HOME/refdata/mm10CAST129SNPs

wget https://zerkalo.curie.fr/partage/temp/snps_CAST_129S1.vcf --no-check-certificate --http-user=guest --http-password=g0est1!

gzip snps_CAST_129S1.vcf

zcat $SNPDIR/snps_CAST_129S1.vcf.gz | grep -v "^#"
# 20,563,466

awk 'BEGIN{OFS="\t"} \
NR==FNR{ snphash[sprintf("%s_%s", $1, $2)]=sprintf("%s_%s\t%s%s", $1, $2, $4, $5) } \
NR>FNR{ key=sprintf("chr%s_%s", $1, $2); if (key in snphash) {printf("%s\t%s_%s\t%s%s\n", snphash[key], $1, $2, $4, $5)} }' \
<( zcat $SNPDIR/snps_CAST_129S1.vcf.gz | grep -v "^#" | grep "PASS" ) \
<( zcat $SNPDIR/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.No129.vcf.gz | grep -v "^#" | grep "PASS" ) \
> $SNPDIR/servantCastSNPoverlap.txt

awk 'BEGIN{OFS="\t"} \
NR==FNR{ snphash[sprintf("%s_%s", $1, $2)]=sprintf("%s_%s\t%s%s", $1, $2, $4, $5) } \
NR>FNR{ key=sprintf("chr%s_%s", $1, $2); if (key in snphash) {printf("%s\t%s_%s\t%s%s\n", snphash[key], $1, $2, $4, $5)} }' \
<( zcat $SNPDIR/snps_CAST_129S1.vcf.gz | grep -v "^#" | grep "PASS" ) \
<( zcat $SNPDIR/129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.NoCast.vcf.gz | grep -v "^#" | grep "PASS" ) \
> $SNPDIR/servant129SNPoverlap.txt

# Must flip the order b/c Servant SNP file is 129/Cast
cat $SNPDIR/servantCastSNPoverlap.txt | \
awk '{ a=substr($4,1,1); b=substr($4,2,1); ba=sprintf("%s%s", b, a); if($2!=ba){ print $0 } }' \
> $SNPDIR/servantCastSNPoverlapDiffSNP.txt

cat $SNPDIR/servant129SNPoverlap.txt | \
awk '$2!=$4{ print $0 }' > $SNPDIR/servant129SNPoverlapDiffSNP.txt

nServant=$( zcat $SNPDIR/snps_CAST_129S1.vcf.gz | grep -v "^#" | grep -c "PASS" )
nCast=$( zcat $SNPDIR/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.No129.vcf.gz | grep -v "^#" | grep -c "PASS" )
n129=$( zcat $SNPDIR/129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.NoCast.vcf.gz | grep -v "^#" | grep -c "PASS" )
nCastOverlap=$( cat $SNPDIR/servantCastSNPoverlap.txt | wc -l )
nCastDiffOverlap=$( cat $SNPDIR/servantCastSNPoverlapDiffSNP.txt | wc -l )
n129Overlap=$( cat $SNPDIR/servant129SNPoverlap.txt | wc -l )
n129DiffOverlap=$( cat $SNPDIR/servant129SNPoverlapDiffSNP.txt | wc -l )
echo ${nServant} ${nCast} ${nCastOverlap} ${nCastDiffOverlap} ${n129} ${n129Overlap} ${n129DiffOverlap}
# 17,516,347 18,215,801 15,449,832 0 2,720,939 2,044,815 0

echo ${nServant} ${nCast} ${nCastOverlap} ${nCastDiffOverlap} ${n129} ${n129Overlap} ${n129DiffOverlap} | \
awk 'BEGIN{ printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", \
"nServant", \
"nCast", "nCastOverlap", "%CastOverlap", "nCastOverlapDiffSNP", "%CastOverlapDiffSNP", \
"n129", "n129Overlap", "%129Overlap", "n129DiffOverlapDiffSNP", "%129OverlapDiffSNP" ) }; \
{ printf("%d\t%d\t%d\t%0.2f\t%d\t%0.2f\t%d\t%d\t%0.2f\t%d\t%0.2f\n", \
$1, \
$2, $3, $3/$2*100, $4, $4/$2*100, \
$5, $6, $6/$5*100, $7, $7/$5*100 ) }' > $SNPDIR/servantCast129SNPoverlapSummary.txt

cat $SNPDIR/servantCast129SNPoverlapSummary.txt
# nServant        nCast   nCastOverlap    %CastOverlap    nCastOverlapDiffSNP     %CastOverlapDiffSNP     n129    n129Overlap    %129Overlap      n129DiffOverlapDiffSNP  %129OverlapDiffSNP
# 17,516,347        18,215,801        15,449,832        84.82   0       0.00    2,720,939 2,044,815 75.15   0       0.00







#################################################
# Get SNPs that overlap Servant SNPs

awk 'BEGIN{OFS="\t"} \
NR==FNR{ snphash[sprintf("%s_%s", $1, $2)]=sprintf("%s_%s\t%s%s", $1, $2, $4, $5) } \
NR>FNR{ key=sprintf("chr%s_%s", $1, $2); if (key in snphash) {print $0} }' \
<( zcat $SNPDIR/snps_CAST_129S1.vcf.gz | grep -v "^#" | grep "PASS" ) \
<( zcat $SNPDIR/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.No129.vcf.gz | grep -v "^#" | grep "PASS" ) \
| gzip -f > $SNPDIR/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.No129.servant.vcf.gz

zcat $SNPDIR/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.No129.servant.vcf.gz | wc -l 
# 15,449,832 

awk 'BEGIN{OFS="\t"} \
NR==FNR{ snphash[sprintf("%s_%s", $1, $2)]=sprintf("%s_%s\t%s%s", $1, $2, $4, $5) } \
NR>FNR{ key=sprintf("chr%s_%s", $1, $2); if (key in snphash) {print $0} }' \
<( zcat $SNPDIR/snps_CAST_129S1.vcf.gz | grep -v "^#" | grep "PASS" ) \
<( zcat $SNPDIR/129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.NoCast.vcf.gz | grep -v "^#" | grep "PASS" ) \
| gzip -f > $SNPDIR/129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.NoCast.servant.vcf.gz

zcat $SNPDIR/129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.NoCast.servant.vcf.gz | wc -l
# 2,044,815








#################################################
#################################################
#################################################
# 20190515
# Filter vcf based on FI field as Nicolas did:
# ##FORMAT=<ID=FI,Number=1,Type=Integer,Description="Whether a sample was a Pass(1) or fail (0) based on FILTER values">

zcat CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz  | grep ":1$" | gzip -f > CAST_EiJ.mgp.v5.snps.dbSNP142.FI.vcf.gz  &
zcat 129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf.gz | grep ":1$" | gzip -f > 129S1_SvImJ.mgp.v5.snps.dbSNP142.FI.vcf.gz &

zcat CAST_EiJ.mgp.v5.snps.dbSNP142.FI.vcf.gz | grep -vc "^#"
# 20,668,274
# 20,668,274

zcat CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz | grep -vc "^#"
# 20,668,274

zcat 129S1_SvImJ.mgp.v5.snps.dbSNP142.FI.vcf.gz | grep -vc "^#"
# 5,173,412

zcat 129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz | grep -vc "^#"
# 5,173,412


##########
# Do PASS and FI always coincide?

zcat CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz | grep ":1$" | grep "PASS" | wc -l 
# 20,668,274

zcat 129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf.gz | grep ":1$" | grep "PASS" | wc -l 
# 5,173,412

diff <( zcat CAST_EiJ.mgp.v5.snps.dbSNP142.FI.vcf.gz ) <( zcat CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz ) > CAST_EiJ.FIvsPASS.diff.txt &

diff <( zcat 129S1_SvImJ.mgp.v5.snps.dbSNP142.FI.vcf.gz ) <( zcat 129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.vcf.gz ) > 129S1_SvImJ.FIvsPASS.diff.txt &

# Yes, yes they do!




















##################################################################################################
##################################################################################################
##################################################################################################
# Random genome selections


#################################################
#################################################
#################################################
# n random sites in genome
# https://www.biostars.org/p/48617/

cd ~/refdata/mm10

module load bedtools/latest

bedtools random -n 100000 -l 1 -g chromInfo.clean.txt > mm10.100000loci.txt










#################################################
#################################################
#################################################
# n random bases in genome
# https://www.biostars.org/p/48617/
# https://sourceforge.net/projects/bbmap/files/

cd ~/refdata/mm10

module load java/latest
module load BBMap/37.66

which randomreads.sh
# /net/gs/vol3/software/modules-sw/BBMap/37.66/Linux/RHEL6/x86_64/randomreads.sh
more /net/gs/vol3/software/modules-sw/BBMap/37.66/Linux/RHEL6/x86_64/randomreads.sh

randomreads.sh -Xmx8g ref=mm10.fa out=mm10.10M10bpReads.fq len=10 reads=10000000 simplenames=t banns=t # samestrand=t

cat mm10.10Mk10bpReads.fq | head
# @0_+_40625053_40625062_chr17
# ATTCGCAAGC
# +
# A??@=<A?AA
# @1_+_77806904_77806913_chr3
# TGGCCATCTA
# +
# BAAA@B>BC=

cat mm10.10Mk10bpReads.fq | grep -A 1 "chr1$" | grep -A 1 _+_ | more
cat mm10.10M10bpReads.fq | grep -A 1 "chr1$" | grep -A 1 _-_ | more

# *** NOTE: APPEARS TO GIVE ZERO BASED POSITIONS ***
# ****      AND REVERSE COMPLEMENTS - STRAND SEQUENCES ***

# So for a list of random bases in genome just take first base from + strand seqs & increment position
# SNP loci given in 1-base position

cat mm10.10M10bpReads.fq | grep -A 1 "_+_" | grep -ve "--" | paste -d '_' - - | grep -v "chr[YM]" | cut -d '_' -f 5,3,6 | tr '_' '\t' | awk 'BEGIN{OFS="\t"}{gsub("chr", "", $2); print $2, $1+1, substr($3,1,1)}' | head -n 1000000 > mm10.1M.randomBases.txt

gzip -f mm10.1M.randomBases.txt

zcat  mm10.1M.randomBases.txt.gz | wc -l
# 1000000

# check:
cat mm10.10M10bpReads.fq | grep -A 1 "_+_" | grep -ve "--" | paste -d '_' - - | grep -v "chr[YM]" | cut -d '_' -f 5,3,6 | tr '_' '\t' | head
# 40625053        chr17   ATTCGCAAGC
# 77806904        chr3    TGGCCATCTA
# 40545755        chr8    TTCCGCTTTC
# 50351022        chr2    AGGAGTGTAT
# 50618433        chr16   CTATGTTGAA
# 108379433       chr8    CCTCCCTTGC
# 83824827        chr12   AATTTTCTTT
# 10412737        chr3    TGTCCTTTTG
# 104742946       chr1    TCCAGTGTCT
# 16130438        chr8    TAAATAAATA

zcat  mm10.1M.randomBases.txt.gz | head
# 17      40625054        A
# 3       77806905        T
# 8       40545756        T
# 2       50351023        A
# 16      50618434        C
# 8       108379434       C
# 12      83824828        A
# 3       10412738        T
# 1       104742947       T
# 8       16130439        T


cat mm10.10M10bpReads.fq | grep -A 1 "_+_" | grep -ve "--" | paste -d '_' - - | grep -v "chr[YM]" | cut -d '_' -f 5,3,6 | tr '_' '\t' | awk 'BEGIN{OFS="\t"}{gsub("chr", "", $2); print $2, $1+1, substr($3,1,1)}' | head -n 4000000 > mm10.4M.randomBases.txt

gzip -f mm10.4M.randomBases.txt

zcat  mm10.4M.randomBases.txt.gz | wc -l
# 4000000

# There are some non-unique reads/bases:
zcat  mm10.4M.randomBases.txt.gz | awk '{printf("%s_%s\n", $1, $2)}' | sort | uniq | wc -l
# 3996866




















##################################################################################################
##################################################################################################
##################################################################################################
# 20191105 N-masked spretus assemblies


#################################################
#################################################
#################################################
# 20191105
# N-masked mm10 assemblies using validated Patski spretus SNPs

mkdir ~/refdata/mm10.Nmasked.spretus
cd ~/refdata/mm10.Nmasked.spretus

~/bin/SNPsplit_v0.3.2/SNPsplit_genome_preparation --help | more
# --skip_filtering              This option skips reading and filtering the VCF file. This assumes that a folder named
#                               'SNPs_<Strain_Name>' exists in the working directory, and that text files with SNP information
#                               are contained therein in the following format:
#                                           SNP-ID     Chromosome  Position    Strand   Ref/SNP
#                               example:   33941939        9       68878541       1       T/G

zcat $HOME/refdata/mm10pseudoSpretus/mgp.v5.merged.snps_all.dbSNP142.vcf.gz | more
# Note: uses Ensemble notation i.e. no 'chr' in front of chromosome #.

# Validated spretus SNPs:
zcat $HOME/refdata/mm10pseudoSpretus/spretus.patski.snps.vcf.gz | head -n 2
# 1       3000185 rs585444580     G       T       999     PASS    AC=2;AN=2       GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI  1/1:114:34:0:285,114,0:255,102,0:2:57:34:0,0,32,2:0:-0.693132:.:1
# 1       3000234 rs579469519     G       A       999     PASS    AC=2;AN=2       GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI  1/1:127:47:0:288,154,0:255,141,0:2:54:47:0,0,41,6:0:-0.693147:.:1

# Generate filtered VCF locally:
# Question: Strand is 1 in all cases right?
#
# Must use strain names in VCF files, even if VCF file is not actually used:
# 'This assumes that a folder named 'SNPs_<Strain_Name>' exists in the working directory'
mkdir ./SNPs_SPRET_EiJ 

# zcat $HOME/refdata/mm10pseudoSpretus/spretus.patski.snps.vcf.gz | \
# awk 'BEGIN{FS="\t"; OFS="\t"}{printf("%s\t%s\t%s\t1\t%s/%s\n", $3, $1, $2, $4, $5)}' > \
# ./SNPs_SPRET_EiJ/spretus.patski.snps.SNPsplitFormat.vcf 
# # Can be named anything apparently as long as in ./SNPs_SPRET_EiJ subfolder <-- *** UH, WRONG!!! ***

# 2019117
zcat $HOME/refdata/mm10pseudoSpretus/spretus.patski.snps.vcf.gz | \
awk -v outfilePrefix="./SNPs_SPRET_EiJ" \
'BEGIN{FS="\t"; OFS="\t"} \
{outfile=sprintf("%s/chr%s.txt", outfilePrefix, $1); printf("%s\t%s\t%s\t1\t%s/%s\n", $3, $1, $2, $4, $5) > outfile}'

ls -ltrh ./SNPs_SPRET_EiJ/

head ./SNPs_SPRET_EiJ/chr1.txt
# rs585444580     1    3000185 1       G/T
# rs579469519     1    3000234 1       G/A
# ..

head ./SNPs_SPRET_EiJ/chr2.txt
# rs585041348     2       3051194 1       G/A
# rs583767423     2       3051213 1       T/G
# ..


# -----------------------------------------------
# REFERENCE ASSEMBLY ISSUE:

# Must ensure that the same version of the genome is used for both VCF annotations and reference genome:
# 'The chromosome name given in the VCF file was '11' and was not found in the reference genome.
#  A rather common mistake might be that the VCF file was downloaded from Ensembl (who use chromosome names such as 1, 2, X, MT)
#  but the genome from UCSC (who use chromosome names such as chr1, chr2, chrX, chrM)'

# This is for sciOmics processing. Check scihic reference assembly
ls -ltrh /net/noble/vol2/home/gbonora/refdata/hg19mm10
# lrwxrwxrwx 1 gbonora noblelab 41 Aug  6 13:58 /net/noble/vol2/home/gbonora/refdata/hg19mm10 -> /net/noble/vol2/home/gurkan/Data/hg19mm10

ls -ltrh /net/noble/vol2/home/gbonora/refdata/hg19mm10/
# total 12G
# -rw-r--r-- 1 gurkan gurkan_g 2.6G Oct 26  2017 mm10.anno.fa
# -rw-r--r-- 1 gurkan gurkan_g 3.0G Oct 26  2017 hg19.anno.fa
# -rw-r--r-- 1 gurkan gurkan_g 5.6G Oct 26  2017 mm10hg19.combo.fa
# drwxr-sr-x 2 gurkan gurkan_g 4.0K Oct 26  2017 bowtie2indexbroken
# drwxr-sr-x 2 gurkan gurkan_g 4.0K Nov  8  2017 bowtie2-index

cat /net/noble/vol2/home/gbonora/refdata/hg19mm10/mm10hg19.combo.fa | grep "^>" | grep -v "_random" | grep -v "chrUn"
# >human_chr1
# ..
# >human_chr21
# >human_chrM
# >chr1
# ..
# >chrY
#
# So scihic pipeline uses UCSC reference w/ 'chr' before chromosome #s

# NOTE: THIS IS NOT THE FASAT THAT WAS USED in the scihic combo assembly
cat /net/noble/vol2/home/gbonora/refdata/hg19mm10/mm10.anno.fa | grep "^>"
# >mouse_chr1
# >mouse_chr10
# ...

cat /net/noble/vol2/home/gbonora/refdata/iGenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa | grep "^>"
# >chr10
# ..
# >chrM
# >chrX
# >chrY

cat /net/noble/vol2/home/gbonora/refdata/iGenomes/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa | grep "^>"
# >10
# 
# >MT
# >X
# >Y

# So use Ensemble assembly to generate N-masked assembly then add 'chr' to chrom names.

# -----------------------------------------------

refFile=/net/noble/vol2/home/gbonora/refdata/iGenomes/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta
altSnpFile=$HOME/refdata/mm10pseudoSpretus/mgp.v5.merged.snps_all.dbSNP142.vcf.gz

# 20191117 Rerun to use correct chr-by-chr SNP files in SNPs_SPRET_EiJ folder
~/bin/SNPsplit_v0.3.2/SNPsplit_genome_preparation \
--nmasking \
--skip_filtering \
--vcf_file ${altSnpFile} \
--reference_genome ${refFile} \
--strain SPRET_EiJ \
--genome_build mm10.Nmasked.spretus 2>&1 | tee SNPsplit_genome_preparation.out

# 2019117 check:
more SPRET_EiJ_genome_preparation_report.txt
# ..
# Processing chromosome 11 (for strain SPRET_EiJ)
# Reading SNPs from file /net/noble/vol4/noble/user/gbonora/refdata/mm10.Nmasked.spretus/SNPs_SPRET_EiJ/chr11.txt
# Writing modified chromosome (N-masking)
# Writing N-masked output to: /net/noble/vol4/noble/user/gbonora/refdata/mm10.Nmasked.spretus/SPRET_EiJ_N-masked/chr11.N-masked.fa
# 1417348 SNPs total for chromosome 11
# 1417348 positions on chromosome 11 were changed to 'N'
# ..
# Summary
# 29159465 Ns were newly introduced into the N-masked genome for strain SPRET_EiJ in total
# ..

ls -ltrh ~/refdata/mm10.Nmasked.spretus/SPRET_EiJ_N-masked
# total 2.6G
# -rw-rw-r-- 1 gbonora noblelab 118M Nov 17 14:04 chr11.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 141M Nov 17 14:04 chr7.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  89M Nov 17 14:04 chrY.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 176M Nov 17 14:04 chr2.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  92M Nov 17 14:04 chr17.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 189M Nov 17 14:04 chr1.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  88M Nov 17 14:04 chr18.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 116M Nov 17 14:04 chr13.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  95M Nov 17 14:04 chr16.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 145M Nov 17 14:05 chr6.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 155M Nov 17 14:05 chr3.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 165M Nov 17 14:05 chrX.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  17K Nov 17 14:05 chrMT.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 121M Nov 17 14:05 chr9.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 116M Nov 17 14:05 chr12.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 121M Nov 17 14:05 chr14.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 101M Nov 17 14:05 chr15.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 125M Nov 17 14:05 chr8.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 151M Nov 17 14:06 chr4.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 126M Nov 17 14:06 chr10.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  60M Nov 17 14:06 chr19.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 147M Nov 17 14:06 chr5.N-masked.fa



#################################################
#################################################
#################################################
# 20191117 Noticed this:

cd ~/refdata/mm10.Nmasked.spretus
more SNPsplit_genome_preparation.out
# ...
# Processing chromosome 11 (for strain SPRET_EiJ)
# Couldn't find SNP file for chromosome '11' '/net/noble/vol4/noble/user/gbonora/refdata/mm10.Nmasked.spretus/SNPs_SPRET_EiJ/chr11
# .txt' didn't exist. Skipping...
# Clearing SNP array...
# Writing modified chromosome (N-masking)
# Writing N-masked output to: /net/noble/vol4/noble/user/gbonora/refdata/mm10.Nmasked.spretus/SPRET_EiJ_N-masked/chr11.N-masked.fa
# 0 SNPs total for chromosome 11
# 0 positions on chromosome 11 were changed to 'N'
# ...

# Cast/129 SNPs stored in a chr-by-chr files (see below):
# *** MUST ENSURE THIS IS THE CASE FOR THE VALIDATE SPRETUS SNPS!!! ***

ls -ltrh SNPs_SPRET_EiJ
# total 813M
# -rw-rw-r-- 1 gbonora noblelab 813M Nov  6 16:03 spretus.patski.snps.SNPsplitFormat.vcf

# Had to fix this... see 2019117 work above.









#################################################
#################################################
#################################################
# 20191105
# Cat new hg19 + Nmasked spretus mm10 combo assembly 
# for scihic pipeline and bowtie2 index

# NOTE: must add 'chr' to mouse chromosome names

mkdir $HOME/refdata/hg19mm10.Nmasked.spretus
cd $HOME/refdata/hg19mm10.Nmasked.spretus

hg19SRCDIR="$HOME/refdata/hg19mm10"
mm10SRCDIR="$HOME/refdata/mm10.Nmasked.spretus/SPRET_EiJ_N-masked"

# Combo fasta
# Must add 'chr' to chrom names
# cat ${mm10SRCDIR}/*.fa | sed 's/^>/>chr/' | sed 's/^>chrMT/>chrM/' | grep "^>"
cat ${hg19SRCDIR}/hg19.anno.fa <( cat ${mm10SRCDIR}/*.fa | sed 's/^>/>chr/' | sed 's/^>chrMT/>chrM/' ) > hg19.mm10.Nmasked.spretus.combo.fa
cat hg19.mm10.Nmasked.spretus.combo.fa | grep "^>"

# Now build bowtie2 index
module load bowtie2/2.2.3
bowtie2-build hg19.mm10.Nmasked.spretus.combo.fa mm10hg19combo

mkdir ./bowtie2-index
mv mm10hg19combo* bowtie2-index

bowtie2-inspect -n ./bowtie2-index/mm10hg19combo
# ..
# human_chr18_gl000207_random
# chr10
# chr11
# chr12
# chr13
# chr14
# chr15
# chr16
# chr17
# chr18
# chr19
# chr1
# chr2
# chr3
# chr4
# chr5
# chr6
# chr7
# chr8
# chr9
# chrM
# chrX
# chrY










#################################################
#################################################
#################################################
# 20191108
# Nmasked spretus assembly for use with sciRNA-seq STAR alignment 

# Look at Jonathan Packer's assembly:
# STAR_INDEX=/net/trapnell/vol1/jspacker/STAR/mouse/GRCm38-primary-assembly

ls -ltrh /net/trapnell/vol1/jspacker/STAR/

# Already copied locally:
ls -ltrh  ~/refdata/STAR

more  ~/refdata/STAR/STAR-genomeGenerate.sh
# #$ -S /bin/bash
# #$ -l mfree=8G
# #$ -l h_rt=3:0:0:0
# #$ -pe serial 12
# 
# module load STAR/2.5.2b
# 
# FASTA=$1
# GTF=$2
# OUTPUT=$3
# WORKING_DIR=$4
# 
# cd $WORKING_DIR
# 
# STAR \
#     --runThreadN 12 \
#     --runMode genomeGenerate \
#     --genomeDir $OUTPUT \
#     --genomeFastaFiles $FASTA \
#     --sjdbGTFfile $GTF \
#     --outTmpDir $WORKING_DIR/tmp


ls -ltrh  ~/refdata/STAR/mouse
zcat ~/refdata/STAR/mouse/GRCm38.primary_assembly.genome.fa.gz | grep -e '^>'
# >chr1 1
# >chr2 2
# >chr3 3
# ..
# >chrX X
# >chrY Y
# >chrM MT
# >GL456210.1 GL456210.1
# ...


# Replace Jonathan's assembly with N-masked one.
# See above for how this was generated:
cd ~/refdata/STAR/mouse

mm10SRCDIR="$HOME/refdata/mm10.Nmasked.spretus/SPRET_EiJ_N-masked"
# cat ${mm10SRCDIR}/*.fa | sed 's/^>\(.*\)/>chr\1 \1/' | sed 's/^>chrMT/>chrM/' | grep "^>"
cat ${mm10SRCDIR}/*.fa | sed 's/^>\(.*\)/>chr\1 \1/' | sed 's/^>chrMT/>chrM/' > ~/refdata/STAR/mouse/GRCm38-primary-assembly.SPRET_EiJ_N-masked.fa

mkdir ~/refdata/STAR/mouse/GRCm38-primary-assembly.SPRET_EiJ_N-masked

# cp ~/refdata/STAR/STAR-genomeGenerate.sh  ~/refdata/STAR/STAR-genomeGenerate.GB.sh 
vi ~/refdata/STAR/STAR-genomeGenerate.GB.sh 

cd ~/refdata/STAR
rm -rf mouse/tmp
qsub -j y -cwd ./STAR-genomeGenerate.GB.sh \
~/refdata/STAR/mouse/GRCm38-primary-assembly.SPRET_EiJ_N-masked.fa \
~/refdata/STAR/mouse/gencode.vM15.primary_assembly.annotation.gtf \
~/refdata/STAR/mouse/GRCm38-primary-assembly.SPRET_EiJ_N-masked \
~/refdata/STAR/mouse

mv mouse/Log.out mouse/Log.spretus.out

ls -ltrh ~/refdata/STAR/mouse/GRCm38-primary-assembly.SPRET_EiJ_N-masked
# Using correctly masked assembly:
# -rw-rw-r-- 1 gbonora noblelab  932 Nov 18 00:42 genomeParameters.txt
# -rw-rw-r-- 1 gbonora noblelab  120 Nov 18 00:44 chrName.txt
# -rw-rw-r-- 1 gbonora noblelab  211 Nov 18 00:44 chrLength.txt
# -rw-rw-r-- 1 gbonora noblelab  236 Nov 18 00:44 chrStart.txt
# -rw-rw-r-- 1 gbonora noblelab  331 Nov 18 00:44 chrNameLength.txt
# -rw-rw-r-- 1 gbonora noblelab  27M Nov 18 01:23 exonGeTrInfo.tab
# -rw-rw-r-- 1 gbonora noblelab 1.1M Nov 18 01:23 geneInfo.tab
# -rw-rw-r-- 1 gbonora noblelab 8.0M Nov 18 01:23 transcriptInfo.tab
# -rw-rw-r-- 1 gbonora noblelab  11M Nov 18 01:23 exonInfo.tab
# -rw-rw-r-- 1 gbonora noblelab 6.8M Nov 18 01:23 sjdbList.fromGTF.out.tab
# -rw-rw-r-- 1 gbonora noblelab 7.7M Nov 18 01:23 sjdbInfo.txt
# -rw-rw-r-- 1 gbonora noblelab 6.8M Nov 18 01:23 sjdbList.out.tab
# -rw-rw-r-- 1 gbonora noblelab 2.6G Nov 18 01:26 Genome
# -rw-rw-r-- 1 gbonora noblelab  21G Nov 18 01:29 SA
# -rw-rw-r-- 1 gbonora noblelab 1.5G Nov 18 01:29 SAindex










#################################################
#################################################
#################################################
# 20191118
# Nmasked spretus bowtie2 for sciATAC pipeline 

# Check sciATAC hg19 reference as an example
ls -ltrh net/trapnell/vol1/genomes/human/hg19/bowtie2/
cat /net/trapnell/vol1/genomes/human/hg19/bowtie2/hg19.fa | head
# >chr1
# NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

cd ~/refdata/mm10.Nmasked.spretus

mkdir bowtie2-index

# Use N-masked assembly:
# See above for how this was generated:
mm10SRCDIR="$HOME/refdata/mm10.Nmasked.spretus/SPRET_EiJ_N-masked"
# cat ${mm10SRCDIR}/*.fa | sed 's/^>\(.*\)/>chr\1 \1/' | sed 's/^>chrMT/>chrM/' | grep "^>"
cat ${mm10SRCDIR}/*.fa | sed 's/^>\(.*\)/>chr\1 \1/' | sed 's/^>chrMT/>chrM/' > ~/refdata/mm10.Nmasked.spretus/bowtie2-index/GRCm38-primary-assembly.spretus_N-masked.fa

# Now build bowtie2 index
cd bowtie2-index
module load bowtie2/2.2.3
bowtie2-build GRCm38-primary-assembly.spretus_N-masked.fa mm10.spretus_Nmasked

bowtie2-inspect -n ./mm10.spretus_Nmasked




















##################################################################################################
##################################################################################################
##################################################################################################
# 20191108 
# N-masked Cast/129 assemblies


#################################################
#################################################
#################################################
# 20191108
# N-masked mm10 assemblies using Cast and 129 SNPs (not validated)

mkdir ~/refdata/mm10.Nmasked.Cast129
cd ~/refdata/mm10.Nmasked.Cast129

~/bin/SNPsplit_v0.3.2/SNPsplit_genome_preparation --help | more
# Dual strain mode:
#    1) The VCF file is read and filtered for high-confidence SNPs in the strain specified with --strain <name>
#    2) The reference genome (given with --reference_genome <genome>) is read into memory, and the filtered high-
#       confidence SNP positions are incorporated as full sequence and optionally as N-masking
#    3) The VCF file is read one more time and filtered for high-confidence SNPs in strain 2 specified with --strain2 <name>
#    4) The filtered high-confidence SNP positions of strain 2 are incorporated as full sequence and optionally as N-masking
#    5) The SNP information of strain and strain 2 relative to the reference genome build are compared, and a new Ref/SNP
#       annotation is constructed whereby the new Ref/SNP information will be Strain/Strain2 (and no longer the standard
#       reference genome strain Black6 (C57BL/6J))
#    6) The full genome sequence given with --strain <name> is read into memory, and the high-confidence SNP positions between
#       Strain and Strain2 are incorporated as full sequence and optionally as N-masking
#
# ...
# 
# --dual_hybrid                 Optional. The resulting genome will no longer relate to the original reference specified with '--reference_genome'.
#                               Instead the new Reference (Ref) is defined by '--strain <strain_name>' and the new SNP genome
#                               is defined by '--strain2 <strain_name>'. '--dual_hybrid' automatically sets '--full_sequence'.
# 
#                               This will invoke a multi-step process:
#                                  1) Read/filter SNPs for first strain (specified with '--strain <strain_name>')
#                                  2) Write full SNP incorporated (and optionally N-masked) genome sequence for first strain
#                                  3) Read/filter SNPs for second strain (specified with '--strain2 <strain_name>')
#                                  4) Write full SNP incorporated (and optionally N-masked) genome sequence for second strain
#                                  5) Generate new Ref/Alt SNP annotations for Strain1/Strain2
#                                  6) Set first strain as new reference genome and construct full SNP incorporated (and optionally
#                                     N-masked) genome sequences for Strain1/Strain2

zcat $HOME/refdata/mm10pseudoSpretus/mgp.v5.merged.snps_all.dbSNP142.vcf.gz | more
# Note: uses Ensemble notation i.e. no 'chr' in front of chromosome #.
# 129S1_SvImJ
# CAST_EiJ

# So use Ensemble assembly to generate N-masked assembly then add 'chr' to chrom names.
# See above for explanation.
refFile=/net/noble/vol2/home/gbonora/refdata/iGenomes/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta
altSnpFile=$HOME/refdata/mm10pseudoSpretus/mgp.v5.merged.snps_all.dbSNP142.vcf.gz


jobfile=SNPsplit_genome_preparation.Cast129.job
echo ${jobfile}

cat << EOF > "${jobfile}"
#!/bin/bash -x
#\$ -l mfree=64G
#\$ -cwd
#\$ -j y

source $HOME/.bashrc
source $HOME/.bash_profile

hostname
printf "\n\nstart: %s\n\n" "\$( date )"

# NOTE: --dual_hybrid mode 
#       and no --skip_filtering
~/bin/SNPsplit_v0.3.2/SNPsplit_genome_preparation \
--dual_hybrid \
--strain CAST_EiJ \
--strain2 129S1_SvImJ \
--nmasking \
--vcf_file ${altSnpFile} \
--reference_genome ${refFile} \
--genome_build mm10.Nmasked.Cast129 

printf "\n\nend: %s\n\n" "\$( date )"
EOF

chmod u+x *.job

qsub ${jobfile}

ls -ltrh ~/refdata/mm10.Nmasked.Cast129
ls -ltrh ls -ltrh ~/refdata/mm10.Nmasked.Cast129/CAST_EiJ_129S1_SvImJ_dual_hybrid.based_on_mm10.Nmasked.Cast129_N-masked
# total 2.6G
# -rw-rw-r-- 1 gbonora noblelab 118M Nov  8 18:23 chr11.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 141M Nov  8 18:23 chr7.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  89M Nov  8 18:23 chrY.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 176M Nov  8 18:23 chr2.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  92M Nov  8 18:23 chr17.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 189M Nov  8 18:23 chr1.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  88M Nov  8 18:23 chr18.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 116M Nov  8 18:23 chr13.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  95M Nov  8 18:23 chr16.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 145M Nov  8 18:23 chr6.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 155M Nov  8 18:24 chr3.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 165M Nov  8 18:24 chrX.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  17K Nov  8 18:24 chrMT.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 121M Nov  8 18:24 chr9.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 116M Nov  8 18:24 chr12.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 121M Nov  8 18:24 chr14.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 101M Nov  8 18:24 chr15.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 125M Nov  8 18:24 chr8.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 151M Nov  8 18:24 chr4.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 126M Nov  8 18:24 chr10.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  60M Nov  8 18:24 chr19.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 147M Nov  8 18:24 chr5.N-masked.fa










#################################################
#################################################
#################################################
# 20191117
# Finally checked the results. They looked ok:

more CAST_EiJ_129S1_SvImJ_dual_hybrid.genome_preparation_report.txt
# Looked at positions from new Reference strain [CAST_EiJ]:               20668547
# Compared positions from new SNP strain [129S1_SvImJ]:           5173413
# ======================================================
# SNPs were the same in Ref and SNP genome (not written out):     2425219
# SNPs were present in both Ref and SNP genome but had a different sequence:      27254
# SNPs were low confidence in one strain and thus ignored:        400802
# SNPs were unique to Ref [CAST_EiJ]:                             18027951
# SNPs were unique to SNP [129S1_SvImJ]:                          2508261
# 
# Changing the genomic reference sequence to the full sequence of strain CAST_EiJ
# 
# Reading and storing all new SNPs with Ref/SNP: CAST_EiJ/129S1_SvImJ from 'all_129S1_SvImJ_SNPs_CAST_EiJ_reference.based_on_mm10.
# Nmasked.Cast129.txt'
# 1010375 positions on chromosome 11 were changed to 'N'
# 1010375 reference positions on chromosome 11 were changed to the SNP alternative base
# 
# 1135083 positions on chromosome 7 were changed to 'N'
# 1135083 reference positions on chromosome 7 were changed to the SNP alternative base
# 
# 215 positions on chromosome Y were changed to 'N'
# 215 reference positions on chromosome Y were changed to the SNP alternative base
# 
# 1292333 positions on chromosome 2 were changed to 'N'
# 1292333 reference positions on chromosome 2 were changed to the SNP alternative base
# 
# ...
# Summary
# 20563466 Ns were newly introduced into the N-masked genome for strainstrain 2 [CAST_EiJ129S1_SvImJ] in total
# 20563466 SNPs were newly introduced into the full sequence genome version for strainstrain 2 [CAST_EiJ/129S1_SvImJ] in total


# NOTE: SNPs are stored in a chr-by-chr files:
# *** MUST ENSURE THIS IS THE CASE FOR THE VALIDATE SPRETUS SNPS!!! ***

ls -tlrh SNPs_CAST_EiJ/
# total 1.8G
# -rw-rw-r-- 1 gbonora noblelab  90M Nov  8 17:22 chr11.txt
# ..
-rw-rw-r-- 1 gbonora noblelab 112M Nov  8 17:22 chr5.txt

ls -tlrh SNPs_129S1_SvImJ
# -rw-rw-r-- 1 gbonora noblelab  18M Nov  8 18:14 chr11.txt
# ..
# -rw-rw-r-- 1 gbonora noblelab  23M Nov  8 18:14 chr5.txt

cat SNPs_129S1_SvImJ/chr11.txt |head
# >11
# 46463467        11      3100685 1       T/G     1/1:60:34:0.0882353:234,60,0,.,.,.:203,48,0,.,.,.:2:34:29:2,3,0,29:17:-0.693079:.:1
# 46463470        11      3100699 1       G/C     1/1:85:25:0.04:241,85,0:214,75,0:2:37:25:0,0,0,25:0:-0.692914:.:1
# 46463471        11      3100701 1       G/C     1/1:79:23:0:244,79,0,.,.,.:217,69,0,.,.,.:2:39:23:0,0,0,23:0:-0.692717:.:1
# 46463504        11      3100850 1       C/A     1/1:5:14:0:145,5,1:129,0,2:2:29:12:2,0,12,0:0:-0.680642:.:1
# 46463715        11      3102606 1       G/T     1/1:127:75:0:292,241,0:255,226,0:2:46:75:0,0,28,47:0:-0.693147:.:1
# 46464748        11      3114433 1       C/A     1/1:26:26:0.0384615:277,26,0:255,18,0:2:48:20:6,0,17,3:0:-0.692067:.:1










#################################################
#################################################
#################################################
# 20191108
# Cat new hg19 + Nmasked Cast/129 mm10 combo assembly 
# for scihic pipeline and bowtie2 index

# NOTE: must add 'chr' to mouse chromosome names

mkdir $HOME/refdata/hg19mm10.Nmasked.Cast129
cd $HOME/refdata/hg19mm10.Nmasked.Cast129

hg19SRCDIR="$HOME/refdata/hg19mm10"
mm10SRCDIR="$HOME/refdata/mm10.Nmasked.Cast129/CAST_EiJ_129S1_SvImJ_dual_hybrid.based_on_mm10.Nmasked.Cast129_N-masked"

# Combo fasta
# Must add 'chr' to chrom names
# cat ${mm10SRCDIR}/*.fa | sed 's/^>/>chr/' | sed 's/^>chrMT/>chrM/' | grep "^>"
cat ${hg19SRCDIR}/hg19.anno.fa <( cat ${mm10SRCDIR}/*.fa | sed 's/^>/>chr/' | sed 's/^>chrMT/>chrM/' ) > hg19.mm10.Nmasked.Cast129.combo.fa
cat hg19.mm10.Nmasked.Cast129.combo.fa | grep "^>"

# Now build bowtie2 index
module load bowtie2/2.2.3
bowtie2-build hg19.mm10.Nmasked.Cast129.combo.fa mm10hg19combo

mkdir ./bowtie2-index
mv mm10hg19combo* bowtie2-index

bowtie2-inspect -n ./bowtie2-index/mm10hg19combo
# ...
# human_chr18_gl000207_random
# chr10
# chr11
# chr12
# chr13
# chr14
# chr15
# chr16
# chr17
# chr18
# chr19
# chr1
# chr2
# chr3
# chr4
# chr5
# chr6
# chr7
# chr8
# chr9
# chrM
# chrX
# chrY










#################################################
#################################################
#################################################
# 20191108
# Nmasked Cast129 assembly for use with sciRNA-seq STAR alignment 

# Look at Jonathan Packer's assembly:
# STAR_INDEX=/net/trapnell/vol1/jspacker/STAR/mouse/GRCm38-primary-assembly

ls -ltrh /net/trapnell/vol1/jspacker/STAR/

# Already copied locally:
ls -ltrh  ~/refdata/STAR

more  ~/refdata/STAR/STAR-genomeGenerate.sh
# #$ -S /bin/bash
# #$ -l mfree=8G
# #$ -l h_rt=3:0:0:0
# #$ -pe serial 12
# 
# module load STAR/2.5.2b
# 
# FASTA=$1
# GTF=$2
# OUTPUT=$3
# WORKING_DIR=$4
# 
# cd $WORKING_DIR
# 
# STAR \
#     --runThreadN 12 \
#     --runMode genomeGenerate \
#     --genomeDir $OUTPUT \
#     --genomeFastaFiles $FASTA \
#     --sjdbGTFfile $GTF \
#     --outTmpDir $WORKING_DIR/tmp

	
ls -ltrh  ~/refdata/STAR/mouse
zcat ~/refdata/STAR/mouse/GRCm38.primary_assembly.genome.fa.gz | grep -e '^>'
# >chr1 1
# >chr2 2
# >chr3 3
# ..
# >chrX X
# >chrY Y
# >chrM MT
# >GL456210.1 GL456210.1
# ...


# Replace Jonathan's assembly with N-masked one.
# See above for how this was generated:
mm10SRCDIR="$HOME/refdata/mm10.Nmasked.Cast129/CAST_EiJ_129S1_SvImJ_dual_hybrid.based_on_mm10.Nmasked.Cast129_N-masked"
cat ${mm10SRCDIR}/*.fa | sed 's/^>\(.*\)/>chr\1 \1/' | sed 's/^>chrMT/>chrM/' | grep "^>"
cat ${mm10SRCDIR}/*.fa | sed 's/^>\(.*\)/>chr\1 \1/' | sed 's/^>chrMT/>chrM/' > ~/refdata/STAR/mouse/GRCm38-primary-assembly.Cast129_N-masked.fa

mkdir ~/refdata/STAR/mouse/GRCm38-primary-assembly.Cast129_N-masked

# cp ~/refdata/STAR/STAR-genomeGenerate.sh  ~/refdata/STAR/STAR-genomeGenerate.GB.sh 
vi ~/refdata/STAR/STAR-genomeGenerate.GB.sh 

cd ~/refdata/STAR
rm -rf mouse/tmp
qsub -j y -cwd ./STAR-genomeGenerate.GB.sh \
~/refdata/STAR/mouse/GRCm38-primary-assembly.Cast129_N-masked.fa \
~/refdata/STAR/mouse/gencode.vM15.primary_assembly.annotation.gtf \
~/refdata/STAR/mouse/GRCm38-primary-assembly.Cast129_N-masked \
~/refdata/STAR/mouse

mv mouse/Log.out mouse/Log.Cast129.out	

ls -ltrh ~/refdata/STAR/mouse/GRCm38-primary-assembly.Cast129_N-masked
# total 25G
# -rw-rw-r-- 1 gbonora noblelab  926 Nov 10 20:00 genomeParameters.txt
# -rw-rw-r-- 1 gbonora noblelab  120 Nov 10 20:02 chrName.txt
# -rw-rw-r-- 1 gbonora noblelab  211 Nov 10 20:02 chrLength.txt
# -rw-rw-r-- 1 gbonora noblelab  236 Nov 10 20:02 chrStart.txt
# -rw-rw-r-- 1 gbonora noblelab  331 Nov 10 20:02 chrNameLength.txt
# -rw-rw-r-- 1 gbonora noblelab  27M Nov 10 20:48 exonGeTrInfo.tab
# -rw-rw-r-- 1 gbonora noblelab 1.1M Nov 10 20:48 geneInfo.tab
# -rw-rw-r-- 1 gbonora noblelab 8.0M Nov 10 20:48 transcriptInfo.tab
# -rw-rw-r-- 1 gbonora noblelab  11M Nov 10 20:48 exonInfo.tab
# -rw-rw-r-- 1 gbonora noblelab 6.8M Nov 10 20:48 sjdbList.fromGTF.out.tab
# -rw-rw-r-- 1 gbonora noblelab 7.7M Nov 10 20:48 sjdbInfo.txt
# -rw-rw-r-- 1 gbonora noblelab 6.8M Nov 10 20:48 sjdbList.out.tab
# -rw-rw-r-- 1 gbonora noblelab 2.6G Nov 10 20:52 Genome
# -rw-rw-r-- 1 gbonora noblelab  21G Nov 10 20:55 SA
# -rw-rw-r-- 1 gbonora noblelab 1.5G Nov 10 20:55 SAindex










#################################################
#################################################
#################################################
# 20191108
# Nmasked Cast129 bowtie2 for sciATAC pipeline 

# Check sciATAC hg19 reference as an example
ls -ltrh net/trapnell/vol1/genomes/human/hg19/bowtie2/
cat /net/trapnell/vol1/genomes/human/hg19/bowtie2/hg19.fa | head
# >chr1
# NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

cd ~/refdata/mm10.Nmasked.Cast129

mkdir bowtie2-index

# Use N-masked assembly:
# See above for how this was generated:
mm10SRCDIR="$HOME/refdata/mm10.Nmasked.Cast129/CAST_EiJ_129S1_SvImJ_dual_hybrid.based_on_mm10.Nmasked.Cast129_N-masked"
# cat ${mm10SRCDIR}/*.fa | sed 's/^>\(.*\)/>chr\1 \1/' | sed 's/^>chrMT/>chrM/' | grep "^>"
cat ${mm10SRCDIR}/*.fa | sed 's/^>\(.*\)/>chr\1 \1/' | sed 's/^>chrMT/>chrM/' > ~/refdata/mm10.Nmasked.Cast129/bowtie2-index/GRCm38-primary-assembly.Cast129_N-masked.fa

# Now build bowtie2 index
cd bowtie2-index
module load bowtie2/2.2.3
bowtie2-build GRCm38-primary-assembly.Cast129_N-masked.fa mm10.Cast129_Nmasked

bowtie2-inspect -n ./mm10.Cast129_Nmasked




















##################################################################################################
##################################################################################################
##################################################################################################
# Get/Send refdata code & data

######
# DOWN

# code
rsync -auvPKLn --include='*.sh' --include='*.py' --include='*.R' --include='*.*awk' --include='*.c' --include='*/' --prune-empty-dirs  --exclude='*' gbonora@nexus.gs.washington.edu:NobleLab/refdata/ ~/NobleLab/refdata/
# data
rsync -auvPKLn --max-size=1M --include='*.bed' --include='*.txt'  --include='*.tsv' --include='*.html'  --include='*/' --prune-empty-dirs  --exclude='*' gbonora@nexus.gs.washington.edu:NobleLab/refdata/ ~/NobleLab/refdata/



######
# UP

# code
rsync -auvPKLn --include='*.sh' --include='*.py' --include='*.R' --include='*.*awk' --include='*.c' --include='*/' --prune-empty-dirs  --exclude='*' ~/NobleLab/refdata/ gs:NobleLab/refdata/
# data
rsync -auvPKLn --max-size=20M --include='*.xlsx' --include='*.rds' --include='*.sqlite' --include='*.pdf' --include='*.bed' --include='*.txt'  --include='*.tsv' --include='*.html'  --include='*/' --prune-empty-dirs  --exclude='*'  ~/NobleLab/refdata/ gs:NobleLab/refdata/
