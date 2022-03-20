# Each organism will have its own directory with same set of files

########################################
# Any files common to multiple organisms
########################################
mkdir -p common_files
mkdir -p common_files/motifs
if [ ! -f common_files/motifs/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.pfm ]; then
    curl http://jaspar.genereg.net/download/CORE/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt > common_files/motifs/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.pfm
fi

########################################
# HG19
########################################
echo 'Prepping hg19...'
mkdir -p hg19
cd hg19

# BOWTIE REFERENCE
echo '   bowtie...'
ln -s /net/shendure/vol10/projects/scATAC/nobackup/genomes/hg19/bowtie bowtie

# WHITELIST REGIONS (All standard chromosomes from end to end)
if [ ! -f whitelist_regions.bed ]; then
    echo '   whitelist regions...'
    cat ~ajh24/common_data/ucsc/goldenPath/hg19/bigZips/hg19.chrom.sizes \
    | grep -v chrM \
    | grep -v random \
    > chromosome_sizes.txt

    cat chromosome_sizes.txt \
    | awk 'BEGIN{OFS="\t"}{ print $1,"1",$2;}' > whitelist_regions.bed


fi

# TSS AND GENE BODY + 2KB DEFINITIONS
if [ ! -f tss.bed.gz ]; then
    echo '   gtf...'
    mkdir -p gtf
    curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh37_mapping/gencode.v29lift37.annotation.gtf.gz > gtf/gencode.v29lift37.annotation.gtf.gz

    Rscript ~ajh24/common_data/sciatac_pipeline_data/generate_tss_file.R gtf/gencode.v29lift37.annotation.gtf.gz tss.temp.bed
    Rscript ~ajh24/common_data/sciatac_pipeline_data/generate_gene_body_file.R gtf/gencode.v29lift37.annotation.gtf.gz gene_bodies.temp.bed

    cat tss.temp.bed \
    | grep -i -E "protein|lincRNA|antisense|TR_V|TR_D|TR_J|TR_C|IG_V|IG_D|IG_J|IG_C|IG_LV" \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,1,$6}' \
    | sort -k1,1V -k2,2n -k3,3n \
    | bedtools intersect -a stdin -b whitelist_regions.bed \
    | uniq \
    | gzip > tss.bed.gz

    cat gene_bodies.temp.bed \
    | grep -i -E "protein|lincRNA|antisense|TR_V|TR_D|TR_J|TR_C|IG_V|IG_D|IG_J|IG_C|IG_LV" \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,1,$6}' \
    | sort -k1,1V -k2,2n -k3,3n \
    | bedtools intersect -a stdin -b whitelist_regions.bed \
    | uniq \
    | gzip > gene_bodies.bed.gz

    bedtools slop -i gene_bodies.bed.gz -s -l 2000 -r 0 -g chromosome_sizes.txt \
    | sort -k1,1V -k2,2n -k3,3n \
    | gzip > gene_bodies.plus_2kb_upstream.bed.gz

    rm tss.temp.bed
    rm gene_bodies.temp.bed
fi

# REFERENCE (FASTA)
if [ ! -f reference/hg19.fa ]; then
    echo '   reference...'
    mkdir -p reference
    zcat ~ajh24/common_data/ucsc/goldenPath/hg19/bigZips/hg19.fa.gz > reference/hg19.fa
    # Also index the fasta
    python -c "import pyfasta; test = pyfasta.Fasta('reference/hg19.fa')"
fi

cd ..

echo 'Done.'

########################################
# MM9
########################################
echo 'Prepping mm9...'
mkdir -p mm9
cd mm9

# BOWTIE REFERENCE
echo '   bowtie...'
ln -s /net/shendure/vol10/projects/scATAC/nobackup/genomes/mm9/tophat/ bowtie

# WHITELIST REGIONS (All standard chromosomes from end to end)
if [ ! -f whitelist_regions.bed ]; then
    echo '   whitelist regions...'
    cat ~ajh24/common_data/ucsc/goldenPath/mm9/bigZips/mm9.chrom.sizes \
    | grep -v chrM \
    | grep -v random \
    > chromosome_sizes.txt

    cat chromosome_sizes.txt \
    |  awk 'BEGIN{OFS="\t"}{ print $1,"1",$2;}' > whitelist_regions.bed
fi

# GTF FILE + TSS DEFINITIONS
if [ ! -f tss.bed.gz ]; then
    echo '   gtf...'
    mkdir -p gtf
    curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz > gtf/gencode.vM1.annotation.gtf.gz

    Rscript ~ajh24/common_data/sciatac_pipeline_data/generate_tss_file.R gtf/gencode.vM1.annotation.gtf.gz tss.temp.bed
    Rscript ~ajh24/common_data/sciatac_pipeline_data/generate_gene_body_file.R gtf/gencode.vM1.annotation.gtf.gz gene_bodies.temp.bed

    cat tss.temp.bed \
    | grep -i -E "protein|lincRNA|antisense|TR_V|TR_D|TR_J|TR_C|IG_V|IG_D|IG_J|IG_C|IG_LV" \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,1,$6}' \
    | sort -k1,1V -k2,2n -k3,3n \
    | bedtools intersect -a stdin -b whitelist_regions.bed \
    | uniq \
    | gzip > tss.bed.gz

    cat gene_bodies.temp.bed \
    | grep -i -E "protein|lincRNA|antisense|TR_V|TR_D|TR_J|TR_C|IG_V|IG_D|IG_J|IG_C|IG_LV" \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,1,$6}' \
    | sort -k1,1V -k2,2n -k3,3n \
    | bedtools intersect -a stdin -b whitelist_regions.bed \
    | uniq \
    | gzip > gene_bodies.bed.gz

    bedtools slop -i gene_bodies.bed.gz -s -l 2000 -r 0 -g chromosome_sizes.txt \
    | sort -k1,1V -k2,2n -k3,3n \
    | gzip > gene_bodies.plus_2kb_upstream.bed.gz

    rm tss.temp.bed
    rm gene_bodies.temp.bed
fi

# REFERENCE (FASTA)
if [ ! -f reference/mm9.fa ]; then
    echo '   reference...'
    mkdir -p reference
    zcat ~ajh24/common_data/ucsc/goldenPath/mm9/chromosomes/mm9.fasta.gz > reference/mm9.fa

    ## Also index the reference
    python -c "import pyfasta; test = pyfasta.Fasta('reference/mm9.fa')"
fi


########################################
# BARNYARD HUMAN/MOUSE
########################################
mkdir -p hg19_mm9

cat <(awk '{ print "hg19" $0 }' hg19/chromosome_sizes.txt) <(awk '{ print "mm9" $0 }' mm9/chromosome_sizes.txt) \
| sort -k1,1V -k2,2n -k3,3n \
> hg19_mm9/chromosome_sizes.txt

cat hg19_mm9/chromosome_sizes.txt \
| awk 'BEGIN{OFS="\t"}{ print $1,"1",$2;}' > hg19_mm9/whitelist_regions.bed

cat <(zcat hg19/tss.bed.gz | awk '{ print "hg19" $0 }') <(zcat mm9/tss.bed.gz | awk '{ print "mm9" $0 }') \
| sort -k1,1V -k2,2n -k3,3n \
| gzip > hg19_mm9/tss.bed.gz


# Bowtie indices
ln -s /net/gs/vol1/home/tdurham/proj/hnd-1_rnaseq/ref_data/WS230 bowtie_indices/WS230
ln -s /net/shendure/vol10/projects/scATAC/nobackup/genomes/dm3 bowtie_indices/dm3
ln -s /net/shendure/vol10/projects/scATAC/nobackup/genomes/mm9/tophat/ bowtie_indices/mm9
ln -s /net/shendure/vol10/projects/scATAC/nobackup/genomes/hg19/bowtie/ bowtie_indices/hg19
ln -s /net/shendure/vol10/projects/scATAC/nobackup/genomes/hg19_mm9/tophat/ bowtie_indices/hg19_mm9
