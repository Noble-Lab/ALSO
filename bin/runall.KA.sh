#!/bin/bash

#  Remote ---------------------------------------------------------------------

#TODO Code for work with faSize
#TODO Code for renaming chromosomes and munging, etc. (chromosome_sizes.txt)
#TODO Code for setting up whitelists
#TODO Checks, e.g., don't run if file present, etc.

genomes_directory="${HOME}/genomes"


# -------------------------------------
#  Get motif database file (such as those from JASPAR) to use for motif searching in peaks
file_JASPAR="JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.pfm"
[[ -f "data/motifs/${file_JASPAR}" ]] ||
    {
        curl "http://jaspar.genereg.net/download/CORE/${file_JASPAR/.pfm/.txt}" \
        > "data/motifs/${file_JASPAR}"
    }


# -------------------------------------
#  Ensembl.129S1-SvImJ, https://www.ebi.ac.uk/ena/browser/view/GCA_001624185.1
name="Ensembl.129S1-SvImJ"
mkdir -p "${genomes_directory}/${name}"
curl "ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/lvx/LVXH01.fasta.gz" \
> "${genomes_directory}/${name}/LVXH01.fasta.gz"
curl "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_001624185.1?download=true&gzip=true" \
> "${genomes_directory}/${name}/GCA_001624185.1.fa.gz"

array_files=(
    CHECKSUMS
    Mus_musculus_129s1svimj.129S1_SvImJ_v1.104.abinitio.gtf.gz
    Mus_musculus_129s1svimj.129S1_SvImJ_v1.104.chr.gtf.gz
    Mus_musculus_129s1svimj.129S1_SvImJ_v1.104.gtf.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gtf"
    mkdir -p "${genomes_directory}/${name}/gtf"
    curl "http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus_129s1svimj/${i}" \
    > "${genomes_directory}/${name}/gtf/${i}"
done


# -------------------------------------
#  Ensembl.C57BL-6NJ, https://www.ebi.ac.uk/ena/browser/view/GCA_001632555.1
name="Ensembl.C57BL-6NJ"
mkdir -p "${genomes_directory}/${name}"
curl "ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/lvx/LVXM01.fasta.gz" \
> "${genomes_directory}/${name}/LVXM01.fasta.gz"
curl "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_001632555.1?download=true&gzip=true" \
> "${genomes_directory}/${name}/GCA_001632555.1.fa.gz"

array_files=(
    CHECKSUMS
    Mus_musculus_c57bl6nj.C57BL_6NJ_v1.104.abinitio.gtf.gz
    Mus_musculus_c57bl6nj.C57BL_6NJ_v1.104.chr.gtf.gz
    Mus_musculus_c57bl6nj.C57BL_6NJ_v1.104.gtf.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gtf"
    mkdir -p "${genomes_directory}/${name}/gtf"
    curl "http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus_c57bl6nj/${i}" \
    > "${genomes_directory}/${name}/gtf/${i}"
done


# -------------------------------------
#  Ensembl.CAST-EiJ, https://www.ebi.ac.uk/ena/browser/view/GCA_001624445.1
name="Ensembl.CAST-EiJ"
mkdir -p "${genomes_directory}/${name}"
curl "ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/lvx/LVXN01.fasta.gz" \
> "${genomes_directory}/${name}/LVXN01.fasta.gz"
curl "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_001624445.1?download=true&gzip=true" \
> "${genomes_directory}/${name}/GCA_001624445.1.fa.gz"

array_files=(
    CHECKSUMS
    Mus_musculus_casteij.CAST_EiJ_v1.104.abinitio.gtf.gz
    Mus_musculus_casteij.CAST_EiJ_v1.104.chr.gtf.gz
    Mus_musculus_casteij.CAST_EiJ_v1.104.gtf.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gtf"
    mkdir -p "${genomes_directory}/${name}/gtf"
    curl "http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus_casteij/${i}" \
    > "${genomes_directory}/${name}/gtf/${i}"
done


# -------------------------------------
#  Ensembl.SPRET-EiJ, https://www.ebi.ac.uk/ena/browser/view/GCA_001624865.1
name="Ensembl.SPRET-EiJ"
mkdir -p "${genomes_directory}/${name}"
curl "ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/lvx/LVXV01.fasta.gz" \
> "${genomes_directory}/${name}/LVXV01.fasta.gz"
curl "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_001624865.1?download=true&gzip=true" \
> "${genomes_directory}/${name}/GCA_001624865.1.fa.gz"

array_files=(
    CHECKSUMS
    Mus_spretus.SPRET_EiJ_v1.104.abinitio.gtf.gz
    Mus_spretus.SPRET_EiJ_v1.104.chr.gtf.gz
    Mus_spretus.SPRET_EiJ_v1.104.gtf.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gtf"
    mkdir -p "${genomes_directory}/${name}/gtf"
    curl "http://ftp.ensembl.org/pub/release-104/gtf/mus_spretus/${i}" \
    > "${genomes_directory}/${name}/gtf/${i}"
done


# -------------------------------------
#  Ensembl.GRCm39
name="Ensembl.GRCm39"
mkdir -p "${genomes_directory}/${name}"
curl "http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" \
> "${genomes_directory}/${name}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"

array_files=(
    CHECKSUMS
    Mus_musculus.GRCm39.104.abinitio.gtf.gz
    Mus_musculus.GRCm39.104.chr.gtf.gz
    Mus_musculus.GRCm39.104.gtf.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gtf"
    mkdir -p "${genomes_directory}/${name}/gtf"
    curl "http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/${i}" > "${genomes_directory}/${name}/gtf/${i}"
done


# -------------------------------------
#  Ensembl.GRCm38
name="Ensembl.GRCm38"
mkdir -p "${genomes_directory}/${name}"
curl "ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz" \
> "${genomes_directory}/${name}/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"

array_files=(
    CHECKSUMS
    Mus_musculus.GRCm38.102.abinitio.gtf.gz
    Mus_musculus.GRCm38.102.chr.gtf.gz
    Mus_musculus.GRCm38.102.gtf.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gtf"
    mkdir -p "${genomes_directory}/${name}/gtf"
    curl "http://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/${i}" > "${genomes_directory}/${name}/gtf/${i}"
done


# -------------------------------------
array_samples=(
    129S1-SvImJ
    C57BL-6NJ
    CAST-EiJ
    SPRET-EiJ
)

#  Get effective genome sizes for samples
for i in "${array_samples[@]}"; do
    fasta_in="$(find "${genomes_directory}/Ensembl.${i}" -name "*.fa" -printf "%P\n")"
    fasta_out="${i}.fa"

    [[ ! -e "${fasta_out}" ]] || rm "${fasta_out}"

    cat "${genomes_directory}/Ensembl.${i}/${fasta_in}" | \
    seqkit grep -vrp "LVX|KV|CM004188|CM004181|CM004277|CM004179" | \
    seqkit rmdup > "${genomes_directory}/Ensembl.${i}/${fasta_out}"
done

for i in "${array_samples[@]}"; do
    fasta_in="${i}.fa"
    count="$(faCount -summary "${genomes_directory}/Ensembl.${i}/${fasta_in}")"
    echo -e "${fasta_in}\n" "${count}\n" >> "data/effective-genome-sizes.txt"
done

#  Index .fa files
for i in "${array_samples[@]}"; do
    fasta_GCA="$(find "${genomes_directory}/Ensembl.${i}" -name "GCA*.fa" -printf "%P\n")"
    fasta_sample="$(find "${genomes_directory}/Ensembl.${i}" -name "${i}.fa" -printf "%P\n")"    
    samtools faidx "${genomes_directory}/Ensembl.${i}/${fasta_GCA}"
    samtools faidx "${genomes_directory}/Ensembl.${i}/${fasta_sample}"
done

#  Prepare TSS, gene-body definitions
for i in "${array_samples[@]}"; do
    Ensembl_directory="${genomes_directory}/Ensembl.${i}"
    gtf="$(find "${Ensembl_directory}/gtf" -name "*.chr.gtf.gz" -printf "%P\n")"
    
    temp_TSS="${Ensembl_directory}/${i}.TSS.temp.bed"
    temp_gene_body="${Ensembl_directory}/${i}.temp.gene_bodies.bed"
    
    whitelist="${Ensembl_directory}/${i}.whitelist_regions.bed"

    final_TSS="${Ensembl_directory}/${i}.TSS.bed"
    final_gene_body="${Ensembl_directory}/${i}.gene_bodies.bed"

    chromosome_sizes="${Ensembl_directory}/${i}.chromosome_sizes.txt"

    final_gene_body_2k="${Ensembl_directory}/${i}.gene_bodies.plus_2kb_upstream.bed"

    [[ ! -e "${temp_TSS}" ]] || rm "${temp_TSS}"
    [[ ! -e "${temp_gene_body}" ]] || rm "${temp_gene_body}"
    [[ ! -e "${final_TSS}" ]] || rm "${final_TSS}"
    [[ ! -e "${final_gene_body}" ]] || rm "${final_gene_body}"
    [[ ! -e "${final_gene_body_2k}" ]] || rm "${final_gene_body_2k}"
    
    Rscript bin/generate_TSS_file.R \
    --infile "${Ensembl_directory}/gtf/${gtf}" \
    --outfile "${temp_TSS}"

    Rscript bin/generate_gene_body_file.R \
    --infile "${Ensembl_directory}/gtf/${gtf}" \
    --outfile "${temp_gene_body}"

    cat "${temp_TSS}" \
    | grep -i -E "protein|lincRNA|antisense|TR_V|TR_D|TR_J|TR_C|IG_V|IG_D|IG_J|IG_C|IG_LV" \
    | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$7,1,$6}' \
    | sort -k1,1V -k2,2n -k3,3n \
    | bedtools intersect -a stdin -b "${whitelist}" \
    | uniq > "${final_TSS}"

    cat "${temp_gene_body}" \
    | grep -i -E "protein|lincRNA|antisense|TR_V|TR_D|TR_J|TR_C|IG_V|IG_D|IG_J|IG_C|IG_LV" \
    | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$7,1,$6}' \
    | sort -k1,1V -k2,2n -k3,3n \
    | bedtools intersect -a stdin -b "${whitelist}" \
    | uniq > "${final_gene_body}"

    bedtools slop -i "${final_gene_body}" -s -l 2000 -r 0 -g "${chromosome_sizes}" \
    | sort -k1,1V -k2,2n -k3,3n > "${final_gene_body_2k}"
done


#  Local ----------------------------------------------------------------------
local="Dropbox/UW/projects-etc"
work_directory="${HOME}/${local}/2021_kga0_4dn-mouse-cross/data/genomes"

cd "${work_directory}" || exit 1

ftp="http://ftp.ensembl.org/pub/release-104/gtf"
wget "${ftp}/mus_musculus_129s1svimj/Mus_musculus_129s1svimj.129S1_SvImJ_v1.104.chr.gtf.gz"
wget "${ftp}//mus_musculus_c57bl6nj/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.104.chr.gtf.gz"
wget "${ftp}//mus_musculus_casteij/Mus_musculus_casteij.CAST_EiJ_v1.104.chr.gtf.gz"
wget "${ftp}//mus_spretus/Mus_spretus.SPRET_EiJ_v1.104.chr.gtf.gz"

cd "${HOME}/${local}/2021_kga0_4dn-mouse-cross" || exit 1

Rscript bin/generate_TSS_file.R \
--infile "${work_directory}/Mus_musculus_129s1svimj.129S1_SvImJ_v1.104.chr.gtf.gz" \
--outfile "${work_directory}/129S1-SvImJ.TSS.temp.bed"

Rscript bin/generate_TSS_file.R \
--infile "${work_directory}/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.104.chr.gtf.gz" \
--outfile "${work_directory}/C57BL-6NJ.TSS.temp.bed"

Rscript bin/generate_TSS_file.R \
--infile "${work_directory}/Mus_musculus_casteij.CAST_EiJ_v1.104.chr.gtf.gz" \
--outfile "${work_directory}/CAST-EiJ.TSS.temp.bed"

Rscript bin/generate_TSS_file.R \
--infile "${work_directory}/Mus_spretus.SPRET_EiJ_v1.104.chr.gtf.gz" \
--outfile "${work_directory}/SPRET-EiJ.TSS.temp.bed"

Rscript bin/generate_gene_body_file.R \
--infile "${work_directory}/Mus_musculus_129s1svimj.129S1_SvImJ_v1.104.chr.gtf.gz" \
--outfile "${work_directory}/129S1-SvImJ.gene_bodies.temp.bed"

Rscript bin/generate_gene_body_file.R \
--infile "${work_directory}/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.104.chr.gtf.gz" \
--outfile "${work_directory}/C57BL-6NJ.gene_bodies.temp.bed"

Rscript bin/generate_gene_body_file.R \
--infile "${work_directory}/Mus_musculus_casteij.CAST_EiJ_v1.104.chr.gtf.gz" \
--outfile "${work_directory}/CAST-EiJ.gene_bodies.temp.bed"

Rscript bin/generate_gene_body_file.R \
--infile "${work_directory}/Mus_spretus.SPRET_EiJ_v1.104.chr.gtf.gz" \
--outfile "${work_directory}/SPRET-EiJ.gene_bodies.temp.bed"

cd "${work_directory}" || exit 1

for i in 129S1-SvImJ C57BL-6NJ CAST-EiJ SPRET-EiJ; do
    cat "${i}.TSS.temp.bed" \
    | grep -i -E "protein|lincRNA|antisense|TR_V|TR_D|TR_J|TR_C|IG_V|IG_D|IG_J|IG_C|IG_LV" \
    | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$7,1,$6}' \
    | sort -k1,1V -k2,2n -k3,3n \
    | bedtools intersect -a stdin -b "${i}.whitelist_regions.bed" \
    | uniq > "${i}.TSS.bed"

    cat "${i}.gene_bodies.temp.bed" \
    | grep -i -E "protein|lincRNA|antisense|TR_V|TR_D|TR_J|TR_C|IG_V|IG_D|IG_J|IG_C|IG_LV" \
    | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$7,1,$6}' \
    | sort -k1,1V -k2,2n -k3,3n \
    | bedtools intersect -a stdin -b "${i}.whitelist_regions.bed" \
    | uniq > "${i}.gene_bodies.bed"

    bedtools slop -i "${i}.gene_bodies.bed" -s -l 2000 -r 0 -g "${i}.chromosome_sizes.txt" \
    | sort -k1,1V -k2,2n -k3,3n > "${i}.gene_bodies.plus_2kb_upstream.bed"
done

gzip ./*.bed

