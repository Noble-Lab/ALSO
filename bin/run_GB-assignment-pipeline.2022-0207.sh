#!/bin/bash

path_base="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
path_tmp="/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora"
cd "${path_tmp}" ||
{
    echo "Directory not found. Exiting."
    return 1 2> /dev/null
    exit 1
}

cp \
"processed.mate-paired.129S1.chrX.bam" \
"processed.mate-paired.CAST.chrX.bam" \
"processed.mate-paired.mm10.chrX.bam" \
"${path_base}/data/kga0"

cd "${path_base}" ||
{
    echo "Directory not found. Exiting."
    return 1 2> /dev/null
    exit 1
}

path_script="./bin/gbonora"
path_data="./data"
bam="kga0/processed.mate-paired.mm10.chrX.bam"
file_ref_SNPs="gbonora/129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.NoCast.vcf.gz"  # ref = 129
file_alt_SNPs="gbonora/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.No129.vcf.gz"  # alt = CAST
thresh_SNP="1"
thresh_Q="30"
path_results="./results/kga0/$(date '+%Y-%m%d')_segregated_reads.thresh_SNP_${thresh_SNP}.thresh_Q_${thresh_Q}"

[[ -d "${path_results}" ]] || mkdir -p "${path_results}"

echo ""
echo "${path_script}"
echo "${path_data}"
echo "${path_data}/${bam}"
echo "${path_data}/${file_ref_SNPs}"
echo "${path_data}/${file_alt_SNPs}"
echo "${path_results}"
echo "${thresh_SNP}"
echo "${thresh_Q}"

echo ""
head "${path_script}/segregate-sciATAC.PE.CIGARaware.refAndAltSNPs.SNPthresh.Qthresh.py"
echo ""
ls "${path_data}/${bam}"
echo ""
samtools view -h "${path_data}/${bam}" | head -50
echo ""
gzip -cd "${path_data}/${file_ref_SNPs}" | head
echo ""
gzip -cd "${path_data}/${file_alt_SNPs}" | head
echo ""

python \
"${path_script}/segregate-sciATAC.PE.CIGARaware.refAndAltSNPs.SNPthresh.Qthresh.py" \
"${path_data}/${bam}" \
"${path_data}/${file_ref_SNPs}" \
"${path_data}/${file_alt_SNPs}" \
"${path_results}" \
"${thresh_SNP}" \
"${thresh_Q}"

mv "${path_results}."*".bam" "${path_results}/"
