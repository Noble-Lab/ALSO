#!/bin/bash

#  Calling GB's allelic-segregation script

path_script="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/doc/gbonora"

bam="./dedup.mm10.chrX.bam"
# bam="./mm10-CAST-129-Nmasked.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.bam"
# bam="mm10-CAST-129-Nmasked.F121-6-CASTx129.undifferentiated.dedup.bam"
ref_SNP_file="./mm10CAST129SNPs/129S1_SvImJ.mgp.v5.snps.dbSNP142.PASS.NoCast.vcf.gz"
alt_SNP_file="./mm10CAST129SNPs/CAST_EiJ.mgp.v5.snps.dbSNP142.PASS.No129.vcf.gz"
thresh_SNP="1"
# thresh_Q="30"
thresh_Q="10"
sub_outdir="./segregated_reads.thresh_SNP_${thresh_SNP}.thresh_Q_${thresh_Q}"
# sub_outdir="segregatedReads.SNPTHRESH${thresh_SNP}.Q${thresh_Q}"

python3 "${path_script}/segregate-sciATAC.PE.CIGARaware.refAndAltSNPs.SNPthresh.Qthresh.py" \
"${bam}" \
"${ref_SNP_file}" \
"${alt_SNP_file}" \
"${sub_outdir}" \
"${thresh_SNP}" \
"${thresh_Q}"
