#!/bin/bash

#  submit_align_Berletch-data_individual_outSAMunmapped-None.sh
#  KA

d_base="/net/noble/vol1/home/kga0/genomes"
d_data="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/data/Berletch_unpublished/fastq"
gtf_SPRET="${d_base}/Ensembl.SPRET-EiJ/gtf/106/Mus_spretus.SPRET_EiJ_v1.106.gtf"
gtf_GRCm39="${d_base}/Ensembl.GRCm39/gtf/106/Mus_musculus.GRCm39.106.gtf"

for i in "${gtf_SPRET}" "${gtf_GRCm39}"; do
    if [[ ${i} == "${gtf_SPRET}" ]]; then
        index="/net/noble/vol1/home/kga0/genomes/Ensembl.SPRET-EiJ/STAR_no-GTF"
        align_to="SPRET-EiJ"
    elif [[ ${i} == "${gtf_GRCm39}" ]]; then
        index="/net/noble/vol1/home/kga0/genomes/Ensembl.GRCm39/STAR_no-GTF"
        align_to="mm11"
    fi

    #  heart run 1
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/181112_NS500488_0731_AHTTCCBGX7_demultiplexed/heart_S1_R1_001.AHTTCCBGX7.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "75"

    #  heart run 2
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/181212_NS500488_0750_AHG5WVBGX9_demultiplexed/heart_S1_R1_001.AHG5WVBGX9.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "75"

    #  heart run 3
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/190205_NS500488_0769_AHCG2FBGX9_demultiplexed/heart_S1_R1_001.AHCG2FBGX9.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "75"


    #  kidney run 1
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/181112_NS500488_0731_AHTTCCBGX7_demultiplexed/kidney_S2_R1_001.AHTTCCBGX7.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "75"

    #  kidney run 2
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/181212_NS500488_0750_AHG5WVBGX9_demultiplexed/kidney_S2_R1_001.AHG5WVBGX9.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "75"

    #  kidney run 3
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/190205_NS500488_0769_AHCG2FBGX9_demultiplexed/kidney_S2_R1_001.AHCG2FBGX9.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "75"
done
