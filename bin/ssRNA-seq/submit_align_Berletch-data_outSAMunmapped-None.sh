#!/bin/bash

#  submit_align_Berletch-data_outSAMunmapped-None.sh
#  KA
#  2022-0721

module load STAR/2.7.6a

d_base="/net/noble/vol1/home/kga0/genomes"
d_data="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/data/Berletch_PRJNA256188/fastq"
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

    #  brain, rep 1
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/brain_rep_1.SRR1525404-05-06.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "35"

    sleep 0.5

    #  brain, rep 1, tech 1
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/brain_rep_1.SRR1525404.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "35"

    sleep 0.5

    #  brain, rep 1, tech 2
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/brain_rep_1.SRR1525405.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "35"

    sleep 0.5
    
    #  brain, rep 1, tech 3
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/brain_rep_1.SRR1525406.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "35"

    sleep 0.5
    
    #  brain, rep 2
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/brain_rep_2.SRR1525407.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "35"

    sleep 0.5
    
    #  spleen, rep 1, combined
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/spleen_rep_1.SRR1525408-09-10.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "35"

    sleep 0.5
    
    #  spleen, rep 1, tech 1
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/spleen_rep_1.SRR1525408.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "35"

    sleep 0.5
    
    #  spleen, rep 1, tech 2
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/spleen_rep_1.SRR1525409.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "35"

    sleep 0.5
    
    #  spleen, rep 1, tech 3
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/spleen_rep_1.SRR1525410.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "35"

    sleep 0.5
    
    #  spleen, rep 2
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/spleen_rep_2.SRR1525411.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "35"

    sleep 0.5
    
    #  ovary, only 1 rep
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/ovary_rep.SRR1525412-13-14.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "35"

    sleep 0.5
    
    #  ovary, rep 1, tech 1
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/ovary_rep.SRR1525412.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "35"

    sleep 0.5
    
    #  ovary, rep 1, tech 2
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/ovary_rep.SRR1525413.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "35"

    sleep 0.5
    
    #  ovary, rep 1, tech 3
    qsub align_Berletch-data_outSAMunmapped-None.sh \
    "${index}" \
    "${d_data}/ovary_rep.SRR1525414.fastq.gz" \
    "${align_to}" \
    "${i}" \
    "35"
done
