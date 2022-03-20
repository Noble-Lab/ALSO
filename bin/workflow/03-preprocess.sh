#!/bin/bash

# preprocess has three steps:
#     + filter reads with MAPQ < 30; 
#     + then remove singleton; 
#     + subread repair.
    
cd /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-08-samtools ||
    {
        echo "Exiting: Directory not found."
        exit 1
    }
module add samtools/1.14


## mm10
ln -s /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-02-21/ mm10
mkdir ../../results/2022-03-09/
#ln -s /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-03-09/ mm10-output
mv /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-03-09 /net/noble/vol4/noble/user/gangliuw/mouse-cross/results/
ln -s /net/noble/vol4/noble/user/gangliuw/mouse-cross/results/2022-03-09/ /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-08-samtools/mm10-output

mkdir mm10-output/MAPQ30
mkdir mm10-output/RmSingleton
mkdir mm10-output/Repair

#mkdir mm10/sorted
for sample_id in {1..22}
do
    echo "${sample_id}"
    bam_input="mm10/get_unique_fragments/Disteche_sample_${sample_id}.dedup.bam"
    #bam_sort="mm10-output/sorted/Disteche_sample_${sample_id}.dedup.sorted.bam"
    #flagstat="mm10/get_unique_fragments/Disteche_sample_${sample_id}.dedup.flagstat.tsv"
    bam_mapq_gte_30="mm10-output/MAPQ30/Disteche_sample_${sample_id}.MAPQ30.bam"
    bam_no_singleton_mapq_gte_30="mm10-output/RmSingleton/Disteche_sample_${sample_id}.RmSingleton.bam"
    bam_ready="mm10-output/Repair/Disteche_sample_${sample_id}.repair.bam"
    
    echo "Preprocessing:"
    #samtools sort -@ 4 "${bam}" -o "${bam_sort}" > "mm10-output/sorted_${sample_id}.output"
    samtools view -h -b -q 30 -@ 4 "${bam_input}" -o "${bam_mapq_gte_30}"
    samtools view -@ 4 -F 0x08 -b ${bam_mapq_gte_30} > ${bam_no_singleton_mapq_gte_30}
    repair -c -i "${bam_no_singleton_mapq_gte_30}" -o "${bam_ready}"

done


## CAST-EiJ
ls /net/noble/vol5/user/gangliuw/mouse-cross/results/
mv /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-03-08 /net/noble/vol5/user/gangliuw/mouse-cross/results/
cd /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-08-samtools ||
    {
        echo "Exiting: Directory not found."
        exit 1
    }
rm CAST-EiJ-output
ln -s /net/noble/vol5/user/gangliuw/mouse-cross/results/2022-03-08/ CAST-EiJ-output

mkdir CAST-EiJ-output/MAPQ30
mkdir CAST-EiJ-output/RmSingleton
mkdir CAST-EiJ-output/Repair

#mkdir CAST-EiJ-output/sorted
for sample_id in {1..22}
do
    echo "${sample_id}"
    bam_input="CAST-EiJ/get_unique_fragments/Disteche_sample_${sample_id}.dedup.bam"
    #bam_sort="CAST-EiJ-output/sorted/Disteche_sample_${sample_id}.dedup.sorted.bam"
    #flagstat="mm10/get_unique_fragments/Disteche_sample_${sample_id}.dedup.flagstat.tsv"
    bam_mapq_gte_30="CAST-EiJ-output/MAPQ30/Disteche_sample_${sample_id}.MAPQ30.bam"
    bam_no_singleton_mapq_gte_30="CAST-EiJ-output/RmSingleton/Disteche_sample_${sample_id}.RmSingleton.bam"
    bam_ready="CAST-EiJ-output/Repair/Disteche_sample_${sample_id}.repair.bam"
    
    echo "Preprocessing:"
    #samtools sort -@ 4 "${bam}" -o "${bam_sort}" > "CAST-EiJ-output/sorted_${sample_id}.output"
    samtools view -h -b -q 30 -@ 4 "${bam_input}" -o "${bam_mapq_gte_30}"
    samtools view -@ 4 -F 0x08 -b ${bam_mapq_gte_30} > ${bam_no_singleton_mapq_gte_30}
    repair -c -i "${bam_no_singleton_mapq_gte_30}" -o "${bam_ready}"
done
