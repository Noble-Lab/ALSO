
#!/bin/bash

cd /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-08-samtools
module add samtools/1.14


## mm10
ln -s /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-02-21/ mm10
mkdir ../../results/2022-03-09/
#ln -s /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-03-09/ mm10-output
mv /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-03-09 /net/noble/vol4/noble/user/gangliuw/mouse-cross/results/
ln -s /net/noble/vol4/noble/user/gangliuw/mouse-cross/results/2022-03-09/ /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-08-samtools/mm10-output

mkdir CAST-EiJ-output/sorted
for sample_id in {1..22}
do
    echo $sample_id
    bam="CAST-EiJ/get_unique_fragments/Disteche_sample_${sample_id}.dedup.bam"
    bam_sort="CAST-EiJ-output/sorted/Disteche_sample_${sample_id}.dedup.sorted.bam"
    #flagstat="mm10/get_unique_fragments/Disteche_sample_${sample_id}.dedup.flagstat.tsv"
    samtools sort -@ 4 "${bam}" -o "${bam_sort}" > CAST-EiJ-output/sorted_${sample_id}.output
done


## CAST-EiJ
ls /net/noble/vol5/user/gangliuw/mouse-cross/results/
mv /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-03-08 /net/noble/vol5/user/gangliuw/mouse-cross/results/
cd /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-08-samtools
rm CAST-EiJ-output
ln -s /net/noble/vol5/user/gangliuw/mouse-cross/results/2022-03-08/ CAST-EiJ-output

mkdir CAST-EiJ-output/sorted
for sample_id in {1..22}
do
    echo $sample_id
    bam="CAST-EiJ/get_unique_fragments/Disteche_sample_${sample_id}.dedup.bam"
    bam_sort="CAST-EiJ-output/sorted/Disteche_sample_${sample_id}.dedup.sorted.bam"
    #flagstat="mm10/get_unique_fragments/Disteche_sample_${sample_id}.dedup.flagstat.tsv"
    samtools sort -@ 4 "${bam}" -o "${bam_sort}" > CAST-EiJ-output/sorted_${sample_id}.output
done
