#!/bin/bash

# preprocess has two steps:
#    + filter reads with MAPQ < 30 and remove singletons,
#    + perform subread repair to pair mates.
    
cd /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-08-samtools ||
    {
        echo "Exiting: Directory not found."
        exit 1
    }
module add samtools/1.14


## intialize mm10 input directory 
# ln -s /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-02-21/ mm10
# mkdir ../../results/2022-03-09/
# #ln -s /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-03-09/ mm10-output
# mv /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-03-09 /net/noble/vol4/noble/user/gangliuw/mouse-cross/results/
# ln -s /net/noble/vol4/noble/user/gangliuw/mouse-cross/results/2022-03-09/ /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-08-samtools/mm10-output

## create output directory
#mkdir mm10-output/MAPQ30
# mkdir mm10-output/MAPQ30_RmSingleton
# mkdir mm10-output/Repair

## for loop is too slow
## let's parralel by samples.

sample_id=$1
echo "${sample_id}"
bam_input="mm10/get_unique_fragments/Disteche_sample_${sample_id}.dedup.bam"
bam_no_singleton_mapq_gte_30="mm10-output/MAPQ30_RmSingleton/Disteche_sample_${sample_id}.MAPQ30_RmSingleton.bam"
bam_ready="mm10-output/Repair/Disteche_sample_${sample_id}.repair.bam"

echo "Preprocessing:"
# samtools view -@ 4 -h -b -f 3 -F 12 -q 30 "${bam_input}" -o "${bam_no_singleton_mapq_gte_30}"
subread/bin/utilities/repair  -d -T 4 -c -i "${bam_no_singleton_mapq_gte_30}" -o "${bam_ready}"
echo "Preprocessing module done."


# #mkdir mm10/sorted
# for sample_id in {1..22}
# do
#     echo "${sample_id}"
#     bam_input="mm10/get_unique_fragments/Disteche_sample_${sample_id}.dedup.bam"
#     #bam_sort="mm10-output/sorted/Disteche_sample_${sample_id}.dedup.sorted.bam"
#     #flagstat="mm10/get_unique_fragments/Disteche_sample_${sample_id}.dedup.flagstat.tsv"
#     #bam_mapq_gte_30="mm10-output/MAPQ30/Disteche_sample_${sample_id}.MAPQ30.bam"
#     bam_no_singleton_mapq_gte_30="mm10-output/MAPQ30_RmSingleton/Disteche_sample_${sample_id}.MAPQ30_RmSingleton.bam"
#     bam_ready="mm10-output/Repair/Disteche_sample_${sample_id}.repair.bam"
    
#     echo "Preprocessing:"
#     #samtools sort -@ 4 "${bam}" -o "${bam_sort}" > "mm10-output/sorted_${sample_id}.output"
#     #samtools view -h -b -q 30 -@ 4 "${bam_input}" -o "${bam_mapq_gte_30}"
#     #samtools view -@ 4 -F 0x08 -b ${bam_mapq_gte_30} > ${bam_no_singleton_mapq_gte_30}
#     #repair -c -i "${bam_no_singleton_mapq_gte_30}" -o "${bam_ready}"
#     samtools view -@ 4 -h -b -f 3 -F 12 -q 30 "${bam_input}" -o "${bam_no_singleton_mapq_gte_30}"
#     repair -d -T 4 -c -i "${bam_no_singleton_mapq_gte_30}" -o "${bam_ready}"
# done

