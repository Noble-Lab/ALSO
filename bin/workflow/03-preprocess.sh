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

#  Check for arguments
if [[ $# -ne 1 ]]; then
    echo ""
    echo "${0}:"
    echo "Take a bam file output by sciatac_pipeline and preprocess it for"
    echo "input into 04-split-index-repair-bam.sh."
    echo ""
    echo "Preprocessing is comprised of the following steps:"
    echo "    1. Filter bam for reads with MAPQ >=30"
    echo "    2. To update flag inormation, sort bam by QNAME and perform"
    echo "       'samtools fixmate'"
    echo "    3. Sort bam by coordinate and then filter with '-f 3 -F 12'"
    echo "    4. To update flag inormation, sort bam by QNAME and perform"
    echo "       'samtools fixmate'"
    echo "Need 1 argument: <RepeatMasker .out.gz>
        <TE threshold (int)> <SINE threshold (int)> <GENCODE .gff.gz>
        <GENCODE .gtf.gz>"
  exit
fi

bam_sciatac="Disteche_sample_6.dedup.bam"
bam_preprocessed="Disteche_sample_6.dedup.preprocessed.bam"

#  Pipeline
samtools view -@ "${parallelize}" -h -b -q 30 "${bam_sciatac}" \
| samtools sort -n -@ "${parallelize}" - \
| samtools fixmate -@ "${parallelize}" - \
| samtools sort -@ "${parallelize}" - \
| samtools view -@ "${parallelize}" -f 3 -F 12 - \
| samtools sort -n -@ "${parallelize}" - \
| samtools fixmate -@ "${parallelize}" - \
| samtools sort -@ "${parallelize}" - -o "${bam_preprocessed}"

## for loop is too slow
## let's  run it in parallel by samples.

sample_id=$1
echo "${sample_id}"
bam_input="CAST-EiJ/get_unique_fragments/Disteche_sample_${sample_id}.dedup.bam"
bam_no_singleton_mapq_gte_30="CAST-EiJ-output/MAPQ30_RmSingleton/Disteche_sample_${sample_id}.RmSingleton.bam"
bam_ready="CAST-EiJ-output/Repair/Disteche_sample_${sample_id}.repair.bam"

echo "Preprocessing:"
samtools view -@ 4 -h -b -f 3 -F 12 -q 30 "${bam_input}" -o "${bam_no_singleton_mapq_gte_30}"
subread/bin/utilities/repair -d -T 4 -c -i "${bam_no_singleton_mapq_gte_30}" -o "${bam_ready}"
echo "Preprocessing module done."
