#!/bin/bash


#  Set up variables, home directory, functions --------------------------------
# shellcheck disable=1091
# shellcheck disable=SC2002

dir_base="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
dir_data="${dir_base}/data/files_bam"
# in="Disteche_sample_13.dedup.CAST.bam"
in="Disteche_sample_13.dedup.mm10.bam"
infile="${dir_data}/${in}"

cd "${dir_base}" || 
    {
        echo "Exiting: cd failed. Check on this."
        exit 1
    }

. ./bin/auxiliary/functions-preprocessing.sh


#  Perform preprocessing ------------------------------------------------------
sortBamByCoordinate 4 "${infile}"
indexBam 4 "${infile/.bam/.sort-c.bam}"

removeLowQualityReads 4 "${infile/.bam/.sort-c.bam}"
indexBam 4 "${infile/.bam/.sort-c.rm.bam}"

runFlagstat 4 "${infile/.bam/.sort-c.rm.bam}"

splitBamByChromosome 4 "${infile/.bam/.sort-c.rm.bam}" "chr19"
indexBam 4 "${infile/.bam/.sort-c.rm.chr19.bam}"

splitBamByChromosome 4 "${infile/.bam/.sort-c.rm.bam}" "chr18"
indexBam 4 "${infile/.bam/.sort-c.rm.chr18.bam}"

splitBamByChromosome 4 "${infile/.bam/.sort-c.rm.bam}" "chr17"
indexBam 4 "${infile/.bam/.sort-c.rm.chr17.bam}"

splitBamByChromosome 4 "${infile/.bam/.sort-c.rm.bam}" "chr16"
indexBam 4 "${infile/.bam/.sort-c.rm.chr16.bam}"

splitBamByChromosome 4 "${infile/.bam/.sort-c.rm.bam}" "chr15"
indexBam 4 "${infile/.bam/.sort-c.rm.chr15.bam}"

splitBamByChromosome 4 "${infile/.bam/.sort-c.rm.bam}" "chr14"
indexBam 4 "${infile/.bam/.sort-c.rm.chr14.bam}"

# shellcheck disable=2043
# for i in 14 15 16 17 18 19; do
for i in 18; do
    Rscript ./bin/generate-qname-lists.R \
    --bam "${infile/.bam/.sort-c.rm.chr${i}.bam}" \
    --bai "${infile/.bam/.sort-c.rm.chr${i}.bam.bai}" \
    --outdir "${dir_data}" \
    --chunk 100000 \
    --mated FALSE \
    --unmated TRUE \
    --ambiguous TRUE \
    --trans TRUE \
    --duplicated TRUE \
    --singleton TRUE \
    --tally TRUE \
    --remove TRUE
    
    catQnameList \
    "${infile/.bam/.sort-c.rm.chr${i}.singleton.txt.gz}" \
    "${infile/.bam/.sort-c.rm.chr${i}.duplicated.txt.gz}" \
    "${infile/.bam/.sort-c.rm.chr${i}.combined.txt.gz}"

    excludeQnameReadsPicard \
    "${infile/.bam/.sort-c.rm.chr${i}.bam}" \
    "${infile/.bam/.sort-c.rm.chr${i}.combined.txt.gz}" \
    "${infile/.bam/.sort-c.rm.chr${i}.filtered.bam}"
done

countLinesBam "${infile/.bam/.sort-c.rm.chr18.bam}"  # 3013449
countLinesBam "${infile/.bam/.sort-c.rm.chr18.filtered.bam}"  # 3013440

runFlagstat 4 "${infile/.bam/.sort-c.rm.chr18.bam}"
runFlagstat 4 "${infile/.bam/.sort-c.rm.chr18.filtered.bam}"

h25 "Disteche_sample_13.dedup.mm10.sort-c.rm.chr18.flagstat.txt"
# 3013449 + 0 in total (QC-passed reads + QC-failed reads)
# 3013449 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 3013449 + 0 mapped (100.00% : N/A)
# 3013449 + 0 primary mapped (100.00% : N/A)
# 3013449 + 0 paired in sequencing
# 1506723 + 0 read1
# 1506726 + 0 read2
# 3013449 + 0 properly paired (100.00% : N/A)
# 3013449 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

h25 "Disteche_sample_13.dedup.mm10.sort-c.rm.chr18.filtered.flagstat.txt"
# 3013440 + 0 in total (QC-passed reads + QC-failed reads)
# 3013440 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 3013440 + 0 mapped (100.00% : N/A)
# 3013440 + 0 primary mapped (100.00% : N/A)
# 3013440 + 0 paired in sequencing
# 1506720 + 0 read1
# 1506720 + 0 read2
# 3013440 + 0 properly paired (100.00% : N/A)
# 3013440 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)


#  mm10_all_with-unique-QNAMEs ------------------------------------------------
Rscript ./bin/generate-qname-lists.R \
--bam "${infile/.bam/.sort-c.rm.bam}" \
--bai "${infile/.bam/.sort-c.rm.bam.bai}" \
--outdir "${dir_data}" \
--chunk 100000 \
--mated FALSE \
--unmated TRUE \
--ambiguous TRUE \
--trans TRUE \
--duplicated TRUE \
--singleton TRUE \
--unique TRUE \
--tally TRUE \
--remove TRUE

combineQnameList \
"${infile/.bam/.sort-c.rm.singleton.txt.gz}" \
"${infile/.bam/.sort-c.rm.duplicated.txt.gz}" \
"${infile/.bam/.sort-c.rm.combined.txt.gz}"

excludeQnameReadsPicard \
"${infile/.bam/.sort-c.rm.bam}" \
"${infile/.bam/.sort-c.rm.combined.txt.gz}" \
"${infile/.bam/.sort-c.rm.filtered.bam}"

countLinesBam "${infile/.bam/.sort-c.rm.bam}"  # 85513712
countLinesBam "${infile/.bam/.sort-c.rm.filtered.bam}"  # 85513504

runFlagstat 4 "${infile/.bam/.sort-c.rm.bam}"
runFlagstat 4 "${infile/.bam/.sort-c.rm.filtered.bam}"

h25 "${infile/.bam/.sort-c.rm.flagstat.txt}"
# 85513712 + 0 in total (QC-passed reads + QC-failed reads)
# 85513712 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 85513712 + 0 mapped (100.00% : N/A)
# 85513712 + 0 primary mapped (100.00% : N/A)
# 85513712 + 0 paired in sequencing
# 42756850 + 0 read1
# 42756862 + 0 read2
# 85513712 + 0 properly paired (100.00% : N/A)
# 85513712 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

h25 "${infile/.bam/.sort-c.rm.filtered.flagstat.txt}"
# 85513504 + 0 in total (QC-passed reads + QC-failed reads)
# 85513504 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 85513504 + 0 mapped (100.00% : N/A)
# 85513504 + 0 primary mapped (100.00% : N/A)
# 85513504 + 0 paired in sequencing
# 42756754 + 0 read1
# 42756750 + 0 read2
# 85513504 + 0 properly paired (100.00% : N/A)
# 85513504 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

echo $(( 85513712 - 85513504 ))  # 208

repairBam 4 "${infile/.bam/.sort-c.rm.filtered.bam}"
# Running repair -d -c on Disteche_sample_13.dedup.mm10.sort-c.rm.filtered.bam...
# Finished scanning the input file. Processing unpaired reads.

# All finished in 2.15 minutes
# Total input reads: 85513446 ; Unpaired reads: 0
sortBamByQnameThenFixmate 4 "${infile/.bam/.sort-c.rm.filtered.bam}"

countLinesBam "${infile/.bam/.sort-c.rm.filtered.bam}"  # 85513504
countLinesBam "${infile/.bam/.sort-c.rm.filtered.repair.bam}"  # 85513504
countLinesBam "${infile/.bam/.sort-c.rm.filtered.fixmate.bam}"  # 85513504

runFlagstat 4 "${infile/.bam/.sort-c.rm.filtered.fixmate.bam}"
h25 "${infile/.bam/.sort-c.rm.filtered.fixmate.flagstat.txt}"
# 85513504 + 0 in total (QC-passed reads + QC-failed reads)
# 85513504 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 85513504 + 0 mapped (100.00% : N/A)
# 85513504 + 0 primary mapped (100.00% : N/A)
# 85513504 + 0 paired in sequencing
# 42756754 + 0 read1
# 42756750 + 0 read2
# 85493298 + 0 properly paired (99.98% : N/A)
# 85513504 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 19120 + 0 with mate mapped to a different chr
# 19120 + 0 with mate mapped to a different chr (mapQ>=5)

#NOTE 1/4 So, filtering with the combined.txt.gz file apparently did not remove
#NOTE 2/4 all of the duplicated entries: Why? Is it because, for duplicated
#NOTE 3/4 QNAMEs, I reduced instances of any QNAME to one, and so only one of
#NOTE 4/4 of the duplicates is being filtered out?

decompressGzip "${infile/.bam/.sort-c.rm.singleton.txt.gz}"
decompressGzip "${infile/.bam/.sort-c.rm.duplicated.txt.gz}"
decompressGzip "${infile/.bam/.sort-c.rm.combined.txt.gz}"

cat "${infile/.bam/.sort-c.rm.singleton.txt}"  # 108
cat "${infile/.bam/.sort-c.rm.duplicated.txt}"  # 24
cat "${infile/.bam/.sort-c.rm.combined.txt}"  # 132
cat "${infile/.bam/.sort-c.rm.combined.txt}" | uniq  # 132


#  mm10_all_with-all-QNAMEs ---------------------------------------------------
Rscript ./bin/generate-qname-lists.R -h

Rscript ./bin/generate-qname-lists.R \
--bam "${infile/.bam/.sort-c.rm.bam}" \
--bai "${infile/.bam/.sort-c.rm.bam.bai}" \
--outdir "${dir_data}" \
--chunk 100000 \
--mated FALSE \
--unmated TRUE \
--ambiguous TRUE \
--trans TRUE \
--duplicated TRUE \
--singleton TRUE \
--unique FALSE \
--tally TRUE \
--remove TRUE

unset list
typeset -a list=(
    "${infile/.bam/.sort-c.rm.ambiguous.txt.gz}"
    "${infile/.bam/.sort-c.rm.ambiguous-tally.txt.gz}"
    "${infile/.bam/.sort-c.rm.singleton.txt.gz}"
    "${infile/.bam/.sort-c.rm.singleton-tally.txt.gz}"
    "${infile/.bam/.sort-c.rm.duplicated.txt.gz}"
    "${infile/.bam/.sort-c.rm.duplicated-tally.txt.gz}"
    "${infile/.bam/.sort-c.rm.unmated.txt.gz}"
    "${infile/.bam/.sort-c.rm.unmated-tally.txt.gz}"
)
for i in "${list[@]}"; do decompressGzip "${i}"; done
for i in "${list[@]}"; do echo "$(basename "${i}" .gz): $(countLines "$(basename "${i}" .gz)")"; done
# Disteche_sample_13.dedup.mm10.sort-c.rm.ambiguous.txt:        4
# Disteche_sample_13.dedup.mm10.sort-c.rm.ambiguous-tally.txt:        4
# Disteche_sample_13.dedup.mm10.sort-c.rm.singleton.txt:      108
# Disteche_sample_13.dedup.mm10.sort-c.rm.singleton-tally.txt:      108
# Disteche_sample_13.dedup.mm10.sort-c.rm.duplicated.txt:       92
# Disteche_sample_13.dedup.mm10.sort-c.rm.duplicated-tally.txt:       24
# Disteche_sample_13.dedup.mm10.sort-c.rm.unmated.txt:      112
# Disteche_sample_13.dedup.mm10.sort-c.rm.unmated-tally.txt:      112

excludeQnameReadsPicard \
"${infile/.bam/.sort-c.rm.bam}" \
"${infile/.bam/.sort-c.rm.duplicated.txt.gz}" \
"${infile/.bam/.sort-c.rm.excl-dup.bam}"

excludeQnameReadsPicard \
"${infile/.bam/.sort-c.rm.excl-dup.bam}" \
"${infile/.bam/.sort-c.rm.singleton.txt.gz}" \
"${infile/.bam/.sort-c.rm.excl-dup-sing.bam}"

countLinesBam "${infile/.bam/.sort-c.rm.bam}"  # 85513712
countLinesBam "${infile/.bam/.sort-c.rm.excl-dup.bam}"  # 85513612
countLinesBam "${infile/.bam/.sort-c.rm.excl-dup-sing.bam}"  # 85513504

runFlagstat 4 "${infile/.bam/.sort-c.rm.bam}"
h25 "${infile/.bam/.sort-c.rm.flagstat.txt}"
# 85513712 + 0 in total (QC-passed reads + QC-failed reads)
# 85513712 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 85513712 + 0 mapped (100.00% : N/A)
# 85513712 + 0 primary mapped (100.00% : N/A)
# 85513712 + 0 paired in sequencing
# 42756850 + 0 read1
# 42756862 + 0 read2
# 85513712 + 0 properly paired (100.00% : N/A)
# 85513712 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

runFlagstat 4 "${infile/.bam/.sort-c.rm.excl-dup.bam}"
h25 "${infile/.bam/.sort-c.rm.excl-dup.flagstat.txt}"
# 85513612 + 0 in total (QC-passed reads + QC-failed reads)
# 85513612 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 85513612 + 0 mapped (100.00% : N/A)
# 85513612 + 0 primary mapped (100.00% : N/A)
# 85513612 + 0 paired in sequencing
# 42756800 + 0 read1
# 42756812 + 0 read2
# 85513612 + 0 properly paired (100.00% : N/A)
# 85513612 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

runFlagstat 4 "${infile/.bam/.sort-c.rm.excl-dup-sing.bam}"
h25 "${infile/.bam/.sort-c.rm.excl-dup-sing.flagstat.txt}"
# 85513504 + 0 in total (QC-passed reads + QC-failed reads)
# 85513504 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 85513504 + 0 mapped (100.00% : N/A)
# 85513504 + 0 primary mapped (100.00% : N/A)
# 85513504 + 0 paired in sequencing
# 42756754 + 0 read1
# 42756750 + 0 read2
# 85513504 + 0 properly paired (100.00% : N/A)
# 85513504 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

#NOTE 01  So, filtering with the combined.txt.gz file apparently did not remove
#NOTE 02  all of the duplicated entries: Why? Is it because, for duplicated
#NOTE 03  QNAMEs, I reduced instances of any QNAME to one, and so only one of
#NOTE 04  of the duplicates is being filtered out?

#NOTE 05  I guess the above is not the issue; might it be that the file needs 
#NOTE 06  to be queryname-sorted prior to filtering, and the list of QNAMEs to
#NOTE 07  filter out needs to be lexogrpahically sorted?

#NOTE 08  How does picard FilterSamReads want the file to be sorted prior to
#NOTE 09  filtering?

sortBamByCoordinate 4 "${infile/.bam/.sort-c.rm.excl-dup-sing.bam}"
indexBam 4 "${infile/.bam/.sort-c.rm.excl-dup-sing.sort-c.bam}"

Rscript ./bin/generate-qname-lists.R \
--bam "${infile/.bam/.sort-c.rm.excl-dup-sing.bam}" \
--bai "${infile/.bam/.sort-c.rm.excl-dup-sing.sort-c.bam}" \
--outdir "${dir_data}" \
--chunk 100000 \
--mated FALSE \
--unmated TRUE \
--ambiguous TRUE \
--trans TRUE \
--duplicated TRUE \
--singleton TRUE \
--unique FALSE \
--tally TRUE \
--remove TRUE
#  [1] "Lines in unmated.txt.gz: 4"

#NOTE 10  It seems that sorting is not an issue... There were some records
#NOTE 11  unique to unmated, i.e., neither present in duplicated nor in
#NOTE 12  singleton, that need to be removed too

#NOTE 13  Next step: combine all unique entries in unmated, duplicated, and
#NOTE 14  singleton, then filter those out of the bam

#NOTE 15  If that works, test with a very, very large file, which should answer
#NOTE 16  whether the approach is applicable to all files

#NOTE 17  Need to understand why some reads are marked as unmated but are
#NOTE 18  neither duplicated nor singleton

#NOTE 19  They seem to from mate pairs in which the other member of the pair
#NOTE 20  was filtered out, leaving behind reads that wuould be marked as
#NOTE 21 duplicates if these were data from single-end sequencing

unset list
typeset -a list=(
    Disteche_sample_13.dedup.mm10.sort-c.rm.excl-dup-sing.unmated-tally.txt.gz
    Disteche_sample_13.dedup.mm10.sort-c.rm.excl-dup-sing.unmated.txt.gz
)
for i in "${list[@]}"; do decompressGzip "${i}"; done

retainQnameReadsPicard \
"${infile/.bam/.sort-c.rm.excl-dup-sing.bam}" \
"${infile/.bam/.sort-c.rm.excl-dup-sing.unmated.txt.gz}" \
"${infile/.bam/.sort-c.rm.excl-dup-sing.unmated.bam}"

samtools view "${infile/.bam/.sort-c.rm.excl-dup-sing.unmated.bam}" | h25
# CGCTCCTTGGGATCATGATAATTGAGGATACTCAGCAACC:2219819760 83  chr5    43781764    42  50M =   43781764    -50 GTCCCGCGAGCGCTCGTCCCAGGCCGCGGCCAGCGCGGGCTCCCGGGACC  F,FFFFF:FFFFF::F,FFF:,,FF:FFFFFF,FFFFFFFFFFF,,,,FF  AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YS:i:-6 YT:Z:CP
# CGCTCCTTGGGATCATGATAATTGAGGATACTCAGCAACC:2219819760 83  chr5    43781764    42  50M =   43781764    -50 GTCCCGCGAGCGCTCGTCCCAGGCCGCGGCCAGCGCGGGCTCCCGGGACC  F,FFFFF:FFFFF::F,FFF:,,FF:FFFFFF,FFFFFFFFFFF,,,,FF  AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YS:i:-6 YT:Z:CP
# GTATGGAGGAATCGGCTTGAGAATGAATAATGCTCAACTC:2214537445 147 chr13   56753228    42  50M =   56753228    -50 GGGCCAGCGCGGCACTTCCCGCCGCCTCCGGCCCGCTAGTGGCGTTGAAC  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF:F  AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YS:i:-4 YT:Z:CP
# GTATGGAGGAATCGGCTTGAGAATGAATAATGCTCAACTC:2214537445 147 chr13   56753228    42  50M =   56753228    -50 GGGCCAGCGCGGCACTTCCCGCCGCCTCCGGCCCGCTAGTGGCGTTGAAC  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF:F  AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YS:i:-4 YT:Z:CP
# TGGTTACGCTTGACCATCAGTGGACCAAGGCCAAGAGCAA:1262738469 83  chr15   93519563    42  50M =   93519563    -50 AGGGGGCGGTGCTGAGCCTCGAACTCGCGGCCCCACACGCTGGGCTAGGG  ,FFF:FF::,:,FFF:F:FF:FF:F,,FFF:F:FFFFFF:FFFFFFF::F  AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YS:i:-3 YT:Z:CP
# TGGTTACGCTTGACCATCAGTGGACCAAGGCCAAGAGCAA:1262738469 83  chr15   93519563    42  50M =   93519563    -50 AGGGGGCGGTGCTGAGCCTCGAACTCGCGGCCCCACACGCTGGGCTAGGG  ,FFF:FF::,:,FFF:F:FF:FF:F,,FFF:F:FFFFFF:FFFFFFF::F  AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YS:i:-3 YT:Z:CP
# ATAGGAGTTCGATCATGATATGGACCAAGGAAGGTACTAA:457284034  83  chr17   87798590    42  50M =   87798590    -50 ACCGCGGACAGCGCATCGGTGCCGCAGCCTGCGGGTGGGAGGCTCCAGAC  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YS:i:-3 YT:Z:CP
# ATAGGAGTTCGATCATGATATGGACCAAGGAAGGTACTAA:457284034  83  chr17   87798590    42  50M =   87798590    -50 ACCGCGGACAGCGCATCGGTGCCGCAGCCTGCGGGTGGGAGGCTCCAGAC  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YS:i:-3 YT:Z:CP

#NOTE 22  These should definitely be removed...

combineQnameList \
"${infile/.bam/.sort-c.rm.singleton.txt.gz}" \
"${infile/.bam/.sort-c.rm.duplicated.txt.gz}" \
"${infile/.bam/.sort-c.rm.combined.txt.gz}"

decompressGzip "${infile/.bam/.sort-c.rm.combined.txt.gz}"
decompressGzip "${infile/.bam/.sort-c.rm.unmated.txt.gz}"

#  ----------------------------------------------------------------------------
#  Get the following code standardized and into a function --------------------
#  ----------------------------------------------------------------------------
getUniqueRecords \
"${infile/.bam/.sort-c.rm.unmated.txt}" \
"${infile/.bam/.sort-c.rm.combined.txt}" \
> "${infile/.bam/.sort-c.rm.unmated-unique.txt}"

cat \
"${infile/.bam/.sort-c.rm.combined.txt}" \
"${infile/.bam/.sort-c.rm.unmated-unique.txt}" \
| uniq \
> "${infile/.bam/.sort-c.rm.combined-final.txt}"
#  ----------------------------------------------------------------------------
#  Get the above code standardized and into a function -----------------------
#  ----------------------------------------------------------------------------

excludeQnameReadsPicard \
"${infile/.bam/.sort-c.rm.bam}" \
"${infile/.bam/.sort-c.rm.combined-final.txt}" \
"${infile/.bam/.sort-c.rm.CORRECTED.bam}"

runFlagstat 4 "${infile/.bam/.sort-c.rm.CORRECTED.bam}"
h25 "${infile/.bam/.sort-c.rm.CORRECTED.flagstat.txt}"
# 85513496 + 0 in total (QC-passed reads + QC-failed reads)
# 85513496 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 85513496 + 0 mapped (100.00% : N/A)
# 85513496 + 0 primary mapped (100.00% : N/A)
# 85513496 + 0 paired in sequencing
# 42756748 + 0 read1
# 42756748 + 0 read2
# 85513496 + 0 properly paired (100.00% : N/A)
# 85513496 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

#  Contrast with...
runFlagstat 4 "${infile/.bam/.sort-c.rm.bam}"
h25 "${infile/.bam/.sort-c.rm.flagstat.txt}"
# 85513712 + 0 in total (QC-passed reads + QC-failed reads)
# 85513712 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 85513712 + 0 mapped (100.00% : N/A)
# 85513712 + 0 primary mapped (100.00% : N/A)
# 85513712 + 0 paired in sequencing
# 42756850 + 0 read1
# 42756862 + 0 read2
# 85513712 + 0 properly paired (100.00% : N/A)
# 85513712 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)


# # -----------------------------------------------------------------------------
# list=(
#     "Disteche_sample_13.dedup.CAST.sort-c.rm.unmated.txt.gz"
#     "Disteche_sample_13.dedup.CAST.sort-c.rm.ambiguous.txt.gz"
#     "Disteche_sample_13.dedup.CAST.sort-c.rm.duplicated.txt.gz"
#     "Disteche_sample_13.dedup.CAST.sort-c.rm.singleton.txt.gz"
# )
# for i in "${list[@]}"; do gzip -dk "${i}"; done
#
# unset list2
# typeset list2=(
#     "Disteche_sample_13.dedup.CAST.sort-c.rm.unmated.txt"
#     "Disteche_sample_13.dedup.CAST.sort-c.rm.ambiguous.txt"
#     "Disteche_sample_13.dedup.CAST.sort-c.rm.duplicated.txt"
#     "Disteche_sample_13.dedup.CAST.sort-c.rm.singleton.txt"
# )
#
# for i in "${list2[@]}"; do echo "${i}: $(countLines "${i}")"; echo ""; done
#
# for i in "${list2[@]}"; do
#     for j in "${list2[@]}"; do
#         {
#             echo -e "----------------------------------------"
#             echo -e "${i}: $(countLines "${i}")"
#             echo -e "${j}: $(countLines "${j}")"
#             echo -e "\n"
#
#             echo -e "--------------------"
#             echo "performDiff: ${i} vs. ${j}"
#             performDiff "${i}" "${j}"
#             echo -e "\n"
#
#             echo -e "--------------------"
#             echo -e "performReverseDiff: ${j} in ${i}"
#             performReverseDiff "${i}" "${j}"
#         } \
#         > "${i%.txt}-vs-${j%.txt}.diff.txt"
#     done
# done
#
# rm "diff.txt"
#
#
# # -----------------------------------------------------------------------------
# cat \
# "${infile/.bam/.sort-c.rm.singleton.txt.gz}" \
# "${infile/.bam/.sort-c.rm.duplicated.txt.gz}" \
# > "${infile/.bam/.sort-c.rm.combined.txt.gz}"
#
# gzip -dk "${infile/.bam/.sort-c.rm.combined.txt.gz}"
#
# excludeQnameReadsPicard \
# "${infile/.bam/.sort-c.rm.bam}" \
# "${infile/.bam/.sort-c.rm.combined.txt}" \
# "${infile/.bam/.sort-c.rm.clean.bam}"
#
# countLinesBam "${infile/.bam/.sort-c.rm.bam}"  # 80638866
# countLinesBam "${infile/.bam/.sort-c.rm.clean.bam}"  # 80638668
# countLines "${infile/.bam/.sort-c.rm.combined.txt}"  # 130
# echo $(( 80638668 + 130 ))  # 80638798
# echo $(( 80638866 - 80638668 ))  # 198
