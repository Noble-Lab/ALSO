#!/bin/bash


dir_base="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
dir_data="${dir_base}/data/files_bam"
in="Disteche_sample_13.dedup.CAST.bam"

cd "${dir_base}" || 
    {
        echo "Exiting: cd failed. Check on this."
        exit 1
    }

cd "${dir_data}" || 
    {
        echo "Exiting: cd failed. Check on this."
        exit 1
    }

sortBamByCoordinate 4 "${in}"
indexBam 4 "${in/.bam/.sort-c.bam}"

removeLowQualityReads 4 "${in/.bam/.sort-c.bam}"
indexBam 4 "${in/.bam/.sort-c.rm.bam}"

splitBamByChromosome 4 "${in/.bam/.sort-c.rm.bam}" "chr19"
indexBam 4 "${in/.bam/.sort-c.rm.chr19.bam}"

splitBamByChromosome 4 "${in/.bam/.sort-c.rm.bam}" "chr18"
indexBam 4 "${in/.bam/.sort-c.rm.chr18.bam}"

splitBamByChromosome 4 "${in/.bam/.sort-c.rm.bam}" "chr17"
indexBam 4 "${in/.bam/.sort-c.rm.chr17.bam}"

cd "${dir_base}" || 
    {
        echo "Exiting: cd failed. Check on this."
        exit 1
    }

echo "${in/.bam/.sort-c.rm.chr19.bam}"
echo "${in/.bam/.sort-c.rm.chr19.bam.bai}"

for i in 17 18 19; do
    Rscript ./bin/preprocess-with-R.R \
    -b "${dir_data}/${in/.bam/.sort-c.rm.chr${i}.bam}" \
    -i "${dir_data}/${in/.bam/.sort-c.rm.chr${i}.bam.bai}" \
    -m FALSE \
    -u TRUE \
    -a TRUE \
    -t TRUE \
    -d TRUE \
    -s TRUE \
    -r FALSE \
    -o "${dir_data}"
done

splitBamByChromosome 4 "${in/.bam/.sort-c.rm.bam}" "chr16"
indexBam 4 "${in/.bam/.sort-c.rm.chr16.bam}"

splitBamByChromosome 4 "${in/.bam/.sort-c.rm.bam}" "chr15"
indexBam 4 "${in/.bam/.sort-c.rm.chr15.bam}"

splitBamByChromosome 4 "${in/.bam/.sort-c.rm.bam}" "chr14"
indexBam 4 "${in/.bam/.sort-c.rm.chr14.bam}"

for i in 14 15 16 17 18 19; do
    Rscript ./bin/preprocess-with-R.R \
    -b "${dir_data}/${in/.bam/.sort-c.rm.chr${i}.bam}" \
    -i "${dir_data}/${in/.bam/.sort-c.rm.chr${i}.bam.bai}" \
    -m FALSE \
    -u TRUE \
    -a TRUE \
    -t TRUE \
    -d TRUE \
    -s TRUE \
    -r FALSE \
    -o "${dir_data}"
done

# -b, --bam         bam infile, including path <chr>
# -i, --bai         bam index, including path <chr>
# -c, --chunk       number of records to read into memory at a single
#                   time <even int> [default: 100000]
# -m, --mated       save mated read QNAME list in a txt.gz outfile
#                   <logical> [default: FALSE]
# -u, --unmated     save unmated read QNAME list in a txt.gz outfile
#                   <logical> [default: TRUE]
# -a, --ambiguous   save ambiguous read QNAME list in a txt.gz outfile
#                   <logical> [default: TRUE]
# -t, --trans       save trans read QNAME list in a txt.gz outfile
#                   <logical> [default: TRUE]
# -d, --duplicated  save duplicate QNAME (i.e., QNAME entries > 2) list
#                   in a txt.gz outfile <logical> [default: TRUE]
# -r, --remove      if present, remove outfiles in the outdirectory
#                   prior to generating outfiles <logical> [default:
#                   TRUE]
# -o, --outdir      directory for saving outfile(s), including path
#                   <chr>


# -----------------------------------------------------------------------------
sortBamByCoordinate 4 "${in}"
indexBam 4 "${in/.bam/.sort-c.bam}"

removeLowQualityReads 4 "${in/.bam/.sort-c.bam}"
indexBam 4 "${in/.bam/.sort-c.rm.bam}"

splitBamByChromosome 4 "${in/.bam/.sort-c.rm.bam}" "chr19"
indexBam 4 "${in/.bam/.sort-c.rm.chr19.bam}"

sortBamByQname 4 "${in/.bam/.sort-c.rm.chr19.bam}"

for i in 19; do
    Rscript ./bin/preprocess-with-R.R \
    -b "${dir_data}/${in/.bam/.sort-c.rm.chr${i}.bam}" \
    -i "${dir_data}/${in/.bam/.sort-c.rm.chr${i}.bam.bai}" \
    -m FALSE \
    -u TRUE \
    -a TRUE \
    -t TRUE \
    -d TRUE \
    -s TRUE \
    -r FALSE \
    -o "${dir_data}"
done

#FIXME Does not work
python ./bin/filter-qname.py \
--bam_in "${dir_data}/${in/.bam/.sort-c.rm.chr19.sort-n.bam}" \
--txt "${dir_data}/${in/.bam/.sort-c.rm.chr19.singleton.txt.gz}" \
--bam_out "${dir_data}/${in/.bam/.sort-c.rm.chr19.sort-n.test_rm-dup-QNAME.bam}"
countLinesBam "${dir_data}/${in/.bam/.sort-c.rm.chr19.sort-n.test_rm-dup-QNAME.bam}"  # 2038987

#FIXME Does not work
python ./bin/filter-qname.py \
--bam_in "${dir_data}/${in/.bam/.sort-c.rm.chr19.sort-n.bam}" \
--txt "${dir_data}/${in/.bam/.sort-c.rm.chr19.singleton.txt}" \
--bam_out "${dir_data}/${in/.bam/.sort-c.rm.chr19.sort-n.test-2_rm-dup-QNAME.bam}"
countLinesBam "${dir_data}/${in/.bam/.sort-c.rm.chr19.sort-n.test-2_rm-dup-QNAME.bam}"  # 2038987

excludeQnameReadsPicard \
"${dir_data}/${in/.bam/.sort-c.rm.chr19.bam}" \
"${dir_data}/${in/.bam/.sort-c.rm.chr19.singleton.txt}" \
"${dir_data}/${in/.bam/.sort-c.rm.chr19.test-3_rm-dup-QNAME.bam}"
countLinesBam "${dir_data}/${in/.bam/.sort-c.rm.chr19.test-3_rm-dup-QNAME.bam}"  # 2038984


# -----------------------------------------------------------------------------
Rscript ./bin/generate-qname-lists.R \
-b "${dir_data}/${in/.bam/.sort-c.rm.bam}" \
-i "${dir_data}/${in/.bam/.sort-c.rm.bam.bai}" \
-m FALSE \
-u TRUE \
-a TRUE \
-t TRUE \
-d TRUE \
-s TRUE \
-r FALSE \
-o "${dir_data}"

# [1] "Started: Using Rsamtools to load in 'Disteche_sample_13.dedup.CAST.sort-c.rm.bam' and 'Disteche_sample_13.dedup.CAST.sort-c.rm.bam.bai', and reading various fields such as 'qname' into memory in chunks of 100,000 records per iteration of while loop."
#
# [1] "Counting the number of records in Disteche_sample_13.dedup.CAST.sort-c.rm.bam..."
# [1] "Number of records: 80,638,866"
#
# 200000  300000  400000  500000  600000  700000  800000  900000  1000000  1100000  1200000  1300000  1400000  1500000  1600000  1700000  1800000  1900000  2000000  2100000  2200000  2300000  2400000  2500000  2600000  2700000  2800000  2900000  3000000  3100000  3200000  3300000  3400000  3500000  3600000  3700000  3800000  3900000  4000000  4100000  4200000  4300000  4400000  4500000  4600000  4700000  4800000  4900000  5000000  5100000  5200000  5300000  5400000  5500000  5600000  5700000  5800000  5900000  6000000  6100000  6200000  6300000  6400000  6500000  6600000  6700000  6800000  6900000  7000000  7100000  7200000  7300000  7400000  7500000  7600000  7700000  7800000  7900000  8000000  8100000  8200000  8300000  8400000  8500000  8600000  8700000  8800000  8900000  9000000  9100000  9200000  9300000  9400000  9500000  9600000  9700000  9800000  9900000  10000000  10100000  10200000  10300000  10400000  10500000  10600000  10700000  10800000  10900000  11000000  11100000  11200000  11300000  11400000  11500000  11600000  11700000  11800000  11900000  12000000  12100000  12200000  12300000  12400000  12500000  12600000  12700000  12800000  12900000  13000000  13100000  13200000  13300000  13400000  13500000  13600000  13700000  13800000  13900000  14000000  14100000  14200000  14300000  14400000  14500000  14600000  14700000  14800000  14900000  15000000  15100000  15200000  15300000  15400000  15500000  15600000  15700000  15800000  15900000  16000000  16100000  16200000  16300000  16400000  16500000  16600000  16700000  16800000  16900000  17000000  17100000  17200000  17300000  17400000  17500000  17600000  17700000  17800000  17900000  18000000  18100000  18200000  18300000  18400000  18500000  18600000  18700000  18800000  18900000  19000000  19100000  19200000  19300000  19400000  19500000  19600000  19700000  19800000  19900000  20000000  20100000  20200000  20300000  20400000  20500000  20600000  20700000  20800000  20900000  21000000  21100000  21200000  21300000  21400000  21500000  21600000  21700000  21800000  21900000  22000000  22100000  22200000  22300000  22400000  22500000  22600000  22700000  22800000  22900000  23000000  23100000  23200000  23300000  23400000  23500000  23600000  23700000  23800000  23900000  24000000  24100000  24200000  24300000  24400000  24500000  24600000  24700000  24800000  24900000  25000000  25100000  25200000  25300000  25400000  25500000  25600000  25700000  25800000  25900000  26000000  26100000  26200000  26300000  26400000  26500000  26600000  26700000  26800000  26900000  27000000  27100000  27200000  27300000  27400000  27500000  27600000  27700000  27800000  27900000  28000000  28100000  28200000  28300000  28400000  28500000  28600000  28700000  28800000  28900000  29000000  29100000  29200000  29300000  29400000  29500000  29600000  29700000  29800000  29900000  30000000  30100000  30200000  30300000  30400000  30500000  30600000  30700000  30800000  30900000  31000000  31100000  31200000  31300000  31400000  31500000  31600000  31700000  31800000  31900000  32000000  32100000  32200000  32300000  32400000  32500000  32600000  32700000  32800000  32900000  33000000  33100000  33200000  33300000  33400000  33500000  33600000  33700000  33800000  33900000  34000000  34100000  34200000  34300000  34400000  34500000  34600000  34700000  34800000  34900000  35000000  35100000  35200000  35300000  35400000  35500000  35600000  35700000  35800000  35900000  36000000  36100000  36200000  36300000  36400000  36500000  36600000  36700000  36800000  36900000  37000000  37100000  37200000  37300000  37400000  37500000  37600000  37700000  37800000  37900000  38000000  38100000  38200000  38300000  38400000  38500000  38600000  38700000  38800000  38900000  39000000  39100000  39200000  39300000  39400000  39500000  39600000  39700000  39800000  39900000  40000000  40100000  40200000  40300000  40400000  40500000
# [1] "Lines in unmated.txt.gz: 114"
# [1] "Lines in ambiguous.txt.gz: 2"
# [1] "Lines in duplicated.txt.gz: 20"
# [1] "Lines in singleton.txt.gz: 110"
#
# [1] "Completed: Using Rsamtools to load in 'Disteche_sample_13.dedup.CAST.sort-c.rm.bam' and 'Disteche_sample_13.dedup.CAST.sort-c.rm.bam.bai', reading into memory various fields such as 'qname', then writing out txt.gz files for QNAMEs."

list=(
    "Disteche_sample_13.dedup.CAST.sort-c.rm.unmated.txt.gz"
    "Disteche_sample_13.dedup.CAST.sort-c.rm.ambiguous.txt.gz"
    "Disteche_sample_13.dedup.CAST.sort-c.rm.duplicated.txt.gz"
    "Disteche_sample_13.dedup.CAST.sort-c.rm.singleton.txt.gz"
)
for i in "${list[@]}"; do gzip -dk "${i}"; done

unset list2
typeset list2=(
    "Disteche_sample_13.dedup.CAST.sort-c.rm.unmated.txt"
    "Disteche_sample_13.dedup.CAST.sort-c.rm.ambiguous.txt"
    "Disteche_sample_13.dedup.CAST.sort-c.rm.duplicated.txt"
    "Disteche_sample_13.dedup.CAST.sort-c.rm.singleton.txt"
)

for i in "${list2[@]}"; do echo "${i}: $(countLines "${i}")"; echo ""; done

for i in "${list2[@]}"; do
    for j in "${list2[@]}"; do
        {
            echo -e "----------------------------------------"
            echo -e "${i}: $(countLines "${i}")"
            echo -e "${j}: $(countLines "${j}")"
            echo -e "\n"

            echo -e "--------------------"
            echo "diff: ${i} vs. ${j}"
            diff "${i}" "${j}"
            echo -e "\n"
            
            echo -e "--------------------"
            echo -e "performReverseDiff: ${j} in ${i}"
            performReverseDiff "${i}" "${j}"
        } \
        > "${i%.txt}-vs-${j%.txt}.diff.txt"
    done
done

rm "diff.txt"


# -----------------------------------------------------------------------------
cat \
"${dir_data}/${in/.bam/.sort-c.rm.singleton.txt.gz}" \
"${dir_data}/${in/.bam/.sort-c.rm.duplicated.txt.gz}" \
> "${dir_data}/${in/.bam/.sort-c.rm.combined.txt.gz}"

gzip -dk "${dir_data}/${in/.bam/.sort-c.rm.combined.txt.gz}"

excludeQnameReadsPicard \
"${dir_data}/${in/.bam/.sort-c.rm.bam}" \
"${dir_data}/${in/.bam/.sort-c.rm.combined.txt}" \
"${dir_data}/${in/.bam/.sort-c.rm.clean.bam}"

countLinesBam "${dir_data}/${in/.bam/.sort-c.rm.bam}"  # 80638866
countLinesBam "${dir_data}/${in/.bam/.sort-c.rm.clean.bam}"  # 80638668
countLines "${dir_data}/${in/.bam/.sort-c.rm.combined.txt}"  # 130
echo $(( 80638668 + 130 ))  # 80638798
echo $(( 80638866 - 80638668 ))  # 198
