#!/bin/bash

h_rt="3:0:0"
mfree="2G"
pe_serial=2
queue="noble-short.q"
outpath="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results/2022-0509_test-preprocessing-module"
path_bam_CAST="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-02-23/get_unique_fragments"
path_bam_mm10="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-02-21/get_unique_fragments"
bam="Disteche_sample_13.dedup.bam"
job_name="filter-qnames"
email="kga0@uw.edu"
grid_jobs_max=5

cd "/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross" ||
    {
        echo "cd failed. Check on this."
    }

for i in "${path_bam_CAST}" "${path_bam_mm10}"; do
    while read -r job; do
        grid_jobs_tally="$( checkJobs | grep -c "${job_name}" )"

        while [[ ${grid_jobs_tally} -ge "${grid_jobs_max}" ]]; do
            sleep 5
            printf "."
            grid_jobs_tally="$( checkJobs | grep -c "${job_name}" )"
        done
        echo ""

        echo "# -------------------------------------------------------"
        printf "Job:\n%s\n\n" "${job}"
        printf "Job contents:\n%s\n\n" "$(cat "${job}")"

        qsub \
        -S /bin/bash \
        -l h_rt="${h_rt}" \
        -l mfree="${mfree}" \
        -l gpgpu=FALSE \
        -pe serial "${pe_serial}" \
        -q "${queue}" \
        -cwd \
        -N "${job_name}" \
        -M "${email}" \
        -o "${outpath}" \
        -e "${outpath}" \
        ./bin/workflow/03-filter-qnames.sh \
        -u FALSE \
        -c TRUE \
        -l TRUE \
        -i "${bam}" \
        -o "${outpath}" \
        -f TRUE \
        -r TRUE \
        -p "${pe_serial}"

        printf "Job submission time: %s\n\n" "$(date)"
        sleep 1
    done < <(find "${i}" -maxdepth 1 -name "*.bam" -type f | sort -V)
done
# -h print this help message and exit
# -u use safe mode: "TRUE" or "FALSE" (logical)
# -c use KA conda environment: "TRUE" or "FALSE" (logical)
# -l run on GS HPC: "TRUE" or "FALSE" (logical)
# -i bam infile, including path (chr)
# -o path for outfiles (chr); path will be made if it does not exist
# -f run samtools flagstat on bams: "TRUE" or "FALSE" (logical)
# -r remove intermediate files: "TRUE" or "FALSE" (logical)
# -p number of cores for parallelization (int >= 1); default: 1


# -----------------------------------------------------------------------------
h_rt="12:0:0"
mfree="2G"
pe_serial=2
queue="noble-long.q"
outpath="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results/2022-0509_test-preprocessing-module_mm10"
path_bam_mm10="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-02-21/get_unique_fragments"
bam="Disteche_sample_11.dedup.bam"
job_name="filter-qnames"
email="kga0@uw.edu"
grid_jobs_max=5

cd "/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross" ||
    {
        echo "cd failed. Check on this."
    }

qsub \
-S /bin/bash \
-l h_rt="${h_rt}" \
-l mfree="${mfree}" \
-l gpgpu=FALSE \
-pe serial "${pe_serial}" \
-q "${queue}" \
-cwd \
-N "${job_name}" \
-M "${email}" \
-o "${outpath}" \
-e "${outpath}" \
./bin/workflow/03-filter-qnames-HPC.sh \
-u FALSE \
-c TRUE \
-l TRUE \
-i "${path_bam_mm10}/${bam}" \
-o "${outpath}" \
-f TRUE \
-r TRUE \
-p "${pe_serial}"
