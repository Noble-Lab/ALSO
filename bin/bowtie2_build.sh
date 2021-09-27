#!/bin/bash
#$ -S /bin/bash
#$ -l h_rt=7:59:59
#$ -l mfree=8G
#$ -l gpgpu=FALSE
#$ -pe serial 1
#$ -q noble-short.q
#$ -cwd

#  KA

#  Description:
#+ See printUsage()
#+
#+ Usage:
#+ See printUsage()
#+
#+ Notes:
#+

printUsage() {
    echo ""
    echo "${0}"
    echo "Run bowtie2-build."
    echo ""
    echo "Arguments:"
    echo "-h <print this help message and exit>"
    echo "-f <unzipped fasta, including directory path (chr)>"
    echo "-i <directory path for bowtie2 indices (chr)>"
    echo "-p <prefix for bowtie2 indices (chr)>"
    exit
}

while getopts "h:f:i:p:" opt; do
    case "${opt}" in
        h) printUsage ;;
        f) fasta="${OPTARG}" ;;
        i) indices="${OPTARG}" ;;
        p) prefix="${OPTARG}" ;;
        *) printUsage ;;
    esac
done

[[ -z "${fasta}" ]] && printUsage
[[ -z "${indices}" ]] && printUsage
[[ -z "${prefix}" ]] && printUsage

#  Check dependency
command -v bowtie2 &>/dev/null ||
        {
            echo "Exiting: bowtie2 not found. Install Bowtie2."
            exit 1
        }

#  Check directories
[[ -d "${fasta}" ]] || 
    {
        echo "Exiting: The following directory does not exist: ${fasta}"
        exit 1
    }
[[ -d "${indices}" ]] || mkdir -p "${indices}"

bowtie2-build "${fasta}" "${indices}/${prefix}"
