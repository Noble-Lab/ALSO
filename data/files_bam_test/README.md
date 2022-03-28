# files_bam_test

```
#!/bin/bash

#  Call from the repo's home directory, 2021_kga0_4dn-mouse-cross
#+ 
#+ User needs an open background connection to UW Nexus

typeset dir_data="./data/files_bam_test"
typeset -a files_bam=(
	Disteche_sample_1.CAST-EiJ.dedup.bam
	Disteche_sample_1.mm10.dedup.bam
)

#  Get files from UW Nexus into a local directory
bash ./bin/get-bams-from-UW-Nexus.sh "${dir_data}"

#  Randomly sample the files
for i in "${files_bam[@]}"; do
	bash ./bin/generate-downsampled-bams.sh \
	-i "${dir_data}/${i}" \
	-o "${dir_data}" \
	-x "${i/.dedup.bam/.downsampled.bam}"
done

#  Remove the initial, unsampled files
for i in "${files_bam[@]}"; do rm -f "${i}"; done

# generate-downsampled-bams.sh:
# -h <print this help message and exit>
# -u <use safe mode: "TRUE" or "FALSE" (logical); default is
#     "FALSE">
# -i <deduplicated bam infile, including path (chr)>
# -o <path for downsampled bam file (chr); path will be made if it
#     does not exist>
# -d <number of paired-end reads to sample down to (even int >= 2);
#     default: 300000>
# -x <prefix for downsampled paired-end bam outfiles (chr, optional);
#     if "-x" is undefined, then prefix is derived from the name of
#     the infile>
# -s <seed number for deterministic random sampling (int >= 1);
#     default: 24>
# -p <number of cores for parallelization (int >= 1); default: 4>
```

`#TODO` Point to `public.html` URL where users can download small bam and bai files
