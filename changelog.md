
# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
- main workflow figure for ALSO.
- readthedocs.


## [0.9.8] - 2022-10-13
### changed
- changelog format from [@GangLiTarheel](https://github.com/Noble-Lab/ALSO/changelog.md).

## [0.9.4] - 2022-06-25
### added
  - Developed Python implementation for filtering bam files by txt-file lists of `QNAMEs`
  - Rewrote/refactored code to handle input for Python-implemented filtering

## [0.9.3] - 2022-06-15
### fixed
  - Update `README.md` and `log.md` files
## [Unreleased]
  - `#TODO` Troubleshoot max memory for JVM when running `picard FilterSamReads`
  - `#TODO` Consolidate shell, R functions into one script for each language
  - `#TODO` Remove harcoded path to GS installation of `picard` from scripts `#DONE`

## [0.9.2] - 2022-06-04
- Pipeline is completed; passed local unit tests with small files; can call it with `driver_allelic-segregation.sh`
- Adding instructions for using `driver_allelic-segregation.sh` to `README`.
## [Unreleased]
- `#TODO` Test on the GS HPC with large files: Do we need to increase max heap memory to the JVM when running `picard`?
- `#TODO` Clean up messages output by the driver
- `#TODO` Determine and list all dependencies

## [0.9.1] - 2022-05-23
- Addressing error in preprocessing pipeline in which some duplicate QNAMEs persist in processed bam.
- Adding instructions for using the correction script, `03-remove-duplicate-qnames.sh`.
- `#TODO` Add corrections in `03-remove-duplicate-qnames.sh` to the initial preprocessing script: `03-filter-problematic-qnames-HPC.sh` `#DONE`


## [0.8.9] - 2022-05-11
### fixed
- Cleaned up old example code from [@Kris] (https://github.com/Noble-Lab/ALSO/XX).

## [Unreleased]
- Will create a pull request to move `03-remove-duplicate-qnames.sh` from the ALSO repo to the Shendure-Lab sciatac pipeline repo after the allele score comparison module is completed
- Kris will work on the allele score comparison module

## [0.8.8] - 2022-05-10
### added
- Kris' new version preprocess module passed tests from both Kris and Gang
- Gang runs the preprocess module on all samples

## [0.8.8] - 2022-05-04
### added
- Gang tested on the one sample from mm10, one sample from CAST
- Kris tested on the largest bam that we have

### fixed
- Bill cleaned up the space of vol6; all future results will be stored in vol6 
- Update the workflow according to Kris's newest preprocess module (4 steps)

## [0.8.7] - 2022-05-02
### added
- Add worflow script `03-filter-qname.sh`, an in-progress shell pipeline to handle the preprocessing
- Write code for handling intermediate files, e.g., deleting, keeping, etc.

### fixed
- Upload/update test code for debugging the preprocess module
- Update `README` for using the script
- Further updates, cleanup of the `README`

## [0.7.9] - 2022-04-13
### fixed
- Upload/update test code for debugging the preprocess module
- Update `README` for information on running the test code

## [0.7.8] - 2022-04-08
### fixed
- Replace the workflow chart
- Upload code to debug preprocess module
  
## [0.7.7] - 2022-03-27
### added
- For `04`, add additional code to remove singletons from split bam files

## [0.7.6] - 2022-03-26
### added
- Add additional options, corrections to `04-split-index-repair-bam.sh`
  - "mm10" mode, which does not output POS and MPOS bed files
  - "strain" mode, which outputs POS and MPOS bed files
  - Additional to sort and index bam infile if necessary
- Update associated test script for new modes

## [0.7.4] - 2022-03-24
### fixed
- Update workflow chart with yellow box (preprocess step)
### added
- Update run script for preprocess step

## [0.7.3] - 2022-03-23
### added
- Add `06-convert-bam-to-df_join-bed_write-rds.R`

## [0.7.2] - 2022-03-24
### removed
- Clean up repo, removing unneeded scripts and data files
### fixed
- Update dependencies listed in `README`

## [0.7.1] - 2022-03-20
### added
- Add 05-lift-strain-to-mm10.sh
- Add script to download and process liftOver chain files: `get-liftOver-chains.sh`
- Add script to downsample bam files: `generate-downsampled-bam.sh`
- Minor changes to workflow scripts `01` and `04`
- Update `README`, including sample-call section

## [0.7.0] - 2022-03-19
### fixed
- Update workflow image
- Update `README` for (filter reads with MAPQ < 30; then removing singleton; subread repair)
- Update code for (filter reads with MAPQ < 30; then removing singleton; subread repair.)

## [0.5.0] - 2022-03-17
### added
- Add new workflow image
### fixed
- CX updated get_unique_fragments.py. Kris will test it on duplicates
- After Shendure lab pipeline, we will first filter reads with MAPQ < 30; then removing singleton; (Kris: no need to sort anymore) subread repair
