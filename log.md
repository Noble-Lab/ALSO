
## News and Updates
* 2022-06-04
 +  Pipeline is completed; can call it with `driver_allelic-segregation.sh`
 + Adding instructions for using `driver_allelic-segregation.sh` to `README`.
 + `#TODO` Tets on the HPC and with large files: Do we need to increase memory to the JVM when running `picard`?
 + `#TODO` Clean up messages output by the driver
 + `#TODO` Determine and list all dependencies

* 2022-05-23
 + Addressing error in preprocessing pipeline in which some duplicate QNAMEs persist in processed bam.
 + Adding instructions for using the correction script, `03-remove-duplicate-qnames.sh`.
 + `#TODO` Add corrections in `03-remove-duplicate-qnames.sh` to the initial preprocessing script: `03-filter-problematic-qnames-HPC.sh` `#DONE`

* 2022-05-11
  + Cleaned up the old example code.
  + Will create a pull request for Shendure lab after allele score comparison..
  + Kris will work on the allele score comparison module.

* 2022-05-10
  + Kris' new version preprocess module passed tests from both Kris and Gang.
  + Gang would run the preprocess module on all samples. 

* 2022-05-04
  + Gang tested on the one sample from mm10, one sample from CAST.
  + Kris tesed test on the largest bam that we have.
  + Bill cleaned the space of vol6, and we would store all the future results in vol6.
  + update the workflow according to Kris's newest preprocess module (4 steps).

* 2022-05-02
  + upload/update test code for debugging the preprocess module
  + add worflow script `03-filter-qname.sh`, an in-progress shell pipeline to handle the preprocessing
  + update `README` for using the script
  + write code for handling intermediate files, e.g., deleting, keeping, etc.
  + further updates, cleanup of the `README`

* 2022-04-13
  + upload/update test code for debugging the preprocess module
  + update `README` for information on running the test code

* 2022-04-08
  + replace the workflow chart
  + upload code to debug preprocess module
  
* 2022-03-27
  + for `04`, add additional code to remove singletons from split bam files

* 2022-03-26
  + add additional options, corrections to `04-split-index-repair-bam.sh`
    * "mm10" mode, which does not output POS and MPOS bed files
    * "strain" mode, which outputs POS and MPOS bed files
    * additional to sort and index bam infile if necessary
  + update associated test script for new modes

* 2022-03-24
  + update workflow chart with yellow box (preprocess step)
  + update run script for preprocess step

* 2022-03-23
  + add `06-convert-bam-to-df_join-bed_write-rds.R`
  + clean up repo, removing unneeded scripts and data files
  + update dependencies listed in `README`

* 2022-03-20
  + add 05-lift-strain-to-mm10.sh
  + add script to download and process liftOver chain files: `get-liftOver-chains.sh`
  + add script to downsample bam files: `generate-downsampled-bam.sh`
  + minor changes to workflow scripts `01` and `04`
  + update `README`, including sample-call section

* 2022-03-19
  + update workflow image
  + update `README` for (filter reads with MAPQ < 30; then removing singleton; subread repair)
  + update code for (filter reads with MAPQ < 30; then removing singleton; subread repair.)

* 2022-03-17
  + add new workflow image
  + CX updated get_unique_fragments.py. Kris will test it on duplicates
  + After Shendure lab pipeline, we will first filter reads with MAPQ < 30; then removing singleton; (Kris: no need to sort anymore) subread repair
