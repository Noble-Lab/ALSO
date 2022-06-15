
## News and Updates

* 2022-05-11
  + Cleaned up old example code
  + Will create a pull request to move `03-remove-duplicate-qnames.sh` from the ALSO repo to the Shendure-Lab sciatac pipeline repo after the allele score comparison module is completed
  + Kris will work on the allele score comparison module

* 2022-05-10
  + Kris' new version preprocess module passed tests from both Kris and Gang
  + Gang runs the preprocess module on all samples

* 2022-05-04
  + Gang tested on the one sample from mm10, one sample from CAST
  + Kris tested on the largest bam that we have
  + Bill cleaned up the space of vol6; all future results will be stored in vol6 
  + Update the workflow according to Kris's newest preprocess module (4 steps)

* 2022-05-02
  + Upload/update test code for debugging the preprocess module
  + Add worflow script `03-filter-qname.sh`, an in-progress shell pipeline to handle the preprocessing
  + Update `README` for using the script
  + Write code for handling intermediate files, e.g., deleting, keeping, etc.
  + Further updates, cleanup of the `README`

* 2022-04-13
  + Upload/update test code for debugging the preprocess module
  + Update `README` for information on running the test code

* 2022-04-08
  + Replace the workflow chart
  + Upload code to debug preprocess module
  
* 2022-03-27
  + For `04`, add additional code to remove singletons from split bam files

* 2022-03-26
  + Add additional options, corrections to `04-split-index-repair-bam.sh`
    + "mm10" mode, which does not output POS and MPOS bed files
    + "strain" mode, which outputs POS and MPOS bed files
    + Additional to sort and index bam infile if necessary
  + Update associated test script for new modes

* 2022-03-24
  + Update workflow chart with yellow box (preprocess step)
  + Update run script for preprocess step

* 2022-03-23
  + Add `06-convert-bam-to-df_join-bed_write-rds.R`
  + Clean up repo, removing unneeded scripts and data files
  + Update dependencies listed in `README`

* 2022-03-20
  + Add 05-lift-strain-to-mm10.sh
  + Add script to download and process liftOver chain files: `get-liftOver-chains.sh`
  + Add script to downsample bam files: `generate-downsampled-bam.sh`
  + Minor changes to workflow scripts `01` and `04`
  + Update `README`, including sample-call section

* 2022-03-19
  + Update workflow image
  + Update `README` for (filter reads with MAPQ < 30; then removing singleton; subread repair)
  + Update code for (filter reads with MAPQ < 30; then removing singleton; subread repair.)

* 2022-03-17
  + Add new workflow image
  + CX updated get_unique_fragments.py. Kris will test it on duplicates
  + After Shendure lab pipeline, we will first filter reads with MAPQ < 30; then removing singleton; (Kris: no need to sort anymore) subread repair
