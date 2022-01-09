#!/usr/bin/env Rscript

parallel::stopCluster(cl = cluster)

packages <- c(
    "doParallel",
    "foreach",
    "parallel",
    "Rsamtools",
    "tidyverse"
)

for(i in packages) {
    suppressPackageStartupMessages(library(i, character.only = TRUE))
}
rm(i, packages)


# -----------------------------------------------------------------------------
#  ...from read_bam.Bonora.4.R

#  Temporary: Set up work directory (location TB∆)
"/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora/segregatedReads.SNPTHRESH1.Q30" %>% setwd()

#  Files are from
#+ /net/noble/vol2/home/gbonora/proj/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouseDiff_Nmasked/segregatedReads.SNPTHRESH1.Q30
#+ 
#+ For more information, see the following script:
#+ /net/noble/vol2/home/gbonora/proj/2019_sciATAC_analysis/results/gbonora/20191105_sciATAC_mouseDiff_Nmasked/20191105_sciATAC_mouseDiff_Nmasked_allelicSegregation_workflow.sh


#  Set up cluster -------------------------------------------------------------
cores <- 2
cluster <- parallel::makeCluster(
    cores, 
    type = "FORK"
)
doParallel::registerDoParallel(cl = cluster)

# print(cluster)

# foreach::getDoParRegistered()
# foreach::getDoParWorkers()


#  Assign files necessary for subsequent commands -----------------------------
file <- list.files(pattern = ".chr.")


#  Run the pipeline, outputting .rds files for munged tibbles -----------------
foreach::foreach(i = 1:length(file)) %dopar% {
    #  Assign .bam information as list
    map_info <- c("qname", "flag", "rname", "strand", "pos", "mapq", "cigar", "mrnm", "mpos", "isize", "seq")
    tag_info <- c("AS", "XS", "NM")
    map_params <- Rsamtools::ScanBamParam(what = map_info, tag = tag_info)
    intermediate <- Rsamtools::scanBam(file[i], param = map_params)
    
    #  Convert .bam information from list to dataframe to tibble
    intermediate <- intermediate %>% as.data.frame() %>% tibble::as_tibble() %>% tibble::rowid_to_column("ID")
    
    #  Reorder rname, mrnm factor levels
    chromosome <- c(paste0("chr", 1:19), "chrX", "chrY", "chrM")
    intermediate$rname <- forcats::fct_relevel(intermediate$rname, chromosome)
    intermediate$mrnm <- forcats::fct_relevel(intermediate$mrnm, chromosome)
    
    #  Drop unused rname factor levels
    intermediate$rname <- intermediate$rname %>% forcats::fct_drop()
    intermediate$mrnm <- intermediate$mrnm %>% forcats::fct_drop()
    
    #  Create and append pos_end, mpos_end vectors
    #+ 
    #+ Use 49 b/c expression is $pos + 50 exclusive, i.e., "[pos, 50)"
    intermediate$pos_end <- intermediate$pos + 49
    intermediate$mpos_end <- intermediate$mpos + 49
    
    #  Reordering and creation of "read" and "coordinate" columns
    #+ 
    #+ Reorder the pos_end and mpos_end columns; create a column of
    #+ concatenated barcode_sequence strings (i.e., "reads"); and then create a
    #+ column of concatenated barcode_fragment strings (i.e., "coordinates,"
    #+ where a single coordinate is rname, pos, and pos_end)
    #+ 
    #+ #NOTE
    #+ - Terminology change: "b_s" is now "read"
    #+ - Terminology change: "b_f" is now "coordinate"
    intermediate <- intermediate %>%
        dplyr::relocate(pos_end, .after = pos) %>% dplyr::relocate(mpos_end, .after = mpos) %>%
        tidyr::unite(read, c("qname", "seq"), sep = "_", remove = FALSE) %>%
        tidyr::unite(coordinate, c("qname", "rname", "pos", "pos_end"), sep = "_", remove = FALSE)
    
    #  Order columns by coordinate, tag.AS, and mapq
    intermediate <- intermediate[order(intermediate$coordinate, -intermediate$tag.AS, -intermediate$mapq), ]
    
    #  Save tibbles as compressed .rds files
    name_prefix <- file[i] %>% strsplit(., "\\.") %>% lapply(., `[[`, 5) %>% unlist() %>% paste0("GB.", .)
    name_suffix <- file[i] %>% strsplit(., "\\.") %>% lapply(., `[[`, 8) %>% unlist()
    name <- paste0(name_prefix, ".", name_suffix)
    saveRDS(intermediate, file = paste0(name, ".rds"), compress = TRUE)
} %>% invisible()

parallel::stopCluster(cl = cluster)
