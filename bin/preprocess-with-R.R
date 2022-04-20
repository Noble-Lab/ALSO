#!/usr/bin/env Rscript

#  preprocess-with-R.R
#  KA


#  Functions ------------------------------------------------------------------
checkLibraryAvailability <- function(x) {
    # Check that library is available in environment; return a message if not
    # 
    # :param x: name (string) for library to check (chr)
    # :return: stop, print message if nzchar() returns empty character vector
    ifelse(
        nzchar(system.file(package = as.character(x))),
        "",
        return(stop(paste0(
            "Library '", x, "' was not found. ",
            "Are you in the correct environment? ",
            "Or do you need to install '", x, "' in the current environment?"
        )))
    )
}


importLibrary <- function(x) {
    # Suppress messages when loading libraries into the R session
    # 
    # :param x: a vector of libraries <chr>
    invisible(
        suppressWarnings(
            suppressMessages(
                lapply(x, library, character.only = TRUE, quietly = TRUE)
            )
        )
    )
}


#  Source libraries, adjust settings ------------------------------------------
libraries <- c("argparser", "pryr", "Rsamtools", "tidyverse")
for(i in 1:length(libraries)) { checkLibraryAvailability(libraries[i]) }
importLibrary(libraries)
rm(i, libraries)

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)


#  Parse arguments ------------------------------------------------------------
script <- "preprocess-with-R.R"
test_in_RStudio <- TRUE  # Hardcode T for testing in RStudio; F for CLI

#  Create a parser
ap <- arg_parser(
    name = script,
    description = "
        Script outputs a QNAME txt.gz to be used when filtering with 'samtools
        view -hN'; the user has the option to output txt.gz tables for mated
        and unmated reads.
        
        #TODO Handle duplicate QNAME issue in which a very small number of
        QNAME entries are >2; \filter those out in this script or prior to
        reading in the bam file? #ANSWER Prior to reading in the bam file.
        
        #TODO Leverage Rsamtools::BamFile yieldSize option to iterate through
        large bam files; see, for example, page 14 of
        bioconductor.org/packages/devel/bioc/manuals/Rsamtools/man/Rsamtools.pdf
        
        See also https://rdrr.io/bioc/GenomicFiles/man/reduceByYield.html
    ",
    hide.opts = TRUE
)

#  Add command line arguments
ap <- add_argument(
    ap,
    short = "-b",
    arg = "--bam",
    type = "character",
    default = NULL,
    help = "bam infile, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-i",
    arg = "--bai",
    type = "character",
    default = NULL,
    help = "bam index, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-m",
    arg = "--mated",
    type = "logical",
    default = FALSE,
    help = "
        save mated read 'qname', 'groupid', and 'mate_status' table in a txt.gz
        outfile <logical>
    "
)
ap <- add_argument(
    ap,
    short = "-u",
    arg = "--unmated",
    type = "logical",
    default = FALSE,
    help = "
        save unmated read 'qname', 'groupid', and 'mate_status' table in a
        txt.gz outfile <logical>
    "
)
ap <- add_argument(
    ap,
    short = "-o",
    arg = "--outdir",
    type = "character",
    default = NULL,
    help = "directory for saving outfile(s), including path <chr>"
)


#  Parse the command line arguments -------------------------------------------
if(isTRUE(test_in_RStudio)) {
    #  RStudio-interactive work
    dir_base <- "/Users/kalavattam/Dropbox/UW/projects-etc"
    dir_proj <- paste0(dir_base, "/", "2021_kga0_4dn-mouse-cross")
    dir_data <- "results/kga0/2022-0416-0418_test-preprocessing-module"
    dir_in_out <- paste0(dir_proj, "/", dir_data)
    bam <- "Disteche_sample_7.CAST.processed.chr1.bam"
    bai <- "Disteche_sample_7.CAST.processed.chr1.bam.bai"
    mated <- TRUE
    unmated <- TRUE
    cl <- c(
        #  Arguments for analysis
        "--bam", paste0(dir_in_out, "/", bam),
        "--bai", paste0(dir_in_out, "/", bai),
        "--mated", mated,
        "--unmated", unmated,
        "--outdir", dir_in_out
    )
    arguments <- parse_args(ap, cl)  
    rm(
        ap, cl,
        dir_base, dir_data, dir_in_out,
        bam, bai, mated, unmated
    )
} else if(isFALSE(test_in_RStudio)) {
    #  Command-line calls
    arguments <- parse_args(ap)
    rm(ap, cl)
} else {
    stop("Stopping: Variable 'test_in_RStudio' is not properly.")
}


#  Check that files exist -----------------------------------------------------
stopifnot(file.exists(arguments$bam))
stopifnot(file.exists(arguments$bai))

#  If outfile directory  does not exist, then create it
dir.create(file.path(arguments$outdir), showWarnings = FALSE)


#  Load in .bam information, including mate information -----------------------
bam <- Rsamtools::BamFile(arguments$bam, index = arguments$bai, asMates = TRUE)
pertinent <- Rsamtools::scanBam(
    bam, param = ScanBamParam(what = c("qname", "groupid", "mate_status"))
)
# pryr::object_size(pertinent)
#  Disteche_sample_7.CAST.processed.chr1.bam; qname, groupid, mate_status: 992.3 MB
pertinent <- Map(as.data.frame, pertinent) %>%
    as.data.frame() %>%
    tibble::as_tibble(column_name = c("qname", "groupid", "mate_status"))
# pryr::object_size(pertinent)
#  Disteche_sample_7.CAST.processed.chr1.bam; qname, groupid, mate_status: 992.3 MB


#  Separate "unmated" reads from "mated" reads --------------------------------
if(isTRUE(arguments$unmated)) {
    #  Split unmated reads into their own tibble
    unmated <- pertinent[pertinent$mate_status == "unmated", ]
    
    #  Write tibble to gzipped txt file
    readr::write_tsv(
        unmated,
        paste0(
            arguments$outdir, "/",
            gsub(".bam", ".unmated-full.txt.gz", basename(arguments$bam))
        )
    )
}

#  Remove reads with $mate_status of "unmated"
pertinent <- pertinent[pertinent$mate_status != "unmated", ]

if(isTRUE(arguments$mated)) {
    #  Write full 'pertinent' tibble to a gzipped txt file
    readr::write_tsv(
        pertinent,
        paste0(
            arguments$outdir, "/",
            gsub(".bam", ".mated-full.txt.gz", basename(arguments$bam))
        )
    )
}

#  Remove all fields except QNAME
pertinent <- dplyr::select(pertinent, -c(groupid, mate_status))
# pryr::object_size(pertinent)
#  Disteche_sample_7.CAST.processed.chr1.bam; qname only: 882.1 MB

#  Remove duplicate QNAME entries, then write tibble to a gzipped txt file
pertinent <- pertinent[seq_len(nrow(pertinent)) %% 2 == 0, ]
# pryr::object_size(pertinent)
#  Disteche_sample_7.CAST.processed.chr1.bam; only deduplicated qname: 826.9 MB

#  Write out final txt file containing QNAMEs to retain in the preprocessed bam
#+ file
outfile <- gsub(".bam", ".mated-QNAMEs-dedup.txt", basename(arguments$bam))
readr::write_tsv(
    pertinent,
    paste0(arguments$outdir, "/", outfile)
)


#  In shell... ----------------------------------------------------------------
# bam="Disteche_sample_7.CAST.processed.chr1.bam"
# keep="Disteche_sample_7.CAST.processed.chr1.mated-QNAMEs-dedup.txt"
# samtools view -@ 4 -hN "${keep}" "${bam}" | samtools view -@ 4 -b - > "${bam/.bam/.mated-QNAMEs-dedup.bam}"
# 
# # keep="Disteche_sample_7.CAST.processed.chr1.mated-QNAMEs.txt"
# # samtools view -@ 4 -hN "${keep}" "${bam}" | samtools view -@ 4 -b - > "${bam/.bam/.mated-QNAMEs.bam}"
# #NOTE Whether *.mated-QNAMEs-dedup.txt or *.mated-QNAMEs.txt, the final
#       output is the same, so go with the option that uses less memory and
#       saves more space
