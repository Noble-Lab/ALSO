#!/usr/bin/env Rscript

#  preprocess-with-R.R
#  KA


#  Functions ------------------------------------------------------------------
assignFilesToVariables <- function(x) {
    # Create a variable based on a file's filename and assign that file to the
    # variable
    #
    # :param x: file, including path (chr)
    # :return filename: variable derived from filename with filename assigned to it
    split <- x %>% strsplit(., "/")
    filename <- split %>% lapply(., `[[`, length(split[[1]]))
    return(filename)
}


checkLibraryAvailability <- function(x) {
    # Check that library is available in environment; return a message if not
    # 
    # :param x:
    ifelse(
        nzchar(system.file(package = as.character(x))),
        "",
        stop(paste0(
            "Library '", x, "' was not found. ",
            "Are you in the correct environment? ",
            "Or do you need to install '", x, "' in the current environment?"
        ))
    )
}


evaluateOperation <- function(x) {
    # Evaluate an operation in string format
    # 
    # :param x: operation string (chr)
    # :return: evaluation of the operation
    return(eval(parse(text = x), envir = .GlobalEnv))
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


makeOperation <- function(x, y) {
    # Create a string for an operation in which a command is assigned to a
    # variable
    # 
    # :param x: variable string (chr)
    # :param y: command string (chr)
    # :return operation: operation string (chr)
    operation <- paste0(x, " <- ", y)
    return(operation)
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

#  Create a parser
ap <- arg_parser(
    name = script,
    description = "
        Script outputs a QNAME txt.gz to be used when filtering with 'samtools
        view -hN'; the user has the option to output txt.gz tables for mated
        and unmated reads.
        
        #TODO Handle duplicate QNAME issue in which a very small number of
        QNAME entries are >2; \filter those out in this script or prior to
        reading in the bam file?
        
        #TODO Leverage Rsamtools::BamFile yieldSize option to iterate through
        large bam files; see, for example, page 14 of
        https://bioconductor.org/packages/devel/bioc/manuals/Rsamtools/man/Rsamtools.pdf
        
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
    short = "-b",
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

#  Parse the command line arguments
directory_base <- "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
directory_data <- "results/kga0/2022-0416-0418_test-preprocessing-module"
directory_out <- paste0(directory_base, "/", directory_data)
bam <- paste0(
    directory_out, "/",
    "Disteche_sample_7.CAST.processed.chr1.bam"
)
bai <- paste0(
    directory_out, "/",
    "Disteche_sample_7.CAST.processed.chr1.bam.bai"
)
mated <- TRUE
unmated <- TRUE
cl <- c(
    #  Arguments for analysis
    "--bam", bam,
    "--bai", bai,
    "--mated", mated,
    "--unmated", unmated,
    "--outdir", directory_out
)
arguments <- parse_args(ap, cl)  # RStudio-interactive work
# arguments <- parse_args(ap)  # Command-line calls

rm(ap, cl, directory_base, directory_data, directory_out, bam, bai, unmated)


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
#  Disteche_sample_7.CAST.processed.chr1.bam; qname, groupid, mate_status: 992,348,800 B
pertinent <- Map(as.data.frame, pertinent) %>%
    as.data.frame() %>%
    tibble::as_tibble(column_name = c("qname", "groupid", "mate_status"))
# pryr::object_size(pertinent)
#  Disteche_sample_7.CAST.processed.chr1.bam; qname, groupid, mate_status: 992,349,224 B

# #  Remove unneeded variables
# rm(file, index, map_params)

# #  Column-bind the standard, AS, and MD fields
# command <- paste0("dplyr::bind_cols(full, tags)")
# operation <- makeOperation(bam, command)
# evaluateOperation(operation)
# # evaluateOperation(paste0("pryr::object_size(", bam, ")"))
# #  Disteche_sample_7.CAST.processed.chr1.bam: 3,716,339,992 B

# #  Remove unneeded variables
# rm(full, tags)
# # pryr::mem_used()
# #  Disteche_sample_7.CAST.processed.chr1.bam: 4.44 GB


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
#  Disteche_sample_7.CAST.processed.chr1.bam; qname only: 882,082,504 B

# readr::write_tsv(
#     pertinent,
#     paste0(
#         arguments$outdir, "/",
#         gsub(".bam", ".mated-QNAMEs.txt", basename(arguments$bam))
#     )
# )

#  Remove duplicate QNAME entries, then write tibble to a gzipped txt file
pertinent <- pertinent[seq_len(nrow(pertinent)) %% 2 == 0, ]
# pryr::object_size(pertinent)
#  Disteche_sample_7.CAST.processed.chr1.bam; only deduplicated qname: 826,948,536 B

readr::write_tsv(
    pertinent,
    paste0(
        arguments$outdir, "/",
        gsub(".bam", ".mated-QNAMEs-dedup.txt", basename(arguments$bam))
    )
)

nrow(pertinent)

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


#  Scraps ---------------------------------------------------------------------
# file <- c(arguments$bam)
# bam <- assignFilesToVariables(file)
# mapply(
#     assign, bam, file, MoreArgs = list(envir = parent.frame())
# )
# 
# file <- c(arguments$bai)
# index <- assignFilesToVariables(file)
# mapply(
#     assign, index, file, MoreArgs = list(envir = parent.frame())
# )

# #  Using Rsamtools, load in standard .bam fields
# command <- paste0(
#     "Rsamtools::BamFile(",
#         bam, ", index = ", index, ", asMates = TRUE",
#     ")", " %>% ",
#     "Rsamtools::scanBam()", " %>% ",
#     "as.data.frame()", " %>% ",
#     "tibble::as_tibble()"
# )
# operation <- makeOperation("full", command)
# evaluateOperation(operation)
# # pryr::object_size(full)
# #  Disteche_sample_7.CAST.processed.chr1.bam: 3,534,162,336 B

# #  Using Rsamtools, load in .bam AS and MD fields
# map_params <- Rsamtools::ScanBamParam(
#     # tag = c(
#     #     "AS", "XS", "XN", "XM", "XO", "XG",
#     #     "NM", "MD", "YS", "YT", "MQ", "MC"
#     # )  #FIXME
#     # what = c("qname", "flag"),
#     tag = c("AS", "MD")
# )
# command <- paste0(
#     bam, " %>% ",
#         "Rsamtools::BamFile(",
#             bam, ", index = ", index, ", asMates = TRUE",
#         ")", " %>% ",
#         "Rsamtools::scanBam(param = map_params)", " %>% ",
#         "as.data.frame()", " %>% ",
#         "tibble::as_tibble()"
# )
# operation <- makeOperation("tags", command)
# evaluateOperation(operation)
# # pryr::object_size(tags)
# #  Disteche_sample_7.CAST.processed.chr1.bam; AS, MD: 182,178,416 B
# #  Disteche_sample_7.CAST.processed.chr1.bam; all tags: #FIXME
# #  Disteche_sample_7.CAST.processed.chr1.bam; qname, flag, AS, MD: 1,119,393,256 B
# 
# #FIXME
# # Warning message:
# # In methods::initRefFields(.Object, classDef, selfEnv, list(...)) :
# #   unnamed arguments to $new() must be objects from a reference class; got an
# #   object of class “character”
