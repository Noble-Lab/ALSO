#!/usr/bin/env Rscript

#  filter-bam-per-qnames-list_BamFile.R
#  KA


script <- "filter-bam-per-qnames-list_BamFile.R"
time_start <- Sys.time()
cat(paste0(script, " started: ", time_start, "\n"))


#  Functions ------------------------------------------------------------------
`%notin%` <- Negate(`%in%`)


check_library_availability <- function(x) {
    # Check that library is available in environment; stop and return a message
    # if not
    # 
    # :param x: name of library to check [chr]
    ifelse(
        nzchar(system.file(package = as.character(x))),
        "",
        stop(paste0(
            "Library '", x, "' was not found. Check on this: ",
            "Did you load a module for ", x,"? ",
            "Are you in the correct environment? ",
            "Do you need to install '", x, "' in the current environment?"
        ))
    )
}


convert_time_HMS <- function(x, y) {
    #TODO Description of function
    #
    # :param x: start time [POSIXct]
    # :param y: end time  [POSIXct]
    # :return: #TODO
    dsec <- as.numeric(difftime(y, x, unit = "secs"))
    hours <- floor(dsec / 3600)
    minutes <- floor((dsec - 3600 * hours) / 60)
    seconds <- dsec - (3600 * hours) - (60 * minutes)
    paste0(
        sapply(
            c(hours, minutes, seconds),
            function(x) formatC(x, width = 2, format = "d", flag = "0")
        ),
        collapse = ":"
    )
}


count_records <- function(x) {
    # Count number of records in txt.gz file
    # 
    # :param x: txt.gz file, including path [chr]
    # :return y: number of records in txt.gz file [int]
    y <- system(
        paste0("zcat ", x, " | wc -l"),
        intern = TRUE
    ) %>%
        as.integer()
    return(y)
}


import_library <- function(x) {
    # Suppress messages when loading libraries into the R session
    # 
    # :param x: a vector of libraries [chr]
    invisible(
        suppressWarnings(
            suppressMessages(
                lapply(x, library, character.only = TRUE, quietly = TRUE)
            )
        )
    )
}


stop_if_not_logical <- function(x) {
    #TODO Description of function
    # 
    # :param x: single-element logical vector to evaluate [logical]
    stopifnot(x == TRUE | x == FALSE)
}


#  Source libraries, adjust settings ------------------------------------------
libraries <- c("argparser", "tidyverse", "Rsamtools")
for(i in 1:length(libraries)) check_library_availability(libraries[i])
import_library(libraries)
rm(i, libraries)

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)


#  Parse arguments ------------------------------------------------------------
#  Create a parser
ap <- arg_parser(
    name = script,
    description = "
        Output indexed bam file filtered to retain QNAMEs in a gzipped list.
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
    help = "bam infile, including path [chr]"
)
ap <- add_argument(
    ap,
    short = "-i",
    arg = "--bai",
    type = "character",
    default = NULL,
    help = "bam index, including path [chr]"
)
ap <- add_argument(
    ap,
    short = "-q",
    arg = "--qnames",
    type = "character",
    default = NULL,
    help = "
        gzipped txt-file list of QNAMEs for filtration, including path [chr];
        QNAMEs should be in a single column
    "
)
ap <- add_argument(
    ap,
    short = "-o",
    arg = "--outdir",
    type = "character",
    default = NULL,
    help = "
        directory for saving outfile bam and index, including path [chr];
        temporary split bams are stored in this directory until removed at the
        conclusion of the script
    "
)
ap <- add_argument(
    ap,
    short = "-r",
    arg = "--retain",
    type = "logical",
    default = TRUE,
    help = "
        generate bam outfile that retains only QNAMEs in gzipped txt-file
        [logical; default: TRUE]; if FALSE, then generate bam outfile that
        excludes all QNAMEs in gzipped txt-file
    "
)
ap <- add_argument(
    ap,
    short = "-s",
    arg = "--strain",
    type = "character",
    default = NULL,
    help = "strain name associated with bam infile [chr]"
)
ap <- add_argument(
    ap,
    short = "-c",
    arg = "--chunk",
    type = "integer",
    default = 100000,
    help = "
        number of gzipped txt-file QNAME records to read into memory at a
        single time [even int > 0; default: 100000]
    "
)
ap <- add_argument(
    ap,
    short = "-p",
    arg = "--progress",
    type = "logical",
    default = FALSE,
    help = "
        show progress bar [logical; default: FALSE]
    "
)


#  Parse the arguments --------------------------------------------------------
test_in_RStudio <- FALSE  # Hardcode TRUE for interactive testing; FALSE: CLI
if(isTRUE(test_in_RStudio)) {
    #  RStudio-interactive work
    dir_base <- "."
    dir_data <- "data"
    dir_in <- paste0(dir_base, "/", dir_data, "/", "files_bam_test")
    dir_out <- paste0(dir_base, "/", dir_data, "/", "files_bam_test")
    bam <- paste0(dir_in, "/CAST/Disteche_sample_11.dedup.corrected.corrected.downsample-6000000.bam")
    bai <- paste0(dir_in, "/CAST/Disteche_sample_11.dedup.corrected.corrected.downsample-6000000.bam.bai")
    qnames <- paste0(dir_in, "/Disteche_sample_11.CAST.combined.txt.gz")
    # qnames <- paste0(dir_in, "/Disteche_sample_11.mm10.combined.txt.gz")
    retain <- TRUE
    strain <- "CAST"
    # strain <- "mm10"
    chunk <- 1000000
    progress <- TRUE
    cl <- c(
        "--bam", bam,
        "--bai", bai,
        "--qnames", qnames,
        "--retain", retain,
        "--strain", strain,
        "--outdir", dir_out,
        "--chunk", chunk,
        "--progress", progress
    )
    arguments <- parse_args(ap, cl)
    rm(
        ap, bam, bai, cl, chunk,
        dir_base, dir_data, dir_in, dir_out,
        progress, qnames, retain, strain
    )
} else if(isFALSE(test_in_RStudio)) {
    #  Command-line calls
    arguments <- parse_args(ap)
    rm(ap)
} else {
    stop(paste0(
        "Stopping: Variable 'test_in_RStudio' is not properly set. ",
        "'test_in_RStudio' should be hardcoded as either TRUE or FALSE."
    ))
}
rm(test_in_RStudio)


#  Check on the arguments that were supplied ----------------------------------
stopifnot(file.exists(arguments$bam))
stopifnot(file.exists(arguments$bai))
stopifnot(file.exists(arguments$qnames))
stop_if_not_logical(arguments$retain)
if(arguments$strain == "") stop("--strain is an empty string.")
if(is.na(arguments$strain)) stop("--strain is NA.")
stopifnot(arguments$chunk != 0)
stopifnot(arguments$chunk %% 1 == 0)
stopifnot(arguments$chunk %% 2 == 0)
stop_if_not_logical(arguments$progress)

#  If it does not exist, then create outfile directory 
dir.create(file.path(arguments$outdir), showWarnings = FALSE)


#  Establish prefix to be used in name of bam outfile -------------------------
#  Add strain name to the bam if it's not already present
if(isFALSE(grepl(arguments$strain, basename(arguments$bam), fixed = TRUE))) {
    outstring <- paste0(".", arguments$strain, ".bam")
    bam_prefix_tmp <- gsub(".bam", outstring, basename(arguments$bam))
    rm(outstring)
    
    #  Add strain name for the gzipped txt file of QNAMEs to the bam outfile
    txt_vector <- sub(".txt.gz", "", basename(arguments$qnames)) %>%
        stringr::str_split("\\.") %>%
        unlist()
    txt_label <- txt_vector[2]
    bam_prefix <- paste0(sub(".bam", "", bam_prefix_tmp), ".", txt_label)
    rm(bam_prefix_tmp, txt_label, txt_vector)
} else {
    bam_prefix <- sub(".bam", "", basename(arguments$bam))
}


#  Tally and record the total number of records in the file of qnames ---------
cat(paste0(
    "Counting the number of records in **", basename(arguments$qnames), "**\n"
))
rec_chunk_txt_qnames <- as.integer(arguments$chunk)
rec_total_txt_qnames <- as.integer(system(
    paste0("zcat ", arguments$qnames, " | wc -l"),
    intern = TRUE
))
cat(paste0("Number of records: ", scales::comma(rec_total_txt_qnames), "\n\n"))


#  Iterate through file of qnames in chunks and write outfiles ----------------
cat(paste0(
    "Started: Filtering **", basename(arguments$bam),
    "** to retain reads with QNAME records listed in **",
    basename(arguments$qnames), "**\n"
))

#  Open connection to txt file of qnames files
connection <- gzfile(arguments$qnames, "r")

#  Open connection to bam file, which will be read in in chunks ---------------
bam <- Rsamtools::BamFile(arguments$bam, index = arguments$bai, asMates = TRUE)
Rsamtools::yieldSize(bam) <- arguments$chunk
open(bam)
#HERE Continue here; need to add nested for loop

#  Set up and run 'for loop'
iteration_txt_qnames <- as.integer(ceiling(
    rec_total_txt_qnames / rec_chunk_txt_qnames
))
#HERE E.g., need ceiling calculation for bam here too

if(isTRUE(arguments$progress)) {
    bar <- utils::txtProgressBar(
        min = 0, max = iteration_txt_qnames, initial = 0, style = 3
    )
}

#  Let the user know what they're in for...
cat(paste0(
    "The script will iterate through **", basename(arguments$qnames),
    "** ", iteration_txt_qnames, " times to filter **", basename(arguments$bam), "**"
))

bam_vector <- vector(mode = "character", length = iteration_txt_qnames)
for(i in 1:iteration_txt_qnames) {
    #  For example, see support.bioconductor.org/p/76906/
    if(isFALSE(arguments$progress)) {
        cat(paste0(
            "Iteration through ", basename(arguments$qnames), ": ",
            i, "/", iteration_txt_qnames, "\n"
        ))
    }
    
    #  Establish name of temporary bam outfile for this iteration
    bam_outfile <- paste0(
        arguments$outdir, "/", bam_prefix, ".",
        i, "-of-", iteration_txt_qnames, ".bam"
    )
    
    #  Read in filtration qnames for this iteration
    wanted_qnames <- read.delim(
        connection,
        nrows = arguments$chunk,
        header = FALSE
    )
    
    #  Set up filter rules needed by Rsamtools::filterBam()
    if(isTRUE(arguments$retain)) {
        filter <- FilterRules(list(
            retain_qname = function(x) x$qname %in% unlist(wanted_qnames)
        ))
    } else {
        filter <- FilterRules(list(
            retain_qname = function(x) x$qname %notin% unlist(wanted_qnames)
        ))
    }
    
    #  Filter bam file for records in 'wanted_qnames', writing out a temporary
    #+ bam file to be subsequently merged
    Rsamtools::filterBam(
        arguments$bam,
        bam_outfile,
        filter = filter,
        param = Rsamtools::ScanBamParam(what = "qname")
    )
    
    bam_vector[i] <- bam_outfile
    if(isTRUE(arguments$progress)) utils::setTxtProgressBar(bar, i)
}

#  Specify path and name for final merged outfile
bam_outfile <- paste0(
    sub(paste0(".1-of-", iteration_txt_qnames, ".bam"), "", bam_vector[1]),
    ".bam"
)

#  Merge iteration outfiles into final merged outfile
cat(paste0(
    "Merging ", iteration_txt_qnames,
    " outfiles into a single, final outfile: ", basename(bam_outfile),
    ". Indexing the final outfile too.\n"
))
Rsamtools::mergeBam(
    files = bam_vector,
    destination = bam_outfile,
    indexDestination = TRUE
)


#  Clean up: Remove the iteration outfiles, which are no longer needed --------
if(isFALSE(file.exists(bam_outfile))) {
    stop(paste0("WARNING: File merged failed. Check on this."))
} else if(is.na(file.exists(bam_outfile))) {
    stop(paste0("WARNING: File merged failed. Check on this."))
} else if(isTRUE(file.exists(bam_outfile))) {
    to_remove <- list.files(
        path = arguments$outdir,
        pattern = paste0(bam_prefix, ".*-of-", iteration_txt_qnames, ".bam"),
        full.names = TRUE
    )
    cat(
        paste0("Removing the following temporary files:\n   "),
        paste0("- ", basename(to_remove), "\n   ")
    )
    invisible(suppressWarnings(
        suppressMessages(file.remove(to_remove))
    ))
}


#  End the script -------------------------------------------------------------
cat(paste0(
    "Completed: Filtering **", basename(arguments$bam),
    "** to retain reads with QNAME records listed in **",
    basename(arguments$qnames), "**\n"
))

time_end <- Sys.time()
cat("\n")
cat(paste0(script, " completed: ", time_end, "\n"))
print(convert_time_HMS(time_start, time_end))
cat("\n\n")
rm(time_start, time_end)
