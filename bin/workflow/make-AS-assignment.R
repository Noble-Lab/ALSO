#!/usr/bin/env Rscript

#  generate-assignment-lists.R
#  KA


script <- "generate-assignment-lists.R"
time_start <- Sys.time()
cat(paste0(script, " started: ", time_start, "\n"))


#  Functions ------------------------------------------------------------------
check_library_availability <- function(x) {
    # Check that library is available in environment; stop and return a message
    # if not
    # 
    # :param x: name of library to check (chr)
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
    # #TODO Description of function
    #
    # :param x: start time (POSIXct)
    # :param y: end time  (POSIXct)
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
    # :param x: txt.gz file, including path (chr)
    # :return y: number of records in txt.gz file (int)
    y <- system(paste0("zcat ", x, " | wc -l"), intern = TRUE) %>%
        as.integer()
    return(y)
}


import_library <- function(x) {
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


make_assignment <- function(x, y) {
    #TODO Description of function
    #
    # :param x: tibble comprised of the following variables: qname (chr), AS
    #           for sample 1 (int), AS for sample 2 (int)
    # :param y: threshold for assignments (int)
    x %>% 
        dplyr::mutate(difference = .[[2]] - .[[3]]) %>% 
        dplyr::mutate(
            assignment = case_when(
                difference >= (-1 * y) & difference <= y ~ "ambiguous",
                difference > y ~ gsub("AS\\.", "", colnames(x[2])),
                difference < (-1 * y) ~ gsub("AS\\.", "", colnames(x[3]))
            )
        )
}


read_line <- function(x, y, z, a) {
    #TODO Description of function
    #
    # :param x: txt/txt.gz file, including path (chr)
    # :param y: line number preceding line of interest (int)
    # :param z: identifier for sample 1 (chr)
    # :param a: identifier for sample 2 (chr)
    # :return: #TODO
    b <- read.delim(
        x,
        header = FALSE,
        skip = y,
        nrows = 1
    ) %>%
        tibble::as_tibble() %>%
        dplyr::rename("qname" = "V1") %>%
        dplyr::rename(!!paste0("AS.", z) := "V2") %>%
        dplyr::rename(!!paste0("AS.", a) := "V3")
    return(b)
}


remove_outfiles <- function(x, y, z, a) {
    # Remove any outfiles present in directory
    # 
    # :param x: directory from which to remove outfiles
    # :param y: prefix for outfiles
    # :param z: string for sample #1
    # :param a: string for sample #2
    system(paste0(
        "rm -f -- ",
        x, "/", y, ".assign-ambiguous.txt ",
        x, "/", y, ".assign-ambiguous.txt.gz ",
        x, "/", y, ".assign-", z,".txt ",
        x, "/", y, ".assign-", z,".txt.gz ",
        x, "/", y, ".assign-", z,"-only.txt ",
        x, "/", y, ".assign-", z,"-only.txt.gz ",
        x, "/", y, ".assign-", a,".txt ",
        x, "/", y, ".assign-", a,".txt.gz ",
        x, "/", y, ".assign-", a,"-only.txt ",
        x, "/", y, ".assign-", a,"-only.txt.gz"
    ))
}


stop_if_not_logical <- function(x) {
    # #TODO Description of function
    # 
    # :param x: single-value logical vector to evaluate (logical)
    stopifnot(x == TRUE | x == FALSE)
}


write_list <- function(x, y, z, a) {
    #TODO Description of function
    #
    # :param x: dataframe or tibble
    # :param y: outdirectory, including path (chr)
    # :param z: outfile prefix (chr)
    # :param a: sample string (chr)
    readr::write_tsv(
        a,
        file = paste0(y, "/", z, ".", a, ".txt.gz"),
        append = TRUE
    )
}


#  Source libraries, adjust settings ------------------------------------------
libraries <- c("argparser", "tidyverse")
for(i in 1:length(libraries)) check_library_availability(libraries[i])
import_library(libraries)
rm(i, libraries)

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)


#  Parse arguments ------------------------------------------------------------
#  Create a parser
ap <- arg_parser(
    name = script,
    description = "",
    hide.opts = TRUE
)

#  Add command line arguments
ap <- add_argument(
    ap,
    short = "-i",
    arg = "--intersection",
    type = "character",
    default = NULL,
    help = "
        tab-separated txt file (gzipped or not), including path, containing
        intersecting query names (QNAME) and alignment scores (AS) for two
        samples <chr>

        For example, the first five rows of a given txt file appears as
        follows, where column 1 is the QNAME, column 2 is the AS for sample 1,
        and column 3 is the AS for sample 3:
        AACCAATAACAAGTCAACGGAATATAGCGGCAAGTCAGCA:1417636611 0   -5
        AACCAATAACAAGTCAACGGAATATAGCGGTTATAATCAA:1988125074 0   0
        AACCAATAACAAGTCAACGGACGCTTCGTCCGAGATAAGA:2423713332 0   -5
        AACCAATAACAAGTCAACGGACGCTTCGTCCTCTGCGATC:1095498111 0   -5
        AACCAATAACAAGTCAACGGACGTACTAGGAAGTTCGCTG:353515658  0   0
    "
)
ap <- add_argument(
    ap,
    short = "-s",
    arg = "--strain_1",
    type = "character",
    default = NULL,
    help = "
        strain name for sample 1; should be the same as the column 2 header in
        the tab-separated txt intersection file <chr>
    "
)
ap <- add_argument(
    ap,
    short = "-u",
    arg = "--strain_2",
    type = "character",
    default = NULL,
    help = "
        strain name for sample 2; should be the same as the column 3 header in
        the tab-separated txt intersection file <chr>
    "
)
ap <- add_argument(
    ap,
    short = "-o",
    arg = "--outdir",
    type = "character",
    default = NULL,
    help = "directory for saving txt.gz outfiles, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-p",
    arg = "--outprefix",
    type = "character",
    default = NULL,
    help = "prefix for outfiles <chr>"
)
ap <- add_argument(
    ap,
    short = "-t",
    arg = "--threshold",
    type = "integer",
    default = 0,
    help = "
        alignment score threshold; the absolute value of the difference in
        alignment scores between samples 1 and 2 must be greater than this
        value in order for a strain-specific assignment to be made; if not, the
        assignment will be \"ambiguous\" (int > 0)

        For example...
    "
)
ap <- add_argument(
    ap,
    short = "-c",
    arg = "--chunk",
    type = "integer",
    default = 1000000,
    help = "number of records to read into memory at a single time (int ≥ 0)"
)
ap <- add_argument(
    ap,
    short = "-r",
    arg = "--remove",
    type = "logical",
    default = TRUE,
    help = "
        if present, remove outfiles in the outdirectory prior to generating
        outfiles <logical>
    "
)


#  Parse the arguments --------------------------------------------------------
test_in_RStudio <- TRUE  # Hardcode T for testing in RStudio; F for CLI
if(isTRUE(test_in_RStudio)) {
    #  RStudio-interactive work
    dir_base <- "."
    dir_in <- paste0(
        dir_base, "/results/kga0/2022-0519_run_join-work"
    )
    dir_out <- paste0(
        dir_base, "/results/kga0/2022-0519_run_make-AS-assignment"
    )
    intersection <- paste0(
        # dir_in, "/", "Disteche_sample_1.dedup.intersection.AS.txt.gz"
        dir_in, "/",
        "Disteche_sample_1.dedup.intersection.AS.no-header.txt.gz"
    )
    strain_1 <- "mm10"
    strain_2 <- "CAST"
    outprefix <- unlist(stringr::str_split(basename(intersection), "\\."))[1]
    threshold <- 0
    chunk <- 1000000L
    remove <- TRUE
    cl <- c(
        "--intersection", intersection,
        "--strain_1", strain_1,
        "--strain_2", strain_2,
        "--outdir", dir_out,
        "--outprefix", outprefix,
        "--threshold", threshold,
        "--chunk", chunk,
        "--remove", remove
    )
    arguments <- parse_args(ap, cl)
    rm(
        ap, cl,
        intersection, strain_1, strain_2,
        dir_base, dir_in, dir_out,
        outprefix, threshold, chunk, remove
    )
} else if(isFALSE(test_in_RStudio)) {
    #  Command-line calls
    arguments <- parse_args(ap)
    rm(ap)
} else {
    stop(paste0(
        "Stopping: Variable 'test_in_RStudio' is not properly set. ",
        "'test_in_RStudio' should be hardcoded as either \"TRUE\" or",
        "\"FALSE\"."
    ))
}
rm(test_in_RStudio)


#  Check on the arguments that were supplied ----------------------------------
stopifnot(file.exists(arguments$intersection))
if(arguments$strain_1 == "") stop("--strain_1 is an empty string.")
if(is.na(arguments$strain_1)) stop("--strain_1 is NA.")
if(arguments$strain_2 == "") stop("--strain_2 is an empty string.")
if(is.na(arguments$strain_2)) stop("--strain_2 is NA.")
if(arguments$outprefix == "") stop("--outprefix is an empty string.")
if(is.na(arguments$outprefix)) stop("--outprefix is NA.")
stopifnot(arguments$threshold %% 1 == 0)
stopifnot(arguments$chunk != 0)
stopifnot(arguments$chunk %% 1 == 0)
stop_if_not_logical(arguments$remove)

#  If it does not exist, then create outfile directory 
dir.create(file.path(arguments$outdir), showWarnings = FALSE)
#TODO Print message if new directory is created

#  Remove already-created outfiles in outdirectory (optional)
if(isTRUE(arguments$remove)) {
    cat(paste0("If present, removing outfiles in the outdirectory...\n"))
    remove_outfiles(
        arguments$outdir,
        arguments$outprefix,
        arguments$strain_1,
        arguments$strain_2
    )
    cat("\n")
}


#  Count the numbers of lines in samples --------------------------------------
#  Tally and record the total number of records in the bam file
cat(paste0(
    "Counting the number of records in ", basename(arguments$intersection),
    "...\n"
))
rec_n <- as.integer(arguments$chunk)
rec_total <- count_records(arguments$intersection)
cat(paste0("Number of records: ", scales::comma(rec_total), "\n"))
cat("\n")

#  Iterate through intersection file in chunks
cat(paste0(
    "Started: Processing ", basename(arguments$intersection), ", writing ",
    "out files based on assignments to \"", arguments$strain_1, "\", \"",
    arguments$strain_2, "\", or \"ambiguous\".\n"
))
n <- ceiling(rec_total / rec_n) %>% as.integer()
bar <- utils::txtProgressBar(min = 0, max = n, initial = 0, style = 3)

denote_outfile <- function(x) {
    paste0(arguments$outdir, "/", arguments$outprefix, ".", x, ".txt.gz")
}

connection <- file(arguments$intersection, "r")
for(i in 1:n) {
    AS <- read.delim(
        connection, nrows = arguments$chunk, header = FALSE
    )
    colnames(AS) <- c("qname", arguments$strain_1, arguments$strain_2)
    AS <- make_assignment(AS, arguments$threshold)
    ambiguous <- AS %>% dplyr::filter(assignment == "ambiguous")
    sample_1 <- AS %>% dplyr::filter(assignment == arguments$strain_1)
    sample_2 <- AS %>% dplyr::filter(assignment == arguments$strain_2)
    
    readr::write_tsv(
        ambiguous[, 1] %>% as.data.frame(),
        file = denote_outfile("ambiguous"),
        append = TRUE
    )
    readr::write_tsv(
        sample_1[, 1] %>% as.data.frame(),
        file = denote_outfile(arguments$strain_1),
        append = TRUE
    )
    readr::write_tsv(
        sample_2[, 1] %>% as.data.frame(),
        file = denote_outfile(arguments$strain_2),
        append = TRUE
    )
    utils::setTxtProgressBar(bar, i)
}
cat(paste0(
    "Completed: Processing ", basename(arguments$intersection), ", writing ",
    "out files based on assignments to \"", arguments$strain_1, "\", \"",
    arguments$strain_2, "\", or \"ambiguous\".\n"
))


#  End the script -------------------------------------------------------------
time_end <- Sys.time()
cat("\n")
cat(paste0(script, " completed: ", time_end, "\n"))
print(convert_time_HMS(time_start, time_end))
cat("\n\n")
rm(time_start, time_end)
