#!/usr/bin/env Rscript

#  preprocess-with-R.R
#  KA


#  Functions ------------------------------------------------------------------
checkLibraryAvailability <- function(x) {
    # Check that library is available in environment; return a message if not
    # 
    # :param x: name (string) for library to check (chr)
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


#FIXME
# convertSecondsToHMS <- function(x) {
#     # #TODO Description of function.
#     #
#     # :param x: number of seconds (int)
#     # :return y: #TODO
#     days <- round(x %/% (60 * 60 * 24))
#     hours <- round((x - days * 60 * 60 * 24) %/% (60 * 60))
#     minutes <- round((x - (days * 60 * 60 * 24) - (hours * 60 * 60)) %/% 60)
#     seconds <- round(x - (days * 60 * 60 * 24) - (hours * 60 * 60) - (minutes * 60))
#     str_days <- ifelse(days == 0, "", paste0(days, "d:"))
#     str_hours <- ifelse((hours == 0 & days == 0), "", paste0(hours, "h:"))
#     str_minutes <- ifelse((minutes == 0 & days == 0 & hours == 0), "", paste0(minutes, "m:"))
#     str_seconds <- paste0(seconds, "s")
#     # y <- paste0(str_days, str_hours, str_minutes, str_seconds)
#     y <- paste0(str_hours, str_minutes, str_seconds)
#     return(y)
# }


evaluateDuplicatedQnames <- function(x) {
    # #TODO Description of function.
    #
    # :param x: tibble containing the variable 'qname'
    # :return y: two-column tibble of duplicate 'qname' entries and associated
    #            tallies
    #TODO Check class(x)
    #TODO Check variables in x
    y <- x$qname %>%
        table() %>%
        dplyr::as_tibble() %>%
        dplyr::rename("qname" = ".") %>%
        dplyr::filter(n > 2)
    
    if(nrow(y) != 0) {
        readr::write_tsv(
            as.data.frame(y$qname),
            paste0(
                arguments$outdir, "/",
                gsub(
                    ".bam",
                    paste0(".duplicated.txt.gz"),
                    basename(arguments$bam)
                )
            ),
            append = TRUE
        )
    }
}


evaluateMates <- function(x, y) {
    # Evaluate mate status for reads in a bam file, reporting if reads with
    # status "mated", "unmated", or "ambiguous" are present; optionally, write
    # out txt.gz tables comprised of 'qname', 'groupid', and 'mate_status' for
    # a given status of "mated", "unmated", or "ambiguous"
    # 
    # :param x: A 3 x 2 tibble or dataframe: the first column is $mate_status
    #           (factor) with levels "mated", "unmated", and "ambiguous"; the
    #           second column is $n with counts (int) for each level
    # :param y: status to test: "mated", "unmated", or "ambiguous" (chr)
    #TODO Check for class(x)...
    if (class(y) != "character") {
        stop("Exiting: Parameter 'y' should be class 'character'.")
    }
    
    if(x[x$mate_status == y, 2] != 0) {
        readr::write_tsv(
            as.data.frame(
                #  Print only unique QNAMEs
                unique(pertinent[pertinent$mate_status == y, ]$qname)
            ),
            paste0(
                arguments$outdir, "/",
                gsub(
                    ".bam",
                    paste0(".", y, ".txt.gz"),
                    basename(arguments$bam)
                )
            ),
            append = TRUE
        )
    }
}


evaluateTransMates <- function(x) {
    # #TODO Description of function.
    # 
    # :param x: tibble with variables 'qname', 'rname', and 'mrnm'
    # :return y: single-column tibble of 'qname' entries
    #TODO Check class(x)
    #TODO Check variables in x
    y <- x[x$rname != x$mrnm, ]
    y <- y %>% dplyr::select(-rname, -mrnm)
    
    readr::write_tsv(
        #  Get only unique QNAMEs
        as.data.frame(unique(y$qname)),
        paste0(
            arguments$outdir, "/",
            gsub(
                ".bam",
                paste0(".trans.txt.gz"),
                basename(arguments$bam)
            )
        ),
        append = TRUE
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


loadFieldsMinimum <- function(x) {
    # #TODO Description of function.
    # 
    # :param x: object of class BamFile
    # :return y: tibble made up of qname, groupid, and mate_status fields
    #TODO Check that BamFile asMates is TRUE
    #TODO Check that BamFile isOpen is TRUE
    y <- Rsamtools::scanBam(  
        x,
        param = ScanBamParam(what = c(
            "qname", "groupid", "mate_status"
        ))
    )
    y <- Map(as.data.frame, y) %>%
        as.data.frame() %>%
        tibble::as_tibble(column_name = c(
            "qname", "groupid", "mate_status"
        ))
    return(y)
}


loadFieldsTrans <- function(x) {
    # #TODO Description of function.
    # 
    # :param x: object of class BamFile
    # :return y: tibble made up of qname, rname, mrnm, groupid, and
    #            mate_status fields
    #TODO Check that BamFile asMates is TRUE
    #TODO Check that BamFile isOpen is TRUE
    y <- Rsamtools::scanBam(
        x,
        param = ScanBamParam(what = c(
            "qname", "rname", "mrnm", "groupid", "mate_status"
        ))
    )
    y <- Map(as.data.frame, y) %>%
        as.data.frame() %>%
        tibble::as_tibble(column_name = c(
            "qname", "rname", "mrnm", "groupid", "mate_status"
        ))
    return(y)
}


#  Source libraries, adjust settings ------------------------------------------
libraries <- c("argparser", "pryr", "Rsamtools", "scales", "tidyverse")
for(i in 1:length(libraries)) checkLibraryAvailability(libraries[i])
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
        view -hN' or filter-qname.py; the user has the options to output txt.gz
        lists for 'mated' 'unmated' and 'ambiguous' reads; outfile namess are
        derived from the name of the bam infile.
        
        #TODO Handle duplicate QNAME issue in which a very small number of
        QNAME entries are >2; filter those out in this script or prior to
        reading in the bam file? #ANSWER Prior to reading in the bam file.
        
        #TODO See also https://rdrr.io/bioc/GenomicFiles/man/reduceByYield.html
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
    short = "-a",
    arg = "--bai",
    type = "character",
    default = NULL,
    help = "bam index, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-c",
    arg = "--chunk",
    type = "integer",
    default = 100000,
    help = "number of records to read into memory at a single time <even int>"
)
ap <- add_argument(
    ap,
    short = "-m",
    arg = "--mated",
    type = "logical",
    default = FALSE,
    help = "save mated read QNAME list in a txt.gz outfile <logical>"
)
ap <- add_argument(
    ap,
    short = "-u",
    arg = "--unmated",
    type = "logical",
    default = TRUE,
    help = "save unmated read QNAME list in a txt.gz outfile <logical>"
)
ap <- add_argument(
    ap,
    short = "-a",
    arg = "--ambiguous",
    type = "logical",
    default = TRUE,
    help = "save ambiguous read QNAME list in a txt.gz outfile <logical>"
)
ap <- add_argument(
    ap,
    short = "-t",
    arg = "--trans",
    type = "logical",
    default = TRUE,
    help = "save trans read QNAME list in a txt.gz outfile <logical>"
)
ap <- add_argument(
    ap,
    short = "-d",
    arg = "--duplicated",
    type = "logical",
    default = TRUE,
    help = "
        save duplicate QNAME (i.e., QNAME entries > 2) list in a txt.gz outfile
        <logical>
    "
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
    # dir_data <- "results/kga0/2022-0416-0418_test-preprocessing-module"
    dir_data <- "data/files_bam"
    dir_in_out <- paste0(dir_proj, "/", dir_data)
    bam <- "Disteche_sample_13.dedup.CAST.rm.fixmate.bam"
    bai <- "Disteche_sample_13.dedup.CAST.rm.sort-c.bam.bai"
    chunk <- 100000
    mated <- FALSE
    unmated <- TRUE
    ambiguous <- FALSE
    trans <- FALSE
    duplicated <- TRUE
    remove <- TRUE
    cl <- c(
        #  Arguments for analysis
        "--bam", paste0(dir_in_out, "/", bam),
        "--bai", paste0(dir_in_out, "/", bai),
        "--chunk", chunk,
        "--mated", mated,
        "--unmated", unmated,
        "--ambiguous", ambiguous,
        "--trans", trans,
        "--duplicated", duplicated,
        "--outdir", dir_in_out,
        "--remove", remove
    )
    arguments <- parse_args(ap, cl)  
    rm(
        ap, cl,
        dir_base, dir_data, dir_in_out,
        bam, bai, chunk, mated, unmated, ambiguous,
        duplicated, trans,
        remove
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


#  Check that files exist -----------------------------------------------------
stopifnot(file.exists(arguments$bam))
stopifnot(file.exists(arguments$bai))
stopifnot(arguments$chunk != 0)
stopifnot((arguments$chunk %% 2) == 0)
stopifnot(arguments$mated == TRUE | arguments$mated == FALSE)
stopifnot(arguments$unmated == TRUE | arguments$unmated == FALSE)
stopifnot(arguments$ambiguous == TRUE | arguments$ambiguous == FALSE)
stopifnot(arguments$duplicated == TRUE | arguments$duplicated == FALSE)
stopifnot(arguments$trans == TRUE | arguments$trans == FALSE)
stopifnot(arguments$remove == TRUE | arguments$remove == FALSE)


#  If outfile directory does not exist, then create it
dir.create(file.path(arguments$outdir), showWarnings = FALSE)
#TODO Print message if new directory is created


#  Set up variables, environment prior to loading in .bam information... ------
#+ ...including mate information
print(paste0(
    "Started: Using Rsamtools to load in '", basename(arguments$bam),
    "' and '", basename(arguments$bai),
    "', and reading various fields such as 'qname' into memort in chunks of ",
    scales::comma(arguments$chunk), " records per iteration of while loop."
))
cat("\n")

#  Record the number of records in a chunk; tally and record the total number
#+ of records in the bam file
print(paste0(
    "Counting the number of records in ", basename(arguments$bam), "..."
))
rec_n <- arguments$chunk
rec_total <- system(
    paste0("samtools view -c ", arguments$bam), intern = TRUE
) %>%
    as.integer()
print(paste0("Number of records: ", scales::comma(rec_total)))
cat("\n")

#  Remove already-created outfiles in outdirectory (optional)
if(isTRUE(arguments$remove)) {
    print(paste0("If present, removing outfiles in the outdirectory..."))
    system(paste0(
        "rm -f -- ",
        arguments$outdir, "/", gsub(
            ".bam", ".mated.txt.gz", basename(arguments$bam)
        ), " ",
        arguments$outdir, "/", gsub(
            ".bam", ".mated.txt", basename(arguments$bam)), " ",
        arguments$outdir, "/", gsub(
            ".bam", ".unmated.txt.gz", basename(arguments$bam)), " ",
        arguments$outdir, "/", gsub(
            ".bam", ".unmated.txt", basename(arguments$bam)), " ",
        arguments$outdir, "/", gsub(
            ".bam", ".ambiguous.txt.gz", basename(arguments$bam)), " ",
        arguments$outdir, "/", gsub(".bam", ".ambiguous.txt", basename(arguments$bam)), " ",
        arguments$outdir, "/", gsub(".bam", ".duplicated.txt.gz", basename(arguments$bam)), " ",
        arguments$outdir, "/", gsub(".bam", ".duplicated.txt", basename(arguments$bam)), " ",
        arguments$outdir, "/", gsub(".bam", ".trans.txt.gz", basename(arguments$bam)), " ",
        arguments$outdir, "/", gsub(".bam", ".trans.txt", basename(arguments$bam))
    ))
}


#  Load in, "open" bam file in user-defined chunks; evaluate via while loop ---
bam <- Rsamtools::BamFile(arguments$bam, index = arguments$bai, asMates = TRUE)
Rsamtools::yieldSize(bam) <- arguments$chunk
open(bam)

time_start <- Sys.time()
while(rec_n < (rec_total / 2) + arguments$chunk) {
    if(isTRUE(arguments$trans)) {
        pertinent <- loadFieldsTrans(bam)
    } else {
        pertinent <- loadFieldsMinimum(bam)
    }
    
    #  Determine and evaluate mate_status levels present in the data
    mate_status <- pertinent$mate_status %>%
        table() %>%
        dplyr::as_tibble() %>%
        dplyr::rename("mate_status" = ".")

    if(isTRUE(arguments$mated)) evaluateMates(mate_status, "mated")
    if(isTRUE(arguments$unmated)) evaluateMates(mate_status, "unmated")
    if(isTRUE(arguments$ambiguous)) evaluateMates(mate_status, "ambiguous")

    #  Determine and evaluate duplicate QNAME (>2) entries in the data
    if(isTRUE(arguments$duplicated)) evaluateDuplicatedQnames(pertinent)
    
    #  Determine and evaluate trans reads present in the data
    if(isTRUE(arguments$trans)) evaluateTransMates(pertinent)

    rec_n <- rec_n + arguments$chunk
    cat(rec_n, " ")
}
close(bam)

if(isTRUE(arguments$mated)) {
    lines <- system(paste0("gunzip - c", arguments$mated, " | wc -l"))
    print(paste0("*.mated.txt.gz: ", lines))
}
if(isTRUE(arguments$unmated)) {
    lines <- system(paste0("gunzip - c", arguments$unmated, " | wc -l"))
    print(paste0("*.unmated.txt.gz: ", lines))
}
if(isTRUE(arguments$ambiguous)) {
    system(paste0("gunzip - c", arguments$ambiguous, " | wc -l"))
    print(paste0("*.ambiguous.txt.gz: ", lines))
}
if(isTRUE(arguments$duplicated)) {
    system(paste0("gunzip - c", arguments$duplicated, " | wc -l"))
    print(paste0("*.duplicated.txt.gz: ", lines))
}
if(isTRUE(arguments$trans)) {
    system(paste0("gunzip - c", arguments$trans, " | wc -l"))
    print(paste0("*.trans.txt.gz: ", lines))
}

time_end <- Sys.time()
cat("\n")
print(paste0(
    "Completed: Using Rsamtools to load in '", basename(arguments$bam),
    "' and '", basename(arguments$bai),
    "', reading into memory various fields such as 'qname', ",
    "then writing out txt.gz files of non-duplicated QNAMEs.\n\n",
    "Run time: ",
    round(
        as.numeric(unlist(stringr::str_split((time_end - time_start), " "))), 2
    ),
    " minutes"
))
cat("\n")
rm(time_start, time_end)
