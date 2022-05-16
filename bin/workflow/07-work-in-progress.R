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


initialize_tibble_2 <- function(x) {
    #TODO Description of function
    #
    # :param x: identifier for sample (chr)
    # :return: #TODO
    y <- dplyr::tibble(
        qname = character(),
        !!paste0("AS.", x) := numeric()
    )
    return(y)
}


initialize_tibble_5 <- function(x, y) {
    #TODO Description of function
    #
    # :param x: identifier for sample #1 (chr)
    # :param x: identifier for sample #2 (chr)
    # :return: #TODO
    z <- dplyr::tibble(
        qname = character(),
        !!paste0("AS.", x) := numeric(),
        !!paste0("AS.", y) := numeric(),
        difference = numeric(),
        assignment = character()
    )
    return(z)
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


read_line <- function(x, y, z) {
    #TODO Description of function
    #
    # :param x: txt/txt.gz file, including path (chr)
    # :param y: line number preceding line of interest (int)
    # :param z: identifier for sample (chr)
    # :return: #TODO
    a <- read.delim(
        x,
        header = FALSE,
        skip = y,
        nrows = 1
    ) %>%
        tibble::as_tibble() %>%
        dplyr::rename("qname" = "V1") %>%
        dplyr::rename(!!paste0("AS.", z) := "V2")
    return(a)
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
    # :param x: tibble
    # :param y: outdirectory, including path (chr)
    # :param z: outfile prefix (chr)
    # :param a: sample string (chr)
    readr::write_tsv(
        dplyr::bind_rows(x),
        paste0(y, "/", z, ".assign-", a, ".txt.gz"),
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
    short = "-1",
    arg = "--sample_1",
    type = "character",
    default = NULL,
    help = "AS.txt.gz file for sample 1, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-2",
    arg = "--sample_2",
    type = "character",
    default = NULL,
    help = "AS.txt.gz file for sample 2, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-3",
    arg = "--strain_1",
    type = "character",
    default = NULL,
    help = "strain name for sample 1 <chr>"
)
ap <- add_argument(
    ap,
    short = "-4",
    arg = "--strain_2",
    type = "character",
    default = NULL,
    help = "strain name for sample 2 <chr>"
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
        threshold alignment score; the absolute value of the difference in
        alignment scores between samples 1 and 2 must be greater than this
        value in order for a strain-specific assignment to be made; if not, the
        assignment will be \"ambiguous\" <positive int>
    "
)
ap <- add_argument(
    ap,
    short = "-b",
    arg = "--objects",
    type = "logical",
    default = FALSE,
    help = "
        store assignments in objects; WARNING: when using large infiles,
        setting this option to TRUE will result in large objects being stored
        in memory, which is likely to cause memory issues <logical>
    "
)
ap <- add_argument(
    ap,
    short = "-l",
    arg = "--lists",
    type = "logical",
    default = FALSE,
    help = "write out assignments as lists <logical>"
)
ap <- add_argument(
    ap,
    short = "-m",
    arg = "--max_alone",
    type = "integer",
    default = 10000,
    help = "
        if --lists is TRUE, then this parameter sets the maximum number of rows
        in the 'alone_1' and 'alone_2' objects; once 'alone_#' increases beyond
        this size, then the first n rows of 'alone_#' (--rows_alone) will be
        written to a txt.gz file and then removed from the 'alone_#' object
        <int>
    "
)
ap <- add_argument(
    ap,
    short = "-n",
    arg = "--rows_alone",
    type = "integer",
    default = 100,
    help = "
        if --lists is TRUE, then this parameter sets the first n rows of
        'alone_#' to be written to a txt.gz file upon surpassing the maximum
        number of rows allowed in 'alone_#' (--max_alone); those rows are then
        removed from the 'alone_#' object <int>
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


#  Parse the arguments --------------------------------------------------------
test_in_RStudio <- TRUE  # Hardcode T for testing in RStudio; F for CLI
if(isTRUE(test_in_RStudio)) {
    #  RStudio-interactive work
    dir_base <- "."
    dir_data <- "data"
    dir_in <- paste0(dir_base, "/", dir_data, "/", "files_bam")
    dir_out <- paste0(dir_base, "/", dir_data, "/", "files_bam")
    sample_1 <- paste0(
        dir_in, "/", "Disteche_sample_1.dedup.mm10.corrected.mm10.AS.txt.gz"
    )
    sample_2 <- paste0(
        dir_in, "/", "Disteche_sample_1.dedup.CAST.corrected.CAST.AS.txt.gz"
    )
    strain_1 <- "mm10"
    strain_2 <- "CAST"
    outprefix <- unlist(stringr::str_split(basename(sample_1), "\\."))[1]
    threshold <- 0
    objects <- FALSE
    lists <- TRUE
    max_alone <- 1000
    max_rows <- 100
    remove <- TRUE
    cl <- c(
        "--sample_1", sample_1,
        "--sample_2", sample_2,
        "--strain_1", strain_1,
        "--strain_2", strain_2,
        "--outdir", dir_out,
        "--outprefix", outprefix,
        "--threshold", threshold,
        "--objects", objects,
        "--lists", lists,
        "--max_alone", max_alone,
        "--rows_alone", max_rows,
        "--remove", remove
    )
    arguments <- parse_args(ap, cl)
    rm(
        ap, cl,
        sample_1, sample_2, strain_1, strain_2,
        dir_base, dir_data, dir_in, dir_out,
        outprefix, threshold, objects, lists, max_alone, max_rows, remove
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
stopifnot(file.exists(arguments$sample_1))
stopifnot(file.exists(arguments$sample_2))
if(arguments$strain_1 == "") stop("--strain_1 is an empty string.")
if(is.na(arguments$strain_1)) stop("--strain_1 is NA.")
if(arguments$strain_2 == "") stop("--strain_2 is an empty string.")
if(is.na(arguments$strain_2)) stop("--strain_2 is NA.")
if(arguments$outprefix == "") stop("--outprefix is an empty string.")
if(is.na(arguments$outprefix)) stop("--outprefix is NA.")
stopifnot(arguments$threshold %% 1 == 0)
stop_if_not_logical(arguments$objects)
stop_if_not_logical(arguments$lists)
stopifnot(arguments$max_alone %% 1 == 0)
stopifnot(arguments$max_rows %% 1 == 0)
stop_if_not_logical(arguments$remove)

#  If it does not exist, then create outfile directory 
dir.create(file.path(arguments$outdir), showWarnings = FALSE)
#TODO Print message if new directory is created


#  Initialize tibbles for storing values --------------------------------------
sample_1 <- sample_2 <- ambiguous <- AS <- salvage_AS <- initialize_tibble_5(
    arguments$strain_1,
    arguments$strain_2
)
line_1 <- alone_1 <- salvage_1 <- initialize_tibble_2(arguments$strain_1)
line_2 <- alone_2 <- salvage_2 <- initialize_tibble_2(arguments$strain_2)


#  Count the numbers of lines in samples --------------------------------------
#  Tally, record the total number of records in the files for samples #1 and #2
cat(paste0(
    "Counting the numbers of records in ", basename(arguments$sample_1),
    " and ", basename(arguments$sample_2), "...\n"
))

count.sample_1 <- count_records(arguments$sample_1)
count.sample_2 <- count_records(arguments$sample_2)

cat(paste0(
    "- ", basename(arguments$sample_1), ": ",
    scales::comma(count.sample_1), " lines\n"
))
cat(paste0(
    "- ", basename(arguments$sample_2), ": ",
    scales::comma(count.sample_2), " lines\n"
))
cat("\n")

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


#  Process lines of interest --------------------------------------------------
n <- pmax(count.sample_1, count.sample_2)
# n <- 10000L  # For tests only

bar <- utils::txtProgressBar(min = 0, max = n, initial = 0, style = 3)
for(i in 1:n) {
    # #####  For tests only  #####
    # i <- 2
    
    #  Read in line of txt.gz file only if that line number is present in the
    #+ txt.gz file; i.e., don't read in a line number greater than the total
    #+ number of lines in the file
    if(i < count.sample_1) {
        line_1 <- read_line(arguments$sample_1, i, arguments$strain_1)
    }
    if(i < count.sample_2) {
        line_2 <- read_line(arguments$sample_2, i, arguments$strain_2)
    }
    
    #  Test: Do QNAMEs match on a per-line basis? If they match, then perform
    #+ logic for AS comparisons (see function make_assignment()); if not, then
    #+ store values in "alone_" variables (further description below)
    if(line_1[1] == line_2[1]) {
        AS <- dplyr::bind_cols(line_1, line_2[2]) %>%
            make_assignment(., arguments$threshold)
        
        if(AS[5] == arguments$strain_1) {
            if(arguments$objects == TRUE) {
                sample_1 <- dplyr::bind_rows(sample_1, AS)
            }
            if(arguments$lists == TRUE) {
                write_list(
                    AS[1],
                    arguments$outdir,
                    arguments$outprefix,
                    arguments$strain_1
                )
            }
        } else if(AS[5] == arguments$strain_2) {
            if(arguments$objects == TRUE) {
                sample_2 <- dplyr::bind_rows(sample_2, AS)
            }
            if(arguments$lists == TRUE) {
                write_list(
                    AS[1],
                    arguments$outdir,
                    arguments$outprefix,
                    arguments$strain_2
                )
            }
        } else if(AS[5] == "ambiguous") {
            if(arguments$objects == TRUE) {
                ambiguous <- dplyr::bind_rows(ambiguous, AS)
            }
            if(arguments$lists == TRUE) {
                write_list(
                    AS[1],
                    arguments$outdir,
                    arguments$outprefix,
                    "ambiguous"
                )
            }
        }
    } else {
        if(i < count.sample_1) alone_1 <- dplyr::bind_rows(alone_1, line_1)
        if(i < count.sample_2) alone_2 <- dplyr::bind_rows(alone_2, line_2)
        
            # #####  For tests only  #####
            # alone_1 <- alone_1 %>%
            #     dplyr::bind_rows(tibble::tibble(
            #         qname = "AACCAATAACAAGTCAACGGACGCTTCGTCCGAGATAAGA:2423713332",
            #         AS.mm10 = -10
            #     ))
            # alone_2

        #  Attempt to "salvage" unmatched entries
        if(sum(alone_1$qname %in% alone_2$qname) == 1) {
            salvage_1 <- alone_1[alone_1$qname %in% alone_2$qname, ]
            salvage_2 <- alone_2[alone_2$qname %in% alone_1$qname, ]
            salvage_AS <- dplyr::bind_cols(salvage_1, salvage_2[2]) %>%
                make_assignment(., arguments$threshold)
            
            if(salvage_AS[5] == arguments$strain_1) {
                if(arguments$objects == TRUE) {
                    sample_1 <- dplyr::bind_rows(sample_1, salvage_AS)
                }
                if(arguments$lists == TRUE) {
                    write_list(
                        salvage_AS[1],
                        arguments$outdir,
                        arguments$outprefix,
                        arguments$strain_1
                    )
                }
            } else if(salvage_AS[5] == arguments$strain_2) {
                if(arguments$objects == TRUE) {
                    sample_2 <- dplyr::bind_rows(sample_2, salvage_AS)
                }
                if(arguments$lists == TRUE) {
                    write_list(
                        salvage_AS[1],
                        arguments$outdir,
                        arguments$outprefix,
                        arguments$strain_2
                    )
                }
            } else if(salvage_AS[5] == "ambiguous") {
                if(arguments$objects == TRUE) {
                    ambiguous <- dplyr::bind_rows(ambiguous, salvage_AS)
                }
                if(arguments$lists == TRUE) {
                    write_list(
                        salvage_AS[1],
                        arguments$outdir,
                        arguments$outprefix,
                        "ambiguous"
                    )
                }
            }
            
            #  To prevent duplicate searches, and since they're no longer
            #+ "alone", remove the salvaged QNAME from alone_1 and alone_2
            alone_1 <- alone_1[!(alone_1$qname %in% salvage_1$qname), ]
            alone_2 <- alone_2[!(alone_2$qname %in% salvage_2$qname), ]
            #QUESTION Does this address the below (potential) issue?
            
            #  An alternative to the above subsetting approach (this is the 
            #+ initial approach I used):
            #+ 
            #+ The dplyr::filter approach is finding all rows not equal to the
            #+ one to exclude, then rewriting those to the tibble; I think this
            #+ will slow down the script as alone_1 and alone_2 increase in
            #+ size
            # alone_1 <- alone_1 %>%
            #     dplyr::filter(qname != as.character(salvage_1[1]))
            # alone_2 <- alone_2 %>%
            #     dplyr::filter(qname != as.character(salvage_2[1]))
            
            # n <- 1000L
            # subset approach
            # [1] "00:00:26"
            # 
            # dplyr::filter approach
            # [1] "00:00:27"
            # 
            # n <- 10000L
            # subset approach
            # [1] "00:02:28"
            # 
            # dplyr::filter approach
            # [1] "00:02:44"
        }
        
        #  Only write out the lists of "alone" QNAMEs (i.e., in one sample but
        #+ not other) after having read in and run logic on all lines; #TODO
        #+ May need to consider a strategy in which we write x number of
        #+ entries and remove them from memory after reaching a certain number
        #+ of lines in the object or after reaching a certain size in RAM;
        #+ e.g., if the "alone" object reaches 1GB in size, then write the
        #+ first x number of lines to disc and remove them from memory (I think
        #+ this approach should be fine because of the lexical-sorted nature
        #+ of the input: earlier entries in one sample are very likely to be
        #+ with earlier entries of other sample and very unlikely to be with
        #+ later entries in that sample)
        if(nrow(alone_1) > arguments$max_alone) {  #NOTE Strategy to deal with growing object size
            write_list(  #TODO Turn this code unit into a function
                alone_1[1:arguments$rows_alone, ][1],
                arguments$outdir,
                arguments$outprefix,
                paste0(arguments$strain_1, "-only")  #TODO Make an argument to write QNAME entries unique to a given sample to a separate list
                # arguments$strain_1
            )
            alone_1 <- tail(alone_1, -arguments$rows_alone)
        }
        if(nrow(alone_2) > arguments$max_alone) {  #NOTE Strategy to deal with growing object size
            write_list(  #TODO Turn this code unit into a function
                alone_2[1:arguments$rows_alone, ][1],
                arguments$outdir,
                arguments$outprefix,
                paste0(arguments$strain_2, "-only")  #TODO Make an argument to write QNAME entries unique to a given sample to a separate list
                # arguments$strain_2
            )
            alone_2 <- tail(alone_2, -arguments$rows_alone)
        }
        if(i == n) {
            if(arguments$lists == TRUE) {
                write_list(
                    alone_1[1],
                    arguments$outdir,
                    arguments$outprefix,
                    paste0(arguments$strain_1, "-only")  #TODO Make an argument to write QNAME entries unique to a given sample to a separate list
                    # arguments$strain_1
                )
                write_list(
                    alone_2[1],
                    arguments$outdir,
                    arguments$outprefix,
                    paste0(arguments$strain_2, "-only")  #TODO Make an argument to write QNAME entries unique to a given sample to a separate list
                    # arguments$strain_2
                )
            }
        }
    }
    utils::setTxtProgressBar(bar, i)
}
# for i in *.assign-*.txt.gz; do decompress_gzip "${i}"; done


#  End the script -------------------------------------------------------------
time_end <- Sys.time()
cat("\n")
cat(paste0(script, " completed: ", time_end, "\n"))
print(convert_time_HMS(time_start, time_end))
cat("\n\n")
rm(time_start, time_end)