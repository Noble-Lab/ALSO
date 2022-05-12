#!/usr/bin/env Rscript

#  work-in-progress.R
#  KA


script <- "work-in-progress.R"
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
    short = "-a",
    arg = "--sample_1",
    type = "character",
    default = NULL,
    help = "*.AS.txt.gz file for sample 1, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-c",
    arg = "--sample_2",
    type = "character",
    default = NULL,
    help = "*.AS.txt.gz file for sample 2, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-b",
    arg = "--strain_1",
    type = "character",
    default = NULL,
    help = "strain name for sample 1 <chr>"
)
ap <- add_argument(
    ap,
    short = "-d",
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


#  Parse the arguments --------------------------------------------------------
test_in_RStudio <- TRUE  # Hardcode T for testing in RStudio; F for CLI
if(isTRUE(test_in_RStudio)) {
    #  RStudio-interactive work
    dir_base <- "."
    dir_data <- "data"
    dir_in <- paste0(dir_base, "/", dir_data, "/", "files_bam")
    dir_out <- paste0(dir_base, "/", dir_data, "/", "files_bam")
    sample_1 <- paste0(
        dir_in, "/", "Disteche_sample_1.dedup.mm10.corrected.AS.2600.txt.gz"
    )
    sample_2 <- paste0(
        dir_in, "/", "Disteche_sample_1.dedup.CAST.corrected.AS.2500.txt.gz"
    )
    strain_1 <- "mm10"
    strain_2 <- "CAST"
    threshold <- 0
    outprefix <- unlist(stringr::str_split(basename(sample_1), "\\."))[1]
    cl <- c(
        "--sample_1", sample_1,
        "--sample_2", sample_2,
        "--strain_1", strain_1,
        "--strain_2", strain_2,
        "--outdir", dir_out,
        "--outprefix", outprefix,
        "--threshold", threshold
    )
    arguments <- parse_args(ap, cl)
    rm(
        sample_1, sample_2, strain_1, strain_2,
        dir_base, dir_data, dir_in, dir_out,
        outprefix, threshold
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


#  Count lines in samples, initialize tibbles for storing values --------------
count.sample_1 <- count_records(arguments$sample_1)
count.sample_2 <- count_records(arguments$sample_2)

sample_1 <- sample_2 <- ambiguous <- AS <- salvage_AS <- initialize_tibble_5(
    arguments$strain_1,
    arguments$strain_2
)

line_1 <- alone_1 <- salvage_1 <- initialize_tibble_2(arguments$strain_1)
line_2 <- alone_2 <- salvage_2 <- initialize_tibble_2(arguments$strain_2)


#  Process lines of interest --------------------------------------------------
n <- pmax(count.sample_1, count.sample_2)

# test.line_1 <- read_line(arguments$sample_1, n - 1, arguments$strain_1)
# n <- 1000L
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
            sample_1 <- dplyr::bind_rows(sample_1, AS)
            write_list(
                AS[1],
                arguments$outdir,
                arguments$outprefix,
                arguments$strain_1
            )
        } else if(AS[5] == arguments$strain_2) {
            sample_2 <- dplyr::bind_rows(sample_2, AS)
            write_list(
                AS[1],
                arguments$outdir,
                arguments$outprefix,
                arguments$strain_2
            )
        } else if(AS[5] == "ambiguous") {
            ambiguous <- dplyr::bind_rows(ambiguous, AS)
            write_list(
                AS[1],
                arguments$outdir,
                arguments$outprefix,
                "ambiguous"
            )
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
                sample_1 <- dplyr::bind_rows(sample_1, salvage_AS)
                write_list(
                    salvage_AS[1],
                    arguments$outdir,
                    arguments$outprefix,
                    arguments$strain_1
                )
            } else if(salvage_AS[5] == arguments$strain_2) {
                sample_2 <- dplyr::bind_rows(sample_2, salvage_AS)
                write_list(
                    salvage_AS[1],
                    arguments$outdir,
                    arguments$outprefix,
                    arguments$strain_2
                )
            } else if(salvage_AS[5] == "ambiguous") {
                ambiguous <- dplyr::bind_rows(ambiguous, salvage_AS)
                write_list(
                    salvage_AS[1],
                    arguments$outdir,
                    arguments$outprefix,
                    "ambiguous"
                )
            }
            
            #  To prevent duplicate searches, and since they're no longer
            #+ "alone", remove QNAME from alone_1 and alone_2: The
            #+ dplyr::filter approach is finding all rows not equal to the one
            #+ to exclude, then rewriting those to the tibble; I think this
            #+ will slow down the script as alone_1 and alone_2 increase in
            #+ size
            alone_1 <- alone_1 %>%
                dplyr::filter(qname != as.character(salvage_1[1]))
            alone_2 <- alone_2 %>%
                dplyr::filter(qname != as.character(salvage_2[1]))
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
        if(i == n) {
            write_list(
                alone_1,
                arguments$outdir,
                arguments$outprefix,
                arguments$strain_1
            )
            write_list(
                alone_2,
                arguments$outdir,
                arguments$outprefix,
                arguments$strain_2
            )
        }
    }
    utils::setTxtProgressBar(bar, i)
}
line_1
pryr::object_size(line_1)
line_2
pryr::object_size(line_2)
AS
pryr::object_size(AS)

alone_1
pryr::object_size(alone_1)
alone_2
pryr::object_size(alone_2)

ambiguous
pryr::object_size(ambiguous)
sample_1
pryr::object_size(sample_1)
sample_2
pryr::object_size(sample_2)

i
pryr::object_size(i)


#  End the script -------------------------------------------------------------
time_end <- Sys.time()
cat("\n")
cat(paste0(script, " completed: ", time_end, "\n"))
print(convert_time_HMS(time_start, time_end))
cat("\n\n")
rm(time_start, time_end)


#  Notes on the logic
# do lexicographic sort of ♂.AS.txt.gz, ♀.AS.txt.gz
# 
# 	♀	♂
# 1	1	1
# 2	2	2
# 3	3	4
# 4	4	5
# 5	5	6
# 6	7	7
# 7	9	10
# 8	10	11
# ...
# 
# line 1: match, ∴ do comparison and store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# line 2: match, ∴ do comparison and store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# line 3: mismatch, ∴ store 3 in tmp.♀, store 4 in tmp.♂; compare tmp.♀ with tmp.♂, find no matches; move on to next line
# line 4: mismatch, ∴ store 4 in tmp.♀, store 5 in tmp.♂; compare tmp.♀ with tmp.♂, find that 4♀ and 4♂ match, store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# line 5: mismatch, ∴ store 5 in tmp.♀, store 6 in tmp.♂; compare tmp.♀ with tmp.♂, find that 5♀ and 5♂ match, store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# line 6: match, ∴ do comparison and store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# line 7:	mismatch, ∴ store 9 in tmp.♀, store 10 in tmp.♂; compare tmp.♀ with tmp.♂, find no matches; move on to next line
# line 8: mismatch, ∴ store 10 in tmp.♀, store 11 in tmp.♂; compare tmp.♀ with tmp.♂, find that 10♀ and 10♂ match, store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line


#  pseudocode takes in two files as input

#  file 1 is, e.g., ♂.AS.txt.gz
#+ (header and first line of ♂.AS.txt.gz)
#+ qname   AS
#+ AACCAATAACAAGTCAACGGAATATAGCGGCAAGTCAGCA:1417636611   -5

#  file 2 is, e.g., ♀.AS.txt.gz
#+ (header and first line of ♀.AS.txt.gz)
#+ qname   AS
#+ AACCAATAACAAGTCAACGGAATATAGCGGCAAGTCAGCA:1417636611   0

#  find the max number of records between ♂.AS.txt.gz and ♀.AS.txt.gz
# n = maximum of ♂.AS.txt.gz and ♀.AS.txt.gz

#  loop over n
# for i in n, do
#     read in line i of ♂.AS.txt.gz
#     read in line i of ♀.AS.txt.gz
#
#     test if ♀.AS QNAME matches ♂.AS QNAME
#     	if so, do
#     		AS comparison
#     			if ♀.AS > ♂.AS, append QNAME to final.♀.txt.gz (write to disc)
#     			if ♀.AS < ♂.AS, append QNAME to final.♂.txt.gz (write to disc)
#     			if ♀.AS = ♂.AS, append QNAME to final.⚥.txt.gz (write to disc)
#
#     	if not, do
#     		append QNAME.♀ to object tmp.♀ (store in memory)
#     		append QNAME.♂ to object tmp.♂ (store in memory)
#
#     		test if QNAME in tmp.♀ matches QNAME in tmp.♂
#     			if so, do
#     				AS comparison
#     					if ♀.AS > ♂.AS,
#                            append QNAME to final.♀.txt.gz (write to disc),
#                            then
#                                remove QNAME from both tmp.♀ and tmp.♂ (remove from memory)
#     					if ♀.AS < ♂.AS,
#                            append QNAME to final.♂.txt.gz (write to disc),
#                            then
#                                remove QNAME from both tmp.♀ and tmp.♂  (remove from memory)
#     					if ♀.AS = ♂.AS,
#                            append QNAME to final.⚥.txt.gz (write to disc),
#                            then
#                                remove QNAME from both tmp.♀ and tmp.♂  (remove from memory)
#
#     			if not, do
#     				nothing

#  Use final lists of QNAMEs with picard filterSamReads to filter bam files
#+ such that they contain only maternal, paternal, or ambiguous reads

# How to handle if lines remain in one file but all lines have been read in the other file?  #DONE, see code
# How to handle if, after comparisons between tmp.♀ and tmp.♂, there are no more matches?  #DONE, see code
