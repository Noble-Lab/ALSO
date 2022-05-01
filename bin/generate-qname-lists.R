#!/usr/bin/env Rscript

#  generate-qname-lists.R
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
#     # #TODO Description of function
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


countLinesOutfile <- function(x, y, z) {
    # #TODO Description of function
    # 
    # :param x: set one of five options to count lines in outfile: m (mated),
    #           u (unmated), a (ambiguous), t (trans), d (duplicated)
    # :param y: outfile directory, including path
    # :param z: bam infile without path
    if(x == "m") status <- "mated"
    if(x == "u") status <- "unmated"
    if(x == "a") status <- "ambiguous"
    if(x == "t") status <- "trans"
    if(x == "d") status <- "duplicated"
    if(x == "s") status <- "singleton"
    #TODO Check if not one of the above options
    outdir <- y
    bam <- z

    if(file.exists(paste0(
        outdir, "/", gsub(".bam", paste0(".", status, ".txt.gz"), bam)
    ))) {
        lines <- system(
            paste0(
                "gunzip -c ", outdir, "/",
                gsub(".bam", paste0(".", status, ".txt.gz"), bam),
                " | wc -l"
            ),
            intern = TRUE
        )
        return(print(paste0("Lines in ", status, ".txt.gz: ", lines)))
    }
}


countRecords <- function(x) {
    # Count number of records in bam file
    # 
    # :param x: bam file, including path (chr)
    # :return y: number of records in bam file (int)
    y <- system(paste0("samtools view -c ", x), intern = TRUE) %>%
        as.integer()
    return(y)
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


loadFields <- function(x, y) {
    # #TODO Description of function
    # 
    # :param x: 'g', load fields needed for general analysis; 't', load fields
    #           needed for analysis of interchromosomal paired reads
    # :param y: object of class BamFile
    # :return z: tibble made up of qname, rname, mrnm, groupid, and
    #            mate_status fields
    #TODO Check that BamFile asMates is TRUE
    #TODO Check that BamFile isOpen is TRUE
    if(x == "g") vars <- c("qname", "groupid", "mate_status")
    if(x == "t") vars <- c("qname", "rname", "mrnm", "groupid", "mate_status")
    #TODO Check if not one of the above options

    z <- Rsamtools::scanBam(y, param = ScanBamParam(what = vars))
    z <- Map(as.data.frame, z) %>%
        as.data.frame() %>%
        tibble::as_tibble(column_name = vars)
    return(z)
}


removeOutfiles <- function(x, y) {
    # Remove any outfiles present in directory
    # 
    # :param x: directory from which to remove outfiles
    # :param y: bam file from which outfiles were derived
    system(paste0(
        "rm -f -- ",
        x, "/", gsub(".bam", ".mated.txt.gz", y), " ",
        x, "/", gsub(".bam", ".mated-tally.txt.gz", y), " ",
        x, "/", gsub(".bam", ".mated.txt", y), " ",
        x, "/", gsub(".bam", ".mated-tally.txt", y), " ",
        x, "/", gsub(".bam", ".unmated.txt.gz", y), " ",
        x, "/", gsub(".bam", ".unmated-tally.txt.gz", y), " ",
        x, "/", gsub(".bam", ".unmated.txt", y), " ",
        x, "/", gsub(".bam", ".unmated-tally.txt", y), " ",
        x, "/", gsub(".bam", ".ambiguous.txt.gz", y), " ",
        x, "/", gsub(".bam", ".ambiguous-tally.txt.gz", y), " ",
        x, "/", gsub(".bam", ".ambiguous.txt", y), " ",
        x, "/", gsub(".bam", ".ambiguous-tally.txt", y), " ",
        x, "/", gsub(".bam", ".duplicated.txt.gz", y), " ",
        x, "/", gsub(".bam", ".duplicated-tally.txt.gz", y), " ",
        x, "/", gsub(".bam", ".duplicated.txt", y), " ",
        x, "/", gsub(".bam", ".duplicated-tally.txt", y), " ",
        x, "/", gsub(".bam", ".singleton.txt.gz", y), " ",
        x, "/", gsub(".bam", ".singleton-tally.txt.gz", y), " ",
        x, "/", gsub(".bam", ".singleton.txt", y), " ",
        x, "/", gsub(".bam", ".singleton-tally.txt", y), " ",
        x, "/", gsub(".bam", ".trans.txt.gz", y), " ",
        x, "/", gsub(".bam", ".trans-tally.txt.gz", y), " ",
        x, "/", gsub(".bam", ".trans.txt", y), " ",
        x, "/", gsub(".bam", ".trans-tally.txt", y)
    ))
}


stopIfNotLogical <- function(x) {
    # #TODO Description of function
    # 
    # :param x: single-value logical vector to evaluate (logical)
    stopifnot(x == TRUE | x == FALSE)
}


writeDuplicatedQnames <- function(x, y, z) {
    # Write out txt.gz outfile(s) for duplicate QNAMEs (>2) in bam file
    #
    # :param x: tibble containing the variable 'qname' (tbl_df)
    # :param y: for list txt.gz outfile, write out only unique QNAMEs (logical)
    # :param z: in addition to list txt.gz outfile, write out table of unique
    #           duplicate QNAME entries and associated tallies (logical)
    # :return: write out list of either all or unique duplicated QNAME entries;
    #          optionally, write out table of unique duplicated QNAME entries
    #          and associated tallies
    #TODO Check class(x)
    #TODO Check variables in x
    #TODO Check that y is either TRUE or FALSE, and nothing else
    #TODO Check that z is either TRUE or FALSE, and nothing else
    a <- x$qname %>%
        table() %>%
        dplyr::as_tibble() %>%
        dplyr::rename("qname" = ".") %>%
        dplyr::filter(n > 2)

    if(isFALSE(y)) {
        b <- x$qname[x$qname %in% a$qname]
        if(length(b) != 0) {
            readr::write_tsv(
                as.data.frame(b),
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
    } else {
        if(nrow(a) != 0) {
            readr::write_tsv(
                as.data.frame(a$qname),
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

    #TODO Function to write out tally
    if(isTRUE(z)) {
        if(nrow(a) != 0) {
            readr::write_tsv(
                a,
                paste0(
                    arguments$outdir, "/",
                    gsub(
                        ".bam",
                        paste0(".duplicated-tally.txt.gz"),
                        basename(arguments$bam)
                    )
                ),
                append = TRUE
            )
        }
    }
}


writeMatedQnames <- function(x, y, z, a) {
    # Evaluate mate status for reads in a bam file, reporting if reads with
    # status "mated", "unmated", or "ambiguous" are present; optionally, write
    # out txt.gz tables comprised of 'qname', 'groupid', and 'mate_status' for
    # a given status of "mated", "unmated", or "ambiguous"
    # 
    # :param x: A 3 x 2 tibble or dataframe: the first column is $mate_status
    #           (factor) with levels "mated", "unmated", and "ambiguous"; the
    #           second column is $n with counts (int) for each level
    # :param y: status to test: "mated", "unmated", or "ambiguous" (chr)
    # :param z: for list txt.gz outfile, write out only unique QNAMEs (logical)
    # :param a: in addition to list txt.gz outfile, write out table of unique
    #           duplicate QNAME entries and associated tallies (logical)
    # :return: write out list of either all or unique QNAME entries for
    #          $mate_status of "mated", "unmated", or "ambiguous"; optionally,
    #          write out table of unique duplicated QNAME entries and
    #          associated tallies for $mate_status of "mated", "unmated", or
    #          "ambiguous"
    #TODO Check for class(x)...
    if (class(y) != "character") {
        stop("Exiting: Parameter 'y' should be class 'character'.")
    }
    #TODO Check(s) for z
    #TODO Check(s) for a
    
    if(x[x$mate_status == y, 2] != 0) {
        qnames <- pertinent[pertinent$mate_status == y, ]$qname
        readr::write_tsv(
            if(isFALSE(z)) {
                as.data.frame(unique(qnames))
            } else {
                as.data.frame(qnames)
            },
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
        
        #TODO Function to write out tally
        if(isTRUE(a)) {
            qnames <- qnames %>%
                table() %>%
                dplyr::as_tibble() %>%
                dplyr::rename("qname" = ".")
            readr::write_tsv(
                qnames,
                paste0(
                    arguments$outdir, "/",
                    gsub(
                        ".bam",
                        paste0(".", y, "-tally.txt.gz"),
                        basename(arguments$bam)
                    )
                ),
                append = TRUE
            )
        }
    }
}


writeSingletonQnames <- function(x, y, z) {
    # Write out txt.gz outfile(s) for singleton QNAMEs (1) in bam file
    #
    # :param x: tibble containing the variable 'qname' (tbl_df)
    # :param y: for list txt.gz outfile, write out only unique QNAMEs (logical)
    # :param z: in addition to list txt.gz outfile, write out table of unique
    #           singleton QNAME entries and associated tallies (logical)
    # :return: write out list of either all or unique singleton QNAME entries;
    #          optionally, write out table of unique singleton QNAME entries
    #          and associated tallies
    #TODO Check class(x)
    #TODO Check variables in x
    a <- x$qname %>%
        table() %>%
        dplyr::as_tibble() %>%
        dplyr::rename("qname" = ".") %>%
        dplyr::filter(n == 1)

    if(nrow(a) != 0) {
        if(isFALSE(y)) {
            b <- x$qname[x$qname %in% a$qname]
            if(length(b) != 0) {
                readr::write_tsv(
                    as.data.frame(b),
                    paste0(
                        arguments$outdir, "/",
                        gsub(
                            ".bam",
                            paste0(".singleton.txt.gz"),
                            basename(arguments$bam)
                        )
                    ),
                    append = TRUE
                )
            }
        } else {
            readr::write_tsv(
                as.data.frame(a$qname),
                paste0(
                    arguments$outdir, "/",
                    gsub(
                        ".bam",
                        paste0(".singleton.txt.gz"),
                        basename(arguments$bam)
                    )
                ),
                append = TRUE
            )
        }
            
        #TODO Function to write out tally
        if(isTRUE(z)) {
            readr::write_tsv(
                a,
                paste0(
                    arguments$outdir, "/",
                    gsub(
                        ".bam",
                        paste0(".singleton-tally.txt.gz"),
                        basename(arguments$bam)
                    )
                ),
                append = TRUE
            )
        }
    }
}


writeTransMates <- function(x, y, z) {
    # Write out txt.gz outfile(s) for interchromosomal QNAMEs in bam file
    # 
    # :param x: tibble containing variables 'qname', 'rname', and 'mrnm'
    #           (tbl_df)
    # :param y: for list txt.gz outfile, write out only unique QNAMEs (logical)
    # :param z: in addition to list txt.gz outfile, write out table of unique
    #           interchromosomal QNAME entries and associated tallies (logical)
    # :return: write out list of either all or unique trans QNAME entries;
    #          optionally, write out table of unique trans QNAME entries
    #          and associated tallies
    #TODO Check class(x)
    #TODO Check variables in x
    a <- x[x$rname != x$mrnm, ]
    
    if(nrow(a) != 0) {
        a <- a$qname %>%
            table() %>%
            dplyr::as_tibble() %>%
            dplyr::rename("qname" = ".")
        
        if(isFALSE(y)) {
            a <- x$qname[x$qname %in% a$qname]
            if(length(a) != 0) {
                readr::write_tsv(
                    as.data.frame(a),
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
        } else {
            readr::write_tsv(
                # as.data.frame(unique(y$qname)),
                as.data.frame(y$qname),
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
            
        #TODO Function to write out tally
        if(isTRUE(z)) {
            readr::write_tsv(
                a,
                paste0(
                    arguments$outdir, "/",
                    gsub(
                        ".bam",
                        paste0(".trans-tally.txt.gz"),
                        basename(arguments$bam)
                    )
                ),
                append = TRUE
            )
        }
    }
}


#  Source libraries, adjust settings ------------------------------------------
libraries <- c("argparser", "Rsamtools", "scales", "tidyverse")
for(i in 1:length(libraries)) checkLibraryAvailability(libraries[i])
importLibrary(libraries)
rm(i, libraries)

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)


#  Parse arguments ------------------------------------------------------------
script <- "generate-qname-lists.R"

#  Create a parser
ap <- arg_parser(
    name = script,
    description = "
        Script outputs QNAME txt.gz list(s) to be used when filtering bam files
        to include or exclude reads with specific QNAMEs; the user has options
        to output 'mated' 'unmated', 'ambiguous', 'duplicated', 'trans'
        (interchromosomal), and 'singleton' QNAME lists; in addition to lists
        of QNAMEs, the user has an option to output lists of unique QNAMEs with
        associated tallies as well.
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
    short = "-o",
    arg = "--outdir",
    type = "character",
    default = NULL,
    help = "directory for saving outfile(s), including path <chr>"
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
        save duplicated QNAME (i.e., QNAME entries > 2) list in a txt.gz
        outfile <logical>
    "
)
ap <- add_argument(
    ap,
    short = "-s",
    arg = "--singleton",
    type = "logical",
    default = TRUE,
    help = "
        save singleton QNAME (i.e., QNAME entries = 1) list in a txt.gz outfile
        <logical>
    "
)
ap <- add_argument(
    ap,
    short = "-q",
    arg = "--unique",
    type = "logical",
    default = FALSE,
    help = "
        for txt.gz outfile(s), save unique QNAME entries, not all QNAME entries
        <logical>
    "
)
ap <- add_argument(
    ap,
    short = "-l",
    arg = "--tally",
    type = "logical",
    default = TRUE,
    help = "
        in addition to list(s) of QNAMEs, generate list(s) with tallies in
        txt.gz outfile(s) <logical>
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


#  Parse the command line arguments -------------------------------------------
test_in_RStudio <- FALSE  # Hardcode T for testing in RStudio; F for CLI
if(isTRUE(test_in_RStudio)) {
    #  RStudio-interactive work
    dir_base <- "/Users/kalavattam/Dropbox/UW/projects-etc"
    dir_proj <- paste0(dir_base, "/", "2021_kga0_4dn-mouse-cross")
    # dir_data <- "results/kga0/2022-0416-0418_test-preprocessing-module"
    dir_data <- "data/files_bam"
    dir_in_out <- paste0(dir_proj, "/", dir_data)
    bam <- "Disteche_sample_13.dedup.mm10.sort-c.rm.bam"
    bai <- "Disteche_sample_13.dedup.mm10.sort-c.rm.bam.bai"
    chunk <- 100000
    mated <- FALSE
    unmated <- TRUE
    ambiguous <- TRUE
    trans <- TRUE
    duplicated <- TRUE
    singleton <- TRUE
    uniq <- TRUE
    tally <- TRUE
    remove <- TRUE
    cl <- c(
        #  Arguments for analysis
        "--bam", paste0(dir_in_out, "/", bam),
        "--bai", paste0(dir_in_out, "/", bai),
        "--outdir", dir_in_out,
        "--chunk", chunk,
        "--mated", mated,
        "--unmated", unmated,
        "--ambiguous", ambiguous,
        "--trans", trans,
        "--duplicated", duplicated,
        "--singleton", singleton,
        "--unique", uniq,
        "--tally", tally,
        "--remove", remove
    )
    arguments <- parse_args(ap, cl)  
    rm(
        ap, cl,
        dir_base, dir_proj, dir_data, dir_in_out,
        bam, bai, chunk, mated, unmated, ambiguous,
        duplicated, trans, singleton,
        uniq, tally, remove
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
stopIfNotLogical(arguments$mated)
stopIfNotLogical(arguments$unmated)
stopIfNotLogical(arguments$ambiguous)
stopIfNotLogical(arguments$trans)
stopIfNotLogical(arguments$duplicated)
stopIfNotLogical(arguments$singleton)
stopIfNotLogical(arguments$unique)
stopIfNotLogical(arguments$tally)
stopIfNotLogical(arguments$remove)

#  If outfile directory does not exist, then create it
dir.create(file.path(arguments$outdir), showWarnings = FALSE)
#TODO Print message if new directory is created


#  Set up variables, environment prior to loading in .bam information... ------
#+ ...including mate information
print(paste0(
    "Started: Using Rsamtools to load in '", basename(arguments$bam),
    "' and reading various fields such as 'qname' into memory (in chunks of ",
    scales::comma(arguments$chunk), " records per iteration)."
))
cat("\n")

#  Record the number of records in a chunk; tally and record the total number
#+ of records in the bam file
print(paste0(
    "Counting the number of records in ", basename(arguments$bam), "..."
))
rec_n <- as.integer(arguments$chunk)
rec_total <- countRecords(arguments$bam)
print(paste0("Number of records: ", scales::comma(rec_total)))
cat("\n")

#  Remove already-created outfiles in outdirectory (optional)
if(isTRUE(arguments$remove)) {
    print(paste0("If present, removing outfiles in the outdirectory..."))
    removeOutfiles(arguments$outdir, basename(arguments$bam))
    cat("\n")
}


#  Load in, "open" bam file in user-defined chunks; evaluate via while loop ---
bam <- Rsamtools::BamFile(arguments$bam, index = arguments$bai, asMates = TRUE)
Rsamtools::yieldSize(bam) <- arguments$chunk
open(bam)

if(isTRUE(arguments$tally)) tally <- TRUE else tally <- FALSE
if(isTRUE(arguments$unique)) uniq <- TRUE else uniq <- FALSE

time_start <- Sys.time()
cat((2 * rec_n), " ")
#  Iterate through bam file in chunks
while(rec_n < (rec_total / 2) + arguments$chunk) {
    #  Load in pertinent bam fields
    if(isTRUE(arguments$trans)) {
        pertinent <- loadFields("t", bam)
    } else {
        pertinent <- loadFields("g", bam)
    }
    
    #  Determine and evaluate mate status levels present in the data; write out
    #+ tables for the mate statuses
    status <- pertinent$mate_status %>%
        table() %>%
        dplyr::as_tibble() %>%
        dplyr::rename("mate_status" = ".")

    if(isTRUE(arguments$mated)) {
        writeMatedQnames(status, "mated", uniq, tally)
    }
    if(isTRUE(arguments$unmated)) {
        writeMatedQnames(status, "unmated", uniq, tally)
    }
    if(isTRUE(arguments$ambiguous)) {
        writeMatedQnames(status, "ambiguous", uniq, tally)
    }

    #  Determine and evaluate duplicate QNAME (>2) entries, writing out tables
    if(isTRUE(arguments$duplicated)) {
        writeDuplicatedQnames(pertinent, uniq, tally)
    }

    #  Determine and evaluate singleton QNAME entries, writing out tables
    if(isTRUE(arguments$singleton)) {
        writeSingletonQnames(pertinent, uniq, tally)
    }
    
    #  Determine and evaluate QNAMEs associated with trans reads, writing out
    #+ tables
    if(isTRUE(arguments$trans)) {
        writeTransMates(pertinent, uniq, tally)
    }

    rec_n <- rec_n + arguments$chunk
    cat((2 * rec_n), " ")
}
close(bam)
cat("\n")


#  Count lines in outfiles, report the counts, end the script -----------------
if(isTRUE(arguments$mated)) {
    countLinesOutfile("m", arguments$outdir, basename(arguments$bam))
}
if(isTRUE(arguments$unmated)) {
    countLinesOutfile("u", arguments$outdir, basename(arguments$bam))
}
if(isTRUE(arguments$ambiguous)) {
    countLinesOutfile("a", arguments$outdir, basename(arguments$bam))
}
if(isTRUE(arguments$duplicated)) {
    countLinesOutfile("d", arguments$outdir, basename(arguments$bam))
}
if(isTRUE(arguments$singleton)) {
    countLinesOutfile("s", arguments$outdir, basename(arguments$bam))
}
if(isTRUE(arguments$trans)) {
    countLinesOutfile("t", arguments$outdir, basename(arguments$bam))
}

time_end <- Sys.time()
cat("\n")
print(paste0("Script completed."))
cat("\n\n")
rm(time_start, time_end)
