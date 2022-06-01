#!/usr/bin/env Rscript

#  generate-qname-lists.R
#  KA


script <- "generate-qname-lists.R"
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


count_lines_outfile <- function(x, y, z) {
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
        return(cat(paste0("Lines in ", status, ".txt.gz: ", lines, "\n")))
    }
}


count_records <- function(x) {
    # Count number of records in bam file
    # 
    # :param x: bam file, including path (chr)
    # :return y: number of records in bam file (int)
    y <- system(paste0("samtools view -c ", x), intern = TRUE) %>%
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


load_fields <- function(x, y) {
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


remove_outfiles <- function(x, y) {
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


stop_if_not_logical <- function(x) {
    # #TODO Description of function
    # 
    # :param x: single-value logical vector to evaluate (logical)
    stopifnot(x == TRUE | x == FALSE)
}


write_duplicated_qnames <- function(x, y, z) {
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

    if(isTRUE(y)) {
        b <- x$qname[x$qname %in% a$qname]
        if(length(b) != 0) {
            readr::write_tsv(
                as.data.frame(unique(b)),  #CHECK Does use of unique() here solve the problem?
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
    } else if(isFALSE(y)) {
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


write_mated_qnames <- function(x, y, z, a) {
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
            if(isTRUE(z)) {
                as.data.frame(unique(qnames))
            } else if(isFALSE(z)) {
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


write_singleton_qnames <- function(x, y, z) {
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
        if(isTRUE(y)) {
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
        } else if(isFALSE(y)) {
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


write_trans_mates <- function(x, y, z) {
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
        
        if(isTRUE(y)) {
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
        } else if(isFALSE(y)) {
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
        Script outputs QNAME txt.gz list(s) to be used when filtering bam files
        to include or exclude reads with specific QNAMEs; the user has options
        to output 'mated' 'unmated', 'ambiguous', 'duplicated', 'trans'
        (interchromosomal), and 'singleton' QNAME lists; in addition to lists
        of QNAMEs, the user has an option to output lists of unique QNAMEs with
        associated tallies.
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
ap <- add_argument(  #FIXME See the comment at line 748
    ap,
    short = "-d",
    arg = "--duplicated",
    type = "logical",
    default = TRUE,
    help = "
        save duplicated QNAME (i.e., QNAME entries > 2) list in a txt.gz
        outfile <logical>

        #BUG A failure occurs with duplicate recognition because of the way
        that the data is read in and analyzed in chunks; duplicates are only
        recognized if they are within the same chunk and, because the data are
        coordinate-sorted, many duplicates are in different chunks; thus, many
        true duplicates are not recognized and persist into outfiles
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
        for txt.gz outfile(s), save unique QNAME entries (if FALSE, save all
        QNAME entries) <logical>
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
    #  Interactive work in RStudio
    dir_base <- "/Users/kalavattam/Dropbox/UW/projects-etc"
    dir_proj <- paste0(dir_base, "/", "2021_kga0_4dn-mouse-cross")
    dir_data <- "data/files_bam"
    dir_in_out <- paste0(dir_proj, "/", dir_data)
    bam <- "Disteche_sample_13.dedup.mm10.sort-c.rm.chr19.bam"
    bai <- "Disteche_sample_13.dedup.mm10.sort-c.rm.chr19.bam.bai"
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
    #  Command line calls
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
stopifnot(arguments$chunk != 0)
stopifnot(arguments$chunk %% 1 == 0)
stopifnot(arguments$chunk %% 2 == 0)
stop_if_not_logical(arguments$mated)
stop_if_not_logical(arguments$unmated)
stop_if_not_logical(arguments$ambiguous)
stop_if_not_logical(arguments$trans)
stop_if_not_logical(arguments$duplicated)
stop_if_not_logical(arguments$singleton)
stop_if_not_logical(arguments$unique)
stop_if_not_logical(arguments$tally)
stop_if_not_logical(arguments$remove)

#  If outfile directory does not exist, then create it
dir.create(file.path(arguments$outdir), showWarnings = FALSE)
#TODO Print message if new directory is created


#  Set up variables, environment prior to loading in .bam information... ------
#+ ...including mate information
cat(paste0(
    "Using Rsamtools to load in '", basename(arguments$bam),
    "' and reading fields such as 'qname' into memory (in chunks of ",
    scales::comma(arguments$chunk), " records per iteration).\n\n"
))

#  Record the number of records in a chunk; tally and record the total number
#+ of records in the bam file
cat(paste0(
    "Counting the number of records in ", basename(arguments$bam), "...\n"
))
rec_n <- as.integer(arguments$chunk)
rec_total <- count_records(arguments$bam)
cat(paste0("Number of records: ", scales::comma(rec_total), "\n\n"))

#  Remove already-created outfiles in outdirectory (optional)
if(isTRUE(arguments$remove)) {
    cat(paste0("If present, removing outfiles in the outdirectory...\n\n"))
    remove_outfiles(arguments$outdir, basename(arguments$bam))
}


#  Load in, "open" bam file in user-defined chunks; evaluate via for loop -----
bam <- Rsamtools::BamFile(arguments$bam, index = arguments$bai, asMates = TRUE)
Rsamtools::yieldSize(bam) <- arguments$chunk
open(bam)

if(isTRUE(arguments$tally)) tally <- TRUE else tally <- FALSE
if(isTRUE(arguments$unique)) uniq <- TRUE else uniq <- FALSE

#  Iterate through bam file in chunks
cat(paste0("Started: Processing ", basename(arguments$bam), "\n"))
n <- ceiling((rec_total / 2) / rec_n) %>% as.integer()
bar <- utils::txtProgressBar(min = 0, max = n, initial = 0, style = 3)

for(i in 1:n) {
# cat((2 * rec_n), " ")
# while(rec_n < (rec_total / 2) + arguments$chunk) {
    #  Load in pertinent bam fields
    if(isTRUE(arguments$trans)) {
        pertinent <- load_fields("t", bam)
    } else {
        pertinent <- load_fields("g", bam)
    }
    
    #  Determine and evaluate mate status levels present in the data; write out
    #+ tables for the mate statuses
    status <- pertinent$mate_status %>%
        table() %>%
        dplyr::as_tibble() %>%
        dplyr::rename("mate_status" = ".")

    if(isTRUE(arguments$mated)) {
        write_mated_qnames(status, "mated", uniq, tally)
    }
    if(isTRUE(arguments$unmated)) {
        write_mated_qnames(status, "unmated", uniq, tally)
    }
    if(isTRUE(arguments$ambiguous)) {
        write_mated_qnames(status, "ambiguous", uniq, tally)
    }

    #  Determine and evaluate duplicate QNAME (>2) entries, writing out tables
    if(isTRUE(arguments$duplicated)) {
        write_duplicated_qnames(pertinent, uniq, tally)
    }
    #FIXME See the note associated with argument 'duplicate' above

    #  Determine and evaluate singleton QNAME entries, writing out tables
    if(isTRUE(arguments$singleton)) {
        write_singleton_qnames(pertinent, uniq, tally)
    }
    
    #  Determine and evaluate QNAMEs associated with trans reads, writing out
    #+ tables
    if(isTRUE(arguments$trans)) {
        write_trans_mates(pertinent, uniq, tally)
    }

    utils::setTxtProgressBar(bar, i)
    # rec_n <- rec_n + arguments$chunk
    # cat((2 * rec_n), " ")
}
close(bam)
rm(i)
cat("\n")
cat(paste0("Completed: Processing ", basename(arguments$bam), "\n\n"))


#  Count lines in outfiles, report the counts, end the script -----------------
if(isTRUE(arguments$mated)) {  #FIXME These calls sometimes throw errors
    count_lines_outfile("m", arguments$outdir, basename(arguments$bam))
}
if(isTRUE(arguments$unmated)) {
    count_lines_outfile("u", arguments$outdir, basename(arguments$bam))
}
if(isTRUE(arguments$ambiguous)) {
    count_lines_outfile("a", arguments$outdir, basename(arguments$bam))
}
if(isTRUE(arguments$duplicated)) {
    count_lines_outfile("d", arguments$outdir, basename(arguments$bam))
}
if(isTRUE(arguments$singleton)) {
    count_lines_outfile("s", arguments$outdir, basename(arguments$bam))
}
if(isTRUE(arguments$trans)) {
    count_lines_outfile("t", arguments$outdir, basename(arguments$bam))
}

time_end <- Sys.time()
cat("\n")
cat(paste0(script, " completed: ", time_end, "\n"))
print(convert_time_HMS(time_start, time_end))
cat("\n\n")
rm(time_start, time_end)
