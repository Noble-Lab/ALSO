#!/usr/bin/env Rscript

#  get-AS-per-qname.R
#  KA


script <- "get-AS-per-qname.R"
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


collapse_mates_into_one_row <- function(x, y) {
    # Input a dataframe or tibble in which mate pairs ('mate 1' and 'mate 2')
    # are in immediately subsequent rows ('row n' and 'row n + 1'); output a
    # tibble in which the mate pairs ('mate n' and 'mate n + 1') are in a
    # single row
    # 
    # :param x: tibble/dataframe input (as class character) <chr>
    # :param y: column in input used to sort tibble output <chr>
    # :return z: tibble output in which mate pairs comprise one row each
    if (class(x) != "character") {
        stop("Exiting: Argument 'x' is not class 'character'.")
    } else {
        z <- eval(parse(text = x))
    }
    
    if (class(y) != "character") {
        stop("Exiting: Argument 'y' should be class 'character'.")
    }

    odd.seq <- seq(1, nrow(z), 2)
    even.seq <- (seq(1, nrow(z), 2)) + 1
    
    odd <- z[odd.seq, ]
    even <- z[even.seq, ]
    
    odd <- odd %>% dplyr::mutate(groupid = row_number())
    even <- even %>% dplyr::mutate(groupid = row_number())
    
    z <- dplyr::full_join(odd, even, by = "groupid") %>%
        dplyr::arrange(y) %>%
        dplyr::rename(groupid.x = groupid)
    z$groupid.y <- z$groupid.x
    
    colnames(z) <- stringr::str_replace_all(colnames(z), "\\.x", "\\.odd")
    colnames(z) <- stringr::str_replace_all(colnames(z), "\\.y", "\\.even")
    
    return(z)
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
    # :param x: object of class BamFile
    # :return y: tibble made up of qname, rname, mrnm, groupid, and
    #            mate_status fields
    #TODO Check that BamFile asMates is TRUE
    #TODO Check that BamFile isOpen is TRUE
    vars <- c(
        "qname", "groupid", "mate_status"
    )
    y <- Rsamtools::scanBam(x, param = ScanBamParam(what = vars, tag = "AS"))
    y <- Map(as.data.frame, y) %>%
        as.data.frame() %>%
        tibble::as_tibble(column_name = vars)
    return(y)
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
        Output a sorted, tab-separated, gzipped table of qnames and AS's for a
        given bam infile.
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
    help = "directory for saving rds outfile, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-s",
    arg = "--strain",
    type = "character",
    default = NULL,
    help = "strain name to be appended to rds outfile columns <chr>"
)
ap <- add_argument(
    ap,
    short = "-c",
    arg = "--chunk",
    type = "integer",
    default = 100000,
    help = "number of records to read into memory at a single time <even int>"
)


#  Parse the arguments --------------------------------------------------------
test_in_RStudio <- FALSE  # Hardcode TRUE for interactive testing; FALSE for CLI
if(isTRUE(test_in_RStudio)) {
    #  RStudio-interactive work
    dir_base <- "."
    dir_data <- "data"
    dir_in <- paste0(dir_base, "/", dir_data, "/", "files_bam")
    dir_out <- paste0(dir_base, "/", dir_data, "/", "files_bam")
    # bam <- paste0(dir_in, "/", "Disteche_sample_1.dedup.CAST.corrected.bam")
    # bai <- paste0(dir_in, "/", "Disteche_sample_1.dedup.CAST.corrected.bam.bai")
    # strain <- "CAST"
    bam <- paste0(dir_in, "/", "Disteche_sample_1.dedup.CAST.corrected.bam")
    bai <- paste0(dir_in, "/", "Disteche_sample_1.dedup.CAST.corrected.bam.bai")
    strain <- "CAST"
    chunk <- 100000
    cl <- c(
        "--bam", bam,
        "--bai", bai,
        "--strain", strain,
        "--outdir", dir_out,
        "--chunk", chunk
    )
    arguments <- parse_args(ap, cl)
    rm(ap, bam, bai, cl, chunk, dir_base, dir_data, dir_in, dir_out, strain)
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
if(arguments$strain == "") stop("--strain is an empty string.")
if(is.na(arguments$strain)) stop("--strain is NA.")
stopifnot(arguments$chunk != 0)
stopifnot(arguments$chunk %% 1 == 0)
stopifnot(arguments$chunk %% 2 == 0)

#  If it does not exist, then create outfile directory 
dir.create(file.path(arguments$outdir), showWarnings = FALSE)
#TODO Print message if new directory is created


#  Set up variables, environment prior to loading in .bam information... ------
#+ ...including mate information

#  Add strain name to the file if it's not already present; save it to
#+ "outname" and "outname_tmp"
if(grepl(arguments$strain, basename(arguments$bam), fixed = TRUE)) {
    outstring <- paste0(stringr::str_remove(basename(arguments$bam), ".bam"))
    
    outname <- paste0(outstring, ".AS.txt.gz")
    outname_tmp <- paste0(outstring, ".AS.tmp.txt.gz")
    
    rm(outstring)
} else {
    outstring <- paste0(".", arguments$strain, ".AS.txt.gz")
    outstring_tmp <- paste0(".", arguments$strain, ".AS.tmp.txt.gz")
    
    outname <- gsub(".bam", outstring, basename(arguments$bam))
    outname_tmp <- gsub(".bam", outstring_tmp, basename(arguments$bam))
    
    rm(outstring, outstring_tmp)
}

cat(paste0(
    "Using Rsamtools to load in '", basename(arguments$bam),
    "' and reading fields such as 'qname' into memory (in chunks of ",
    scales::comma(arguments$chunk), " records per iteration).\n"
))
cat("\n")

#  Tally and record the total number of records in the bam file
cat(paste0(
    "Counting the number of records in ", basename(arguments$bam), "...\n"
))
rec_n <- as.integer(arguments$chunk)
rec_total <- count_records(arguments$bam)
cat(paste0("Number of records: ", scales::comma(rec_total), "\n"))
cat("\n")

#  Iterate through bam file in chunks
bam <- Rsamtools::BamFile(arguments$bam, index = arguments$bai, asMates = TRUE)
Rsamtools::yieldSize(bam) <- arguments$chunk
open(bam)

cat(paste0("Started: Processing ", basename(arguments$bam), "\n"))
n <- ceiling((rec_total / 2) / rec_n) %>% as.integer()
bar <- utils::txtProgressBar(min = 0, max = n, initial = 0, style = 3)

for(i in 1:n) {
    #  Using Rsamtools, load in qname, mate_status, and AS fields
    pertinent <- load_fields(bam)
    
    #  Determine which mate_status levels are present in the data, then create
    #+ objects for them
    mate_status <- pertinent$mate_status %>%
        table() %>%
        dplyr::as_tibble() %>%
        dplyr::rename("mate_status" = ".")
    
    #  Throw warning if "unmated" reads are found
    if(mate_status[mate_status$mate_status == "unmated", 2] != 0) {
        cat(paste0(
            "WARNING: Iteration ", i, "/", n, ": ",
            mate_status[mate_status$mate_status == "unmated", 2],
            " 'unmated' reads detected.\n"
        ))
    }
    
    #  Throw warning if "ambiguous" reads are found
    if(mate_status[mate_status$mate_status == "ambiguous", 2] != 0) {
        cat(paste0(
            "WARNING: Iteration ", i, "/", n, ": ",
            mate_status[mate_status$mate_status == "ambiguous", 2],
            " 'ambiguous' reads detected.\n"
        ))
    }
    
    #  Throw warning if no "mated" reads are found; otherwise, collect mated
    #+ read fields into a tibble
    if(mate_status[mate_status$mate_status == "mated", 2] == 0) {
        cat(paste0(
            "WARNING: Iteration ", i, "/", n, ": ",
            "No 'mated' reads detected.\n"
        ))
    } else if(mate_status[mate_status$mate_status == "mated", 2] > 0) {
        pertinent <- pertinent[pertinent$mate_status == "mated", ]
        
        #  Clean up the tibble
        pertinent <- dplyr::select(pertinent, -mate_status)
        colnames(pertinent)[-1] <- paste0(
            colnames(pertinent)[-1], ".", arguments$strain
        )
        
        #  Organize mated reads into one row per mate pair
        pertinent <- collapse_mates_into_one_row("pertinent", "pos") %>%
            dplyr::select(-c(qname.even, groupid.odd, groupid.even)) %>%
            dplyr::rename(qname = qname.odd) %>%
            dplyr::select(-c(colnames(.)[3], colnames(.)[5]))
            
        #  Take the AS min for each mate pair: Taking the AS min via pmin()
        #+ gives us the worst AS score per read pair; our logic is that the
        #+ read pair should be represented by the worst of its two alignment
        #+ scores
        pertinent$pmin <- pmin(pertinent[, 2][[1]], pertinent[, 3][[1]])
        
        #  Rename, reorder, remove columns
        string <- colnames(pertinent)[[2]] %>%
            stringr::str_split("\\.") %>%
            unlist()
        colnames(pertinent)[4] <- paste0(string[1], ".", string[2], ".pmin")
        pertinent <- pertinent %>% dplyr::select(-c(
            colnames(pertinent)[2], colnames(pertinent)[3]
        ))

        #  Write out QNAMEs and AS's to txt.gz file
        readr::write_tsv(
            pertinent,
            paste0(arguments$outdir, "/", outname),
            append = TRUE
        )
    }
    utils::setTxtProgressBar(bar, i)
}


cat(paste0("Completed: Processing ", arguments$bam, "\n"))
cat(paste0("Have written out ", arguments$outdir, "/", outname, "\n"))


#  End the script -------------------------------------------------------------
time_end <- Sys.time()
cat("\n")
cat(paste0(script, " completed: ", time_end, "\n"))
print(convert_time_HMS(time_start, time_end))
cat("\n\n")
rm(time_start, time_end)
