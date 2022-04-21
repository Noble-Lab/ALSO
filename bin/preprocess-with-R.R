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


evaluateMateStatus <- function(x, y, z) {
    # Evaluate mate status for reads in a bam file, reporting if reads with
    # status "mated", "unmated", or "ambiguous" are present; optionally, write
    # out txt.gz tables comprised of 'qname', 'groupid', and 'mate_status' for
    # a given status of "mated", "unmated", or "ambiguous"
    # 
    # :param x: A 3 x 2 tibble or dataframe: the first column is $mate_status
    #           (factor) with levels "mated", "unmated", and "ambiguous"; the
    #           second column is $n with counts (int) for each level
    # :param y: "--mated", "--unmated", or "--ambiguous" argument (logical)
    # :param z: status to test: "mated", "unmated", or "ambiguous" (chr)
    #TODO Check for class(x)...
    
    if (class(y) != "logical") {
        stop("Exiting: Argument 'y' should be class 'logical'.")
    }
    
    if (class(z) != "character") {
        stop("Exiting: Argument 'z' should be class 'character'.")
    }
    
    if(isTRUE(y)) {
        if(x[x$mate_status == z, 2] == 0) {
            print(paste0("No ", z, " reads found."))
        } else {
            print(paste0(
                stringr::str_to_title(z),
                " reads found: ",
                scales::comma(unlist(x[x$mate_status == z, 2])),
                "; writing out table."
            ))
            #  Write tibble to gzipped txt file
            readr::write_tsv(
                pertinent[pertinent$mate_status == z, ],
                paste0(
                    arguments$outdir, "/",
                    gsub(
                        ".bam",
                        paste0(".", z, "-full.txt.gz"),
                        basename(arguments$bam)
                    )
                )
            )
        }
    } else if(isFALSE(y)) {
        if(x[x$mate_status == z, 2] == 0) {
            print(paste0("No ", z, " reads found."))
        } else {
            print(paste0(
                stringr::str_to_title(z),
                " reads found: ",
                scales::comma(unlist(x[x$mate_status == z, 2])),
                "; not writing out table."
            ))
        }
    } else {
        stop(paste0(
            "Stopping: --", z, " is ", y,
            " but must be either TRUE or FALSE."
        ))
    }
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
libraries <- c("argparser", "pryr", "Rsamtools", "scales", "tidyverse")
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
        view -hN'; the user has the options to output txt.gz tables for 'mated'
        reads and, if present in the bam file, 'unmated' and 'ambiguous' reads;
        outfiles are derived from the name of the bam infile.
        
        #TODO Handle duplicate QNAME issue in which a very small number of
        QNAME entries are >2; filter those out in this script or prior to
        reading in the bam file? #ANSWER Prior to reading in the bam file.
        
        #TODO Leverage Rsamtools::BamFile yieldSize option to iterate through
        large bam files; see, for example, page 14 of
        bioconductor.org/packages/devel/bioc/manuals/Rsamtools/man/Rsamtools.pdf
        
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
    short = "-a",
    arg = "--ambiguous",
    type = "logical",
    default = FALSE,
    help = "
        save ambiguous read 'qname', 'groupid', and 'mate_status' table in a
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
    bam <- "Disteche_sample_6.dedup.CAST.sort-c.bam"
    # bam <- "Disteche_sample_7.CAST.processed.chr1.bam"
    # bam <- "Disteche_sample_7.mm10.processed.chr1.bam"
    bai <- "Disteche_sample_6.dedup.CAST.sort-c.bam.bai"
    # bai <- "Disteche_sample_7.CAST.processed.chr1.bam.bai"
    # bai <- "Disteche_sample_7.mm10.processed.chr1.bam.bai"
    mated <- FALSE
    unmated <- FALSE
    ambiguous <- FALSE
    cl <- c(
        #  Arguments for analysis
        "--bam", paste0(dir_in_out, "/", bam),
        "--bai", paste0(dir_in_out, "/", bai),
        "--mated", mated,
        "--unmated", unmated,
        "--ambiguous", ambiguous,
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
    rm(ap)
} else {
    stop(paste0(
        "Stopping: Variable 'test_in_RStudio' is not properly set. ",
        "'test_in_RStudio' should be hardcoded as either TRUE or FALSE."
    ))
}


#  Check that files exist -----------------------------------------------------
stopifnot(file.exists(arguments$bam))
stopifnot(file.exists(arguments$bai))
stopifnot(arguments$mated == TRUE | arguments$mated == FALSE)
stopifnot(arguments$unmated == TRUE | arguments$unmated == FALSE)
stopifnot(arguments$ambiguous == TRUE | arguments$ambiguous == FALSE)

#  If outfile directory  does not exist, then create it
dir.create(file.path(arguments$outdir), showWarnings = FALSE)
#TODO Print message if new directory is created


#  Load in .bam information, including mate information -----------------------
print(paste0(
    "Started: Using Rsamtools to load in '", basename(arguments$bam),
    "' and '", basename(arguments$bai),
    "', and reading 'qname', 'groupid', and 'mate_status' into memory."
))

start_time <- Sys.time()

bam <- Rsamtools::BamFile(arguments$bam, index = arguments$bai, asMates = TRUE)
pertinent <- Rsamtools::scanBam(
    bam, param = ScanBamParam(what = c("qname", "groupid", "mate_status"))
)
# pryr::object_size(pertinent)
#  Disteche_sample_7.CAST.processed.chr1.bam; qname, groupid, mate_status: 992.3 MB
#  Disteche_sample_7.mm10.processed.chr1.bam; qname, groupid, mate_status: 937.8 MB
pertinent <- Map(as.data.frame, pertinent) %>%
    as.data.frame() %>%
    tibble::as_tibble(column_name = c("qname", "groupid", "mate_status"))
# pryr::object_size(pertinent)
#  Disteche_sample_7.CAST.processed.chr1.bam; qname, groupid, mate_status: 992.3 MB
#  Disteche_sample_7.CAST.processed.chr1.bam; qname, groupid, mate_status: 937.8 MB

end_time <- Sys.time()

print(paste0(
    "Completed: Using Rsamtools to load in '", basename(arguments$bam),
    "' and '", basename(arguments$bai),
    "', and reading 'qname', 'groupid', and 'mate_status' into memory. ",
    "Memory in use: ",
    scales::comma(as.integer(pryr::object_size(pertinent))), " bytes. ",
    "Run time: ",
    round(as.numeric(unlist(stringr::str_split((end_time - start_time), " "))), 3),
    " seconds."
))  #TODO Report on time spent for the process too...
cat("\n")
rm(start_time, end_time)


#  Determine and evaluate mate_status levels present in the data --------------
#+ ...then create objects for them
print(paste0(
    "Started: Evaluating and reporting mate_status levels ('mated', ",
    "'unmated', 'ambiguous') present in '", basename(arguments$bam), "'."
))

mate_status <- pertinent$mate_status %>%
    table() %>%
    dplyr::as_tibble() %>%
    dplyr::rename("mate_status" = ".")
# pryr::object_size(mate_status)
#  Disteche_sample_7.CAST.processed.chr1.bam; mate_status: 1,216 B
#  Disteche_sample_7.mm10.processed.chr1.bam; mate_status: 1,216 B

evaluateMateStatus(mate_status, arguments$mated, "mated")
evaluateMateStatus(mate_status, arguments$unmated, "unmated")
evaluateMateStatus(mate_status, arguments$unmated, "ambiguous")

print(paste0(
    "Completed: Evaluating and reporting mate_status levels ('mated', ",
    "'unmated', 'ambiguous') present in '", basename(arguments$bam), "'."
))
cat("\n")


#  Separate "unmated" and "ambiguous" reads from "mated" reads ----------------
print(paste0(
    "Started: Writing output txt file comprised of 'qname' entries for reads ",
    "to retain in '", basename(arguments$bam), "'."
))

#  Remove reads with $mate_status of "unmated"
if(mate_status[mate_status$mate_status == "unmated", 2] != 0) {
    pertinent <- pertinent[pertinent$mate_status != "unmated", ]
}
# pryr::object_size(pertinent)
#  Disteche_sample_7.CAST.processed.chr1.bam; qname, groupid, mate_status: 937.8 MB

#  Remove reads with $mate_status of "ambiguous"
if(mate_status[mate_status$mate_status == "ambiguous", 2] != 0) {
    pertinent <- pertinent[pertinent$mate_status != "ambiguous", ]
}

#  Remove all fields except QNAME
pertinent <- dplyr::select(pertinent, -c(groupid, mate_status))
# pryr::object_size(pertinent)
#  Disteche_sample_7.CAST.processed.chr1.bam; qname only: 882.1 MB
#  Disteche_sample_7.mm10.processed.chr1.bam; qname only: 833.6 MB

#  Remove duplicate QNAME entries, then write tibble to a gzipped txt file
pertinent <- pertinent[seq_len(nrow(pertinent)) %% 2 == 0, ]
# pryr::object_size(pertinent)
#  Disteche_sample_7.CAST.processed.chr1.bam; only deduplicated qname: 826.9 MB
#  Disteche_sample_7.mm10.processed.chr1.bam; only deduplicated qname: 781.5 MB

#  Write out final txt file containing qnames to retain in the preprocessed bam
#+ file
outfile <- gsub(".bam", ".mated-qname-dedup.txt", basename(arguments$bam))
readr::write_tsv(
    pertinent,
    paste0(arguments$outdir, "/", outfile)
)

print(paste0(
    "Completed: Writing output txt file comprised of 'qname' entries for ",
    "reads to retain in '", basename(arguments$bam),
    "'. Number of 'qname' entires: ", scales::comma(nrow(pertinent)), "."
))
cat("\n")


#  In shell... ----------------------------------------------------------------
# bam="Disteche_sample_7.CAST.processed.chr1.bam"
# keep="Disteche_sample_7.CAST.processed.chr1.mated-qname-dedup.txt"
# samtools view -@ 2 -hN "${keep}" "${bam}" | samtools view -@ 2 -b - > "${bam/.bam/.mated-qname-dedup.bam}"
# samtools flagstat -@ 2 "${bam/.bam/.mated-qname-dedup.bam}" > "${bam/.bam/.mated-qname-dedup.flagstat.txt}"

# bam="Disteche_sample_7.mm10.processed.chr1.bam"
# keep="Disteche_sample_7.mm10.processed.chr1.mated-qname-dedup.txt"
# samtools view -@ 2 -hN "${keep}" "${bam}" | samtools view -@ 2 -b - > "${bam/.bam/.mated-qname-dedup.bam}"
# samtools flagstat -@ 2 "${bam/.bam/.mated-qname-dedup.bam}" > "${bam/.bam/.mated-qname-dedup.flagstat.txt}"
