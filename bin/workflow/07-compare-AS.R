#!/usr/bin/env Rscript

#  compare-AS.R
#  KA


#  Functions ------------------------------------------------------------------
evaluateOperation <- function(x) {
    # Evaluate a string as an R operation
    # 
    # :param x: operation string to be evaluated <chr>
    # :retrun x: evaluated operation string
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
    # Make a string that will be evaluated as an R operation
    # 
    # :param x: for operation string, variable to be assigned to <chr>
    # :param y: for operation string, assignment <chr>
    # :return operation: an R operation
    operation <- paste0(x, " <- ", y)
    return(operation)
}


mungeMatesIntoTwoRows <- function(x, y) {
    # Input a tibble in which mate pairs are in a single row; output a tibble
    # in which in which the mate pairs are in immediately subsequent rows
    # ('row n' and 'row n + 1')
    # 
    # :param x: tibble/dataframe input (as class character) <chr>
    # :param y: column in input used to sort tibble output <chr>
    # :return z: tibble output in which mate pairs comprise immediately
    #            subsequent rows
    
    #  All function arguments should be of class "character"
    if (class(x) != "character") {
        stop("Exiting: Argument 'x' is not class 'character'.")
    } else {
        z <- eval(parse(text = x))
    }

    if (class(y) != "character") {
        stop("Exiting: Argument 'y' should be class 'character'.")
    }

    #  Sort the "uniline tibbles" by...
    z <- z %>% dplyr::arrange(y)

    #  Split the "uniline tibbles" by odd or even status
    odd <- z[stringr::str_subset(colnames(z), "\\.odd")]
    even <- z[stringr::str_subset(colnames(z), "\\.even")]

    #  Strip suffixes from column names
    colnames(odd) <- stringr::str_replace_all(colnames(odd), "\\.odd", "")
    colnames(even) <- stringr::str_replace_all(colnames(even), "\\.even", "")

    #  Set odd tibble to tibble 1, even tibble to tibble 2
    odd$tibble <- "1"
    even$tibble <- "2"

    #  Interleave the odd/even rows, then arrange them by group ID and
    #+ tibble number; return the interleaved "biline tibble"
    z <- odd %>%
        dplyr::mutate(groupid = row_number()) %>%
        dplyr::bind_rows(even %>% mutate(groupid = row_number())) %>%
        dplyr::arrange(groupid, tibble) %>%
        dplyr::select(-tibble)

    return(z)
}


#  Source libraries, adjust settings ------------------------------------------
importLibrary(c("argparser", "Rsamtools", "tidyverse"))

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)


#  Parse arguments ------------------------------------------------------------
#  Create a parser
ap <- arg_parser(
    name = "compare-AS.R",
    description = "",
    hide.opts = TRUE
)

#  Add command line arguments
ap <- add_argument(
    ap,
    short = "-1",
    arg = "--rds_1",
    type = "character",
    default = NULL,
    help = "first of two rds infiles to compare, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-2",
    arg = "--rds_2",
    type = "character",
    default = NULL,
    help = "second of two rds infiles to compare, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-t",
    arg = "--threshold",
    type = "character",
    default = 0,
    help = "positive integer for thresholding when assigning categories <int>"
)
ap <- add_argument(
    ap,
    short = "-o",
    arg = "--outfile",
    type = "character",
    default = NULL,
    help = "name of txt tensor outfile, including path <chr>"
)

#  Parse the command line arguments
dir_base <- "."
dir_data <- "data"
dir_in <- paste0(
    dir_base, "/", dir_data, "/", "2022-0404_prepare_test-datasets_mm10-CAST"
)
dir_out <- paste0(
    dir_base, "/", dir_data, "/", "2022-0404_prepare_test-datasets_mm10-CAST"
)
rds_1 <- paste0(dir_in, "/", "Disteche_sample_6.dedup.mm10.prepro.chr19.repair.mated.rds")
rds_2 <- paste0(dir_in, "/", "Disteche_sample_6.dedup.CAST.prepro.chr19.repair.mated.rds")
threshold <- 0
outfile <- paste0(dir_out, "/", "Disteche_sample_6.AS.mm10-vs-CAST.txt")
cl <- c(
    "--rds_1", rds_1,
    "--rds_2", rds_2,
    "--threshold", threshold,
    "--outfile", outfile
)
arguments <- parse_args(ap, cl)  # RStudio-interactive work
# arguments <- parse_args(ap)  # Command-line calls

rm(
    ap, cl, dir_base, dir_data, dir_in, dir_out,
    outfile, rds_1, rds_2, threshold
)


#  Check that files exist -----------------------------------------------------
stopifnot(file.exists(arguments$rds_1))
stopifnot(file.exists(arguments$rds_2))
stopifnot(arguments$threshold >= 0)
# stopifnot(any(arguments$strain %in% c("CAST", "mm10", "SPRET", "CAROLI")))

#  If it does not exist, then create outfile directory 
dir.create(file.path(dirname(arguments$outfile)), showWarnings = FALSE)
#TODO Print message if new directory is created


#  Load in rds files ----------------------------------------------------------
rds_1 <- readRDS(arguments$rds_1)
rds_2 <- readRDS(arguments$rds_2)


#  Take the AS min for each mate pair -----------------------------------------
#  Taking the AS min via pmin() gives us the worst AS score per read pair; our
#+ logic is that the read pair should be represented by the worst of its two
#+ AS's

#  rds_1
AS_1 <- rds_1[, grep("AS\\.", colnames(rds_1))]
AS_1$AS_pmin <- pmin(AS_1[, 1][[1]], AS_1[, 2][[1]])

string <- colnames(AS_1)[[1]] %>% stringr::str_split("\\.") %>% unlist()
colnames(AS_1)[3] <- paste0(colnames(AS_1)[3], ".", string[2])

AS_1$qname <- rds_1$qname.odd

rds_1[, length(colnames(rds_1)) + 1] <- AS_1[, 3]

#  rds_2
AS_2 <- rds_2[, grep("AS\\.", colnames(rds_2))]
AS_2$AS_pmin <- pmin(AS_2[, 1][[1]], AS_2[, 2][[1]])

string <- colnames(AS_2)[[1]] %>% stringr::str_split("\\.") %>% unlist()
colnames(AS_2)[3] <- paste0(colnames(AS_2)[3], ".", string[2])

AS_2$qname <- rds_2$qname.odd

rds_2[, length(colnames(rds_2)) + 1] <- AS_2[, 3]


#  Compare AS min for each rds, assigning strain based on superior score ------
AS_compare <- dplyr::full_join(AS_1[, c(3, 4)], AS_2[, c(3, 4)]) %>%
    dplyr::select("qname", dplyr::everything())

#  Initially assign categories to the difference based on values with regards
#+ to the variable int (see below)
int <- arguments$threshold  #TODO Make into an argument
AS_compare <- AS_compare %>% 
    dplyr::mutate(difference = .[[2]] - .[[3]]) %>% 
    dplyr::mutate(
        assignment.initial = case_when(
            difference >= (-1 * int) & difference <= int ~ "Ambiguous",
            difference > int ~ gsub("AS_pmin\\.", "", colnames(AS_compare[2])),
            difference < (-1 * int) ~ gsub("AS_pmin\\.", "", colnames(AS_compare[3]))
        )
    )
rm(int)

AS_compare$assignment <- ifelse(
    is.na(AS_compare$assignment.initial),
    ifelse(
        !is.na(
            eval(parse(text = paste0("AS_compare$", colnames(AS_compare[, 2]))))
        ),
        gsub("AS_pmin\\.", "", colnames(AS_compare[2])),
        ifelse(
            !is.na(
                eval(parse(text = paste0("AS_compare$", colnames(AS_compare[, 3]))))
            ),
            gsub("AS_pmin\\.", "", colnames(AS_compare[3])),
            AS_compare$assignment.initial
        )
    ),
    AS_compare$assignment.initial
) %>%
    forcats::as_factor()

AS_decision <- AS_compare[, c(1, 6)]
table(AS_decision$assignment, useNA = "always")


#  Write out text file for AS decisions
readr::write_tsv(
    AS_decision,
    arguments$outfile,
    col_names = TRUE
)

print("Script has completed.")
