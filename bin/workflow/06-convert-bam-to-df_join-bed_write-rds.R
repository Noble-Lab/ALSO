#!/usr/bin/env Rscript

#  convert-bam-to-df_join-bed_write-rds.R
#  KA


#  Functions ------------------------------------------------------------------
assignFileToVariables <- function(x) {
    # Assign a file to a variable that derives its name from the file
    # 
    # :param x: file <chr>
    # :return name: #TODO
    split <- x %>% strsplit(., "/")
    name <- split %>% lapply(., `[[`, length(split[[1]])) %>% unlist()
    return(name)
}


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


mungeMatesIntoOneRow <- function(x, y) {
    # Input a dataframe/tibble in which mate pairs ('mate 1' and 'mate 2')
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
    
    colnames(z) <- str_replace_all(colnames(z), "\\.x", "\\.odd")
    colnames(z) <- str_replace_all(colnames(z), "\\.y", "\\.even")
    
    return(z)
}


#TODO Function is not used here; move function to script where it is used
# mungeMatesIntoTwoRows <- function(x, y) {
#     #  All function arguments should be of class "character"
#     if (class(x) != "character") {
#         stop("Exiting: Argument 'x' is not class 'character'.")
#     } else {
#         z <- eval(parse(text = x))
#     }
#    
#     if (class(y) != "character") {
#         stop("Exiting: Argument 'y' should be class 'character'.")
#     }
#    
#     #  Sort the "uniline tibbles" by...
#     z <- z %>% dplyr::arrange(y)
#
#     #  Split the "uniline tibbles" by odd or even status
#     odd <- z[stringr::str_subset(colnames(z), "\\.odd")]
#     even <- z[stringr::str_subset(colnames(z), "\\.even")]
#
#     #  Strip suffixes from column names
#     colnames(odd) <- stringr::str_replace_all(colnames(odd), "\\.odd", "")
#     colnames(even) <- stringr::str_replace_all(colnames(even), "\\.even", "")
#
#     #  Set odd tibble to tibble 1, even tibble to tibble 2
#     odd$tibble <- "1"
#     even$tibble <- "2"
#
#     #  Interleave the odd/even rows, then arrange them by group ID and
#     #+ tibble number; return the interleaved "biline tibble"
#     z <- odd %>%
#         dplyr::mutate(groupid = row_number()) %>%
#         dplyr::bind_rows(even %>% mutate(groupid = row_number())) %>%
#         dplyr::arrange(groupid, tibble) %>%
#         dplyr::select(-tibble)
#
#     return(z)
# }


#  Source libraries, adjust settings ------------------------------------------
importLibrary(c("argparser", "Rsamtools", "tidyverse"))

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)


# #  Set up working directory ---------------------------------------------------
# project <- "2021_kga0_4dn-mouse-cross"
# current <- stringr::str_split(getwd(), "/")[[1]][
#     length(stringr::str_split(getwd(), "/")[[1]])
# ]
# if(current != project) {
#     if(dir.exists(project)) {
#         setwd(project)
#     } else {
#         setwd(
#             readline(
#                 prompt = paste0(
#                     "Enter path to and including the project directory, ",
#                     project, ":"
#                 )
#             )  #TODO Check again: current == project?
#         )
#     }
# }
# rm(project, current)


#  Parse arguments ------------------------------------------------------------
#  Create a parser
ap <- arg_parser(
    name = "convert-bam-to-df_join-bed_write-rds.R",
    description = "",
    hide.opts = TRUE
)

#  Add command line arguments
ap <- add_argument(
    ap,
    short = "-i",
    arg = "--bam",
    type = "character",
    default = NULL,
    help = "bam infile, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-b",
    arg = "--bai",
    type = "character",
    default = NULL,
    help = "bam index, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-p",
    arg = "--pos",
    type = "character",
    default = NULL,
    help = "liftOver bed file for \"POS\", including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-m",
    arg = "--mpos",
    type = "character",
    default = NULL,
    help = "liftOver bed file for \"MPOS\", including path <chr>"
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
    short = "-r",
    arg = "--rds",
    type = "character",
    default = NULL,
    help = "name of rds outfile to be saved in outdir <chr>"
)

# #  Parse the command line arguments
# dir_base <- "."
# dir_data <- "data"
# dir_in <- paste0(
#     dir_base, "/", dir_data, "/", "2022-0320_test_04-05_all"
# ) 
# dir_out <- paste0(
#     dir_base, "/", dir_data, "/", "2022-0323_test_06_chr19"
# )
# bam <- paste0(dir_in, "/", "test.300000.chr19.bam")
# bai <- paste0(dir_in, "/", "test.300000.chr19.bam.bai")
# pos <- paste0(dir_in, "/", "test.300000.chr19.pos.liftOver.CAST-EiJ.bed")
# mpos <- paste0(dir_in, "/", "test.300000.chr19.mpos.liftOver.CAST-EiJ.bed")
# rds <- "test.300000.chr19.rds"
# cl <- c(
#     "--bam", bam,
#     "--bai", bai,
#     "--pos", pos,
#     "--mpos", mpos,
#     "--outdir", dir_out,
#     "--rds", rds
# )
# arguments <- parse_args(ap, cl)  # RStudio-interactive work
arguments <- parse_args(ap)  # Command-line calls

rm(dir_base, dir_data, dir_in, dir_out, bam, bai, pos, mpos, rds, ap)


#  Check that files exist -----------------------------------------------------
stopifnot(file.exists(arguments$bam))
stopifnot(file.exists(arguments$bai))
stopifnot(file.exists(arguments$pos))
stopifnot(file.exists(arguments$mpos))

#  If it does not exist, then create outfile directory 
dir.create(file.path(arguments$outdir), showWarnings = FALSE)


#  Load in bam information, including mate information ------------------------
bam <- assignFileToVariables(arguments$bam)
mapply(
    assign, bam, arguments$bam, MoreArgs = list(envir = parent.frame())
)

index <- assignFileToVariables(arguments$bai)
mapply(
    assign, index, arguments$bai, MoreArgs = list(envir = parent.frame())
)

#  Using Rsamtools, load in standard bam fields
command <- paste0(
    "Rsamtools::BamFile(",
        bam, ", index = ", index, ", asMates = TRUE",
    ")", " %>% ",
    "Rsamtools::scanBam()", " %>% ",
    "as.data.frame()", " %>% ",
    "tibble::as_tibble()"
)
operation <- makeOperation("full", command)
evaluateOperation(operation)

#  Using Rsamtools, load in bam AS and MD tags
map_params <- Rsamtools::ScanBamParam(tag = c("AS", "MD"))
command <- paste0(
    bam, " %>% ",
        "Rsamtools::BamFile(",
            bam, ", index = ", index, ", asMates = TRUE",
        ")", " %>% ",
        "Rsamtools::scanBam(param = map_params)", " %>% ",
        "as.data.frame()", " %>% ",
        "tibble::as_tibble()"
)
operation <- makeOperation("AS_MD", command)
evaluateOperation(operation)

rm(map_params)

#  Column-bind the AS and MD tags to the bam-file tibble
command <- paste0(
    "dplyr::bind_cols(", "full, ", "AS_MD", ")"
)
operation <- makeOperation(bam, command)
evaluateOperation(operation)

#  Munge the bam-file tibble
command <- paste0(
    "dplyr::select(", bam, ", ", "-c(groupid, mate_status)", ")"
)
operation <- makeOperation(bam, command)
evaluateOperation(operation)

#  Remove unneeded variables
operation <- paste0("rm(", "full, ", "AS_MD", ")")
evaluateOperation(operation)

rm(index)

#  Assign chromosome variable needed for subsequent steps
command <- paste0(bam, "$rname[1] %>% as.character()")
operation <- makeOperation("chromosome", command)
evaluateOperation(operation)


#  Join the POS and MPOS bed files to the bam-file tibble ---------------------
pos <- readr::read_tsv(
    arguments$pos,
    col_names = c(
        "lO_rname", "lO_pos", "lO_pos_end", "qname", "lO_pos_status"
    )
)

mpos <- readr::read_tsv(
    arguments$mpos,
    col_names = c(
        "lO_mrnm", "lO_mpos", "lO_mpos_end", "qname", "lO_mpos_status"
    )
)

pos_mpos <- dplyr::full_join(pos, mpos, by = "qname") %>%
    dplyr::rename("qname.odd" = "qname")

#  Remove unneeded variables
rm(pos, mpos)

#  Get the mates into a tibble with one row per pair
command <- paste0("mungeMatesIntoOneRow(bam, \"pos\")")
operation <- operation <- makeOperation("rows_one", command)
evaluateOperation(operation)

#  Join the liftOver POS/MPOS information to the one-row-per-pair tibble; then,
#+ reorganize the tibble for easier readability 
join <- dplyr::full_join(rows_one, pos_mpos, by = "qname.odd") %>%
    dplyr::select(-c(lO_pos_end, lO_mpos_end)) %>%
    dplyr::mutate(
        lO_rname.even = lO_rname,
        lO_pos.even = lO_pos,
        lO_pos_status.even = lO_pos_status,
        lO_mrnm.even = lO_mrnm,
        lO_mpos.even = lO_mpos,
        lO_mpos_status.even = lO_mpos_status
    ) %>%
    dplyr::rename(
        lO_rname.odd = lO_rname,
        lO_pos.odd = lO_pos,
        lO_pos_status.odd = lO_pos_status,
        lO_mrnm.odd = lO_mrnm,
        lO_mpos.odd = lO_mpos,
        lO_mpos_status.odd = lO_mpos_status
    ) %>%
    dplyr::relocate(
        c(
            lO_rname.odd, lO_pos.odd, lO_pos_status.odd,
            lO_mrnm.odd, lO_mpos.odd, lO_mpos_status.odd,
            rname.odd, pos.odd, mrnm.odd, mpos.odd, isize.odd
        ),
        .after = strand.odd
    ) %>%
    dplyr::relocate(
        c(
            lO_rname.even, lO_pos.even, lO_pos_status.even,
            lO_mrnm.even, lO_mpos.even, lO_mpos_status.even,
            rname.even, pos.even, mrnm.even, mpos.even, isize.even
        ),
        .after = strand.even
    )

#  Remove unneeded variables
rm(pos_mpos, rows_one)

operation <- paste0("rm(", bam, ")")
evaluateOperation(operation)


#  Retain pairs with successful liftOver status and matching chromosomes ------
retain <- join %>%
    subset(
        lO_pos_status.odd == "liftOver successful" &
        lO_mpos_status.odd == "liftOver successful" &
        lO_pos_status.even == "liftOver successful" &
        lO_mpos_status.even == "liftOver successful"
    ) %>%
    subset(
        lO_rname.odd == chromosome &
        lO_mrnm.odd == chromosome &
        lO_rname.even == chromosome &
        lO_mrnm.even == chromosome
    )

#TODO 1/2 Save this information for, for example, reporting metrics on reads
#TODO 2/2 that were lost due to liftOver?
# lost_1 <- join %>%
#     subset(
#         .,
#         lO_pos_status.odd != "liftOver successful" |
#         lO_mpos_status.odd != "liftOver successful" |
#         lO_pos_status.even != "liftOver successful" |
#         lO_mpos_status.even != "liftOver successful"
#     )
# lost_2 <- join %>%
#     subset(
#         lO_rname.odd != chromosome |
#         lO_mrnm.odd != chromosome |
#         lO_rname.even != chromosome |
#         lO_mrnm.even != chromosome
#     )

#  Remove unneeded variables, etc.
rm(join, bam, chromosome, cl, command, operation, script)


#  Write out rds file for munged pairs information ----------------------------
saveRDS(retain, file = paste0(arguments$outdir, "/", arguments$rds))
