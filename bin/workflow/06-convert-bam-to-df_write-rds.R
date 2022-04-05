#!/usr/bin/env Rscript

#  convert-bam-to-df_write-rds.R
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


#  Source libraries, adjust settings ------------------------------------------
importLibrary(c("argparser", "Rsamtools", "tidyverse"))

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)


#  Parse arguments ------------------------------------------------------------
#  Create a parser
ap <- arg_parser(
    name = "convert-bam-to-df_write-rds.R",
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
# ap <- add_argument(
#     ap,
#     short = "-p",
#     arg = "--pos",
#     type = "character",
#     default = NULL,
#     help = "liftOver bed file for \"POS\", including path <chr>"
# )
# ap <- add_argument(
#     ap,
#     short = "-m",
#     arg = "--mpos",
#     type = "character",
#     default = NULL,
#     help = "liftOver bed file for \"MPOS\", including path <chr>"
# )
ap <- add_argument(
    ap,
    short = "-s",
    arg = "--strain",
    type = "character",
    default = NULL,
    help = "
    strain name to be appended to rds outfile columns; current options are
    'mm10', 'CAST', 'SPRET', and 'CAROLI' <chr>
    "
)
ap <- add_argument(
    ap,
    short = "-o",
    arg = "--outdir",
    type = "character",
    default = NULL,
    help = "directory for saving rds outfile, including path <chr>"
)
# ap <- add_argument(
#     ap,
#     short = "-r",
#     arg = "--rds",
#     type = "character",
#     default = NULL,
#     help = "name of rds outfile to be saved in outdir <chr>"
# )
ap <- add_argument(
    ap,
    short = "-m",
    arg = "--mated",
    type = "logical",
    default = TRUE,
    help = "write rds file for mated reads <logical>"
)
ap <- add_argument(
    ap,
    short = "-u",
    arg = "--unmated",
    type = "logical",
    default = FALSE,
    help = "write rds file for unmated reads <logical>"
)
ap <- add_argument(
    ap,
    short = "-a",
    arg = "--ambiguous",
    type = "logical",
    default = FALSE,
    help = "write rds file for ambiguous reads <logical>"
)

# #  Parse the command line arguments
# dir_base <- "."
# dir_data <- "data"
# dir_in <- paste0(
#     dir_base, "/", dir_data, "/", "2022-0404_prepare_test-datasets_mm10-CAST"
# )
# dir_out <- paste0(
#     dir_base, "/", dir_data, "/", "2022-0404_prepare_test-datasets_mm10-CAST"
# )
# bam <- paste0(dir_in, "/", "Disteche_sample_6.dedup.mm10.prepro.chr19.repair.bam")
# bai <- paste0(dir_in, "/", "Disteche_sample_6.dedup.mm10.prepro.chr19.bam.bai")
# # pos <- paste0(dir_in, "/", "")
# # mpos <- paste0(dir_in, "/", "")
# strain <- "CAST"
# mated <- TRUE
# unmated <- TRUE
# ambiguous <- TRUE
# # rds <- "Disteche_sample_6.dedup.mm10.prepro.chr19.rds"
# cl <- c(
#     "--bam", bam,
#     "--bai", bai,
#     # "--pos", pos,
#     # "--mpos", mpos,
#     "--strain", strain,
#     "--outdir", dir_out,
#     # "--rds", rds
#     "--mated", mated,
#     "--unmated", unmated,
#     "--ambiguous", ambiguous
# )
# arguments <- parse_args(ap, cl)  # RStudio-interactive work
arguments <- parse_args(ap)  # Command-line calls

# # rm(dir_base, dir_data, dir_in, dir_out, bam, bai, rds, ap)
# rm(
#     ambiguous, ap, bam, bai, cl, dir_base, dir_data,
#     dir_in, dir_out, mated, strain, unmated
# )


#  Check that files exist -----------------------------------------------------
stopifnot(file.exists(arguments$bam))
stopifnot(file.exists(arguments$bai))
stopifnot(any(arguments$strain %in% c("CAST", "mm10", "SPRET", "CAROLI")))
stopifnot(arguments$mated == TRUE | arguments$mated == FALSE)
stopifnot(arguments$unmated == TRUE | arguments$unmated == FALSE)
stopifnot(arguments$ambiguous == TRUE | arguments$ambiguous == FALSE)

#  If it does not exist, then create outfile directory 
dir.create(file.path(arguments$outdir), showWarnings = FALSE)
#TODO Print message if new directory is created


#  Load in bam information, including mate information ------------------------
# bam <- assignFileToVariables(arguments$bam)
# mapply(
#     assign, bam, arguments$bam, MoreArgs = list(envir = parent.frame())
# )
# 
# index <- assignFileToVariables(arguments$bai)
# mapply(
#     assign, index, arguments$bai, MoreArgs = list(envir = parent.frame())
# )

file_bam <- arguments$bam
file_index <- arguments$bai

#  Using Rsamtools, load in standard bam fields; include ticks around objects
#+ to handle non-standard characters such as dashes
# command <- paste0(
#     "Rsamtools::BamFile(",
#         bam, ", index = ", index, ", asMates = TRUE",
#     ")", " %>% ",
#     "Rsamtools::scanBam()", " %>% ",
#     "as.data.frame()", " %>% ",
#     "tibble::as_tibble()"
# )
# operation <- makeOperation("full", command)
# evaluateOperation(operation)

#  Using Rsamtools, load in bam AS and MD tags
# map_params <- Rsamtools::ScanBamParam(tag = c("AS", "MD"))
# command <- paste0(
#     bam, " %>% ",
#         "Rsamtools::BamFile(",
#             bam, ", index = ", index, ", asMates = TRUE",
#         ")", " %>% ",
#         "Rsamtools::scanBam(param = map_params)", " %>% ",
#         "as.data.frame()", " %>% ",
#         "tibble::as_tibble()"
# )
# operation <- makeOperation("AS_MD", command)
# suppressWarnings(evaluateOperation(operation))
# #FIXME 1/2 "unnamed arguments to $new() must be objects from a reference class;
# #FIXME 2/2 got an object of class 'chatacter'"

bam <- Rsamtools::BamFile(file_bam, index = file_index, asMates = TRUE) %>% 
    Rsamtools::scanBam() %>%
    as.data.frame() %>%
    tibble::as_tibble()

map_params <- Rsamtools::ScanBamParam(tag = c("AS", "MD"))
AS_MD <- Rsamtools::BamFile(file_bam, index = file_index, asMates = TRUE) %>%
    Rsamtools::scanBam(param = map_params) %>%
    as.data.frame() %>%
    tibble::as_tibble()

#  Column-bind the AS and MD tags to the bam-file tibble
bam <- dplyr::bind_cols(bam, AS_MD)

rm(AS_MD, file_bam, file_index, map_params)

#  Determine which mate_status levels are present in the data, then create
#+ objects for those found
mate_status <- bam$mate_status %>%
    table() %>%
    dplyr::as_tibble() %>%
    dplyr::rename("mate_status" = ".")

if(mate_status[mate_status$mate_status == "mated", 2] == 0) {
    print(paste0"No mated reads found.")
} else if(mate_status[mate_status$mate_status == "mated", 2] > 0) {
    print(paste0("Mated reads found: ", length(mate_status[mate_status$mate_status == "mated", 2])))
    mated <- bam[bam$mate_status == "mated", ]
    
    #NOTE 1/2 Following step assumes user is working with bam file that has been 
    #NOTE 2/2 split by chromosome
    #  Assign chromosome variable needed for subsequent steps
    chromosome <- mated$rname[1] %>% as.character()
    
    #  Clean up the tibble
    mated <- dplyr::select(mated, -c(groupid, mate_status))
    colnames(mated)[-1] <- paste0(colnames(mated)[-1], ".", arguments$strain)
    colnames(mated) <- gsub("tag.", "", colnames(mated))
    
    #  Organize mated reads into one row per mate pair
    mated.rows_one <- mungeMatesIntoOneRow("mated", "pos")
    
    if(arguments$mated == TRUE) {
        saveRDS(
            mated.rows_one,
            file = paste0(
                arguments$outdir, "/",
                gsub(".bam", ".mated.rds", basename(arguments$bam))
            )
        )
        
        print(
            paste0(
                "Outfile '",
                gsub(".bam", ".mated.rds", basename(arguments$bam)),
                "' saved to '",
                arguments$outdir, "'"
            )
        )
    } else {
        print("However, no mated reads written to rds file.")
    }
}

if(mate_status[mate_status$mate_status == "unmated", 2] == 0) {
    print("No unmated reads found.")
} else if(mate_status[mate_status$mate_status == "unmated", 2] > 0) {
    print("Unmated reads found.")
    unmated <- bam[bam$mate_status == "unmated", ]
    
    #  Clean up the tibble
    unmated <- dplyr::select(unmated, -c(groupid, mate_status))
    colnames(unmated)[-1] <- paste0(colnames(unmated)[-1], ".", arguments$strain)
    colnames(unmated) <- gsub("tag", "", colnames(unmated))
    
    if(arguments$unmated == TRUE) {
        saveRDS(
            unmated,
            file = paste0(
                arguments$outdir, "/",
                gsub(".bam", ".unmated.rds", basename(arguments$bam))
            )
        )
        
        print(
            paste0(
                "Outfile '",
                gsub(".bam", ".unmated.rds", basename(arguments$bam)),
                "' saved to '",
                arguments$outdir, "'"
            )
        )
    } else {
        print("However, no unmated reads written to rds file.")
    }
}

if(mate_status[mate_status$mate_status == "ambiguous", 2] == 0) {
    print("No ambiguous reads found.")
} else if(mate_status[mate_status$mate_status == "ambiguous", 2] > 0) {
    print("Ambiguous reads found.")
    ambiguous <- bam[bam$mate_status == "ambiguous", ]
    
    #  Clean up the tibble
    ambiguous <- dplyr::select(ambiguous, -c(groupid, mate_status))
    colnames(ambiguous)[-1] <- paste0(colnames(ambiguous)[-1], ".", arguments$strain)
    colnames(ambiguous) <- gsub("tag", "", colnames(ambiguous))
    
    if(arguments$ambiguous == TRUE) {
        saveRDS(
            ambiguous,
            file = paste0(
                arguments$outdir, "/",
                gsub(".bam", ".ambiguous.rds", basename(arguments$bam))
            )
        )
        
        print(
            paste0(
                "Script completed: outfile '",
                gsub(".bam", ".ambiguous.rds", basename(arguments$bam)),
                "' saved to '",
                arguments$outdir, "'"
            )
        )
    } else {
        print("However, no ambiguous reads written to rds file.")
    }
}

#  Remove unneeded objects
rm(bam, chromosome, mate_status, mated, mated.rows_one, unmated)

print("Script completed.")


# #  Join the POS and MPOS bed files to the bam-file tibble ---------------------
# pos <- readr::read_tsv(
#     arguments$pos,
#     col_names = c(
#         "lO_rname", "lO_pos", "lO_pos_end", "qname", "lO_pos_status"
#     )
# )
# 
# mpos <- readr::read_tsv(
#     arguments$mpos,
#     col_names = c(
#         "lO_mrnm", "lO_mpos", "lO_mpos_end", "qname", "lO_mpos_status"
#     )
# )
# 
# pos_mpos <- dplyr::full_join(pos, mpos, by = "qname") %>%
#     dplyr::rename("qname.odd" = "qname")
# 
# #  Remove unneeded variables
# rm(pos, mpos)


#  Join the liftOver POS/MPOS information to the one-row-per-pair tibble ------
# #  Then, reorganize the tibble for easier readability 
# join <- dplyr::full_join(rows_one, pos_mpos, by = "qname.odd") %>%
#     dplyr::select(-c(lO_pos_end, lO_mpos_end)) %>%
#     dplyr::mutate(
#         lO_rname.even = lO_rname,
#         lO_pos.even = lO_pos,
#         lO_pos_status.even = lO_pos_status,
#         lO_mrnm.even = lO_mrnm,
#         lO_mpos.even = lO_mpos,
#         lO_mpos_status.even = lO_mpos_status
#     ) %>%
#     dplyr::rename(
#         lO_rname.odd = lO_rname,
#         lO_pos.odd = lO_pos,
#         lO_pos_status.odd = lO_pos_status,
#         lO_mrnm.odd = lO_mrnm,
#         lO_mpos.odd = lO_mpos,
#         lO_mpos_status.odd = lO_mpos_status
#     ) %>%
#     dplyr::relocate(
#         c(
#             lO_rname.odd, lO_pos.odd, lO_pos_status.odd,
#             lO_mrnm.odd, lO_mpos.odd, lO_mpos_status.odd,
#             rname.odd, pos.odd, mrnm.odd, mpos.odd, isize.odd
#         ),
#         .after = strand.odd
#     ) %>%
#     dplyr::relocate(
#         c(
#             lO_rname.even, lO_pos.even, lO_pos_status.even,
#             lO_mrnm.even, lO_mpos.even, lO_mpos_status.even,
#             rname.even, pos.even, mrnm.even, mpos.even, isize.even
#         ),
#         .after = strand.even
#     )


#  Remove unneeded variables
# rm(pos_mpos, rows_one)
# 
# operation <- paste0("rm(", bam, ")")
# evaluateOperation(operation)


# #  Retain pairs with successful liftOver status and matching chromosomes ------
# retain <- join %>%
#     subset(
#         lO_pos_status.odd == "liftOver successful" &
#         lO_mpos_status.odd == "liftOver successful" &
#         lO_pos_status.even == "liftOver successful" &
#         lO_mpos_status.even == "liftOver successful"
#     ) %>%
#     subset(
#         lO_rname.odd == chromosome &
#         lO_mrnm.odd == chromosome &
#         lO_rname.even == chromosome &
#         lO_mrnm.even == chromosome
#     )

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

# #  Remove unneeded variables, etc.
# rm(bam, chromosome, cl, command, operation, script)


#  Write out rds file for munged pairs information ----------------------------
# saveRDS(retain, file = paste0(arguments$outdir, "/", arguments$rds))
# saveRDS(rows_one, file = paste0(arguments$outdir, "/", arguments$rds))
# 
# print(
#     paste0("Script completed: outfile '", arguments$rds, "' saved to '", arguments$outdir, "'")
# )
