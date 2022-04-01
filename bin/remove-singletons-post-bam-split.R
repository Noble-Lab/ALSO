#!/usr/bin/env Rscript

#  remove-singletons-post-bam-split.R
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


#  Source libraries, adjust settings ------------------------------------------
importLibrary(c("argparser", "Rsamtools", "tidyverse"))

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)


#  Parse arguments ------------------------------------------------------------
#  Create a parser
ap <- arg_parser(
    name = "remove-singletons-post-bam-split.R",
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
    short = "-d",
    arg = "--header",
    type = "character",
    default = NULL,
    help = "bam header, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-o",
    arg = "--outdir",
    type = "character",
    default = NULL,
    help = "directory for saving rds outfile, including path <chr>"
)

#  Parse the command line arguments
dir_base <- "."
dir_data <- "data"
dir_in <- paste0(
    dir_base, "/",
    dir_data, "/",
    "2022-0330_troubleshoot_pipeline-singletons-issue"
)
dir_out <- paste0(
    dir_base, "/",
    dir_data, "/",
    "2022-0330_troubleshoot_pipeline-singletons-issue"
)
bam <- paste0(
    dir_in, "/",
    "Disteche_sample_6.dedup.process.sort_name.fixmate.sort_coord.chr11.bam"
)
bai <- paste0(
    dir_in, "/",
    "Disteche_sample_6.dedup.process.sort_name.fixmate.sort_coord.chr11.bam.bai"
)
header <- paste0(
    dir_in, "/",
    "Disteche_sample_6.dedup.process.sort_name.fixmate.sort_coord.header.txt"
)

cl <- c(
    "--bam", bam,
    "--bai", bai,
    "--header", header,
    "--outdir", dir_out
)
arguments <- parse_args(ap, cl)  # RStudio-interactive work
# arguments <- parse_args(ap)  # Command-line calls

rm(dir_base, dir_data, dir_in, dir_out, bam, bai, header, ap, cl)


#  Check that files exist -----------------------------------------------------
stopifnot(file.exists(arguments$bam))
stopifnot(file.exists(arguments$bai))
stopifnot(file.exists(arguments$header))

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
map_params <- Rsamtools::ScanBamParam(
    # tag = c(
    #     "AS", "XS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT", "MQ", "MC"
    # )
    tag = c("AS", "MD")
)
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

#  Remove reads with $mate_status of "unmated"
command <- paste0(
    bam, "[", bam, "$mate_status != \"unmated\", ]"
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

#  Add string to $tag.AS 
command <- paste0(
    "paste0(\"AS:i:\", ", bam, "$tag.AS)"
)
operation <- makeOperation(paste0(bam, "$tag.AS"), command)
evaluateOperation(operation)

# #  Add string to $tag.XS
# command <- paste0(
#     "paste0(\"XS:i:\", ", bam, "$tag.XS)"
# )
# operation <- makeOperation(paste0(bam, "$tag.XS"), command)
# evaluateOperation(operation)
# 
# #  Add string to $tag.XN
# command <- paste0(
#     "paste0(\"XN:i:\", ", bam, "$tag.XN)"
# )
# operation <- makeOperation(paste0(bam, "$tag.XN"), command)
# evaluateOperation(operation)
# 
# #  Add string to $tag.XM
# command <- paste0(
#     "paste0(\"XM:i:\", ", bam, "$tag.XM)"
# )
# operation <- makeOperation(paste0(bam, "$tag.XM"), command)
# evaluateOperation(operation)
# 
# #  Add string to $tag.XO
# command <- paste0(
#     "paste0(\"XO:i:\", ", bam, "$tag.XO)"
# )
# operation <- makeOperation(paste0(bam, "$tag.XO"), command)
# evaluateOperation(operation)
# 
# #  Add string to $tag.XG
# command <- paste0(
#     "paste0(\"XG:i:\", ", bam, "$tag.XG)"
# )
# operation <- makeOperation(paste0(bam, "$tag.XG"), command)
# evaluateOperation(operation)
# 
# #  Add string to $tag.NM
# command <- paste0(
#     "paste0(\"NM:i:\", ", bam, "$tag.NM)"
# )
# operation <- makeOperation(paste0(bam, "$tag.NM"), command)
# evaluateOperation(operation)

#  Add string to $tag.MD
command <- paste0(
    "paste0(\"MD:Z:\", ", bam, "$tag.MD)"
)
operation <- makeOperation(paste0(bam, "$tag.MD"), command)
evaluateOperation(operation)

# #  Add string to $tag.YS
# command <- paste0(
#     "paste0(\"YS:i:\", ", bam, "$tag.YS)"
# )
# operation <- makeOperation(paste0(bam, "$tag.YS"), command)
# evaluateOperation(operation)
# 
# #  Add string to $tag.YT
# command <- paste0(
#     "paste0(\"YT:Z:\", ", bam, "$tag.YT)"
# )
# operation <- makeOperation(paste0(bam, "$tag.YT"), command)
# evaluateOperation(operation)
# 
# #  Add string to $tag.MQ
# command <- paste0(
#     "paste0(\"MQ:i:\", ", bam, "$tag.MQ)"
# )
# operation <- makeOperation(paste0(bam, "$tag.MQ"), command)
# evaluateOperation(operation)
# 
# #  Add string to $tag.MC
# command <- paste0(
#     "paste0(\"MC:Z:\", ", bam, "$tag.MC)"
# )
# operation <- makeOperation(paste0(bam, "$tag.MC"), command)
# evaluateOperation(operation)


#  Write out bam file ---------------------------------------------------------
#  Write out sam file
readr::write_tsv(
    x = eval(parse(text = bam)),
    file = paste0(
        arguments$outdir, "/", stringr::str_remove(bam, ".bam"), ".tmp.sam"
    ),
    col_names = FALSE
)

#  Add header to sam file
shell_code <- paste0(
    "cat ", arguments$outdir, "/", stringr::str_remove(bam, ".bam"), ".tmp.sam",
    " >> ", arguments$header,
    " && ",
    "mv -f ", arguments$header, " ",
    arguments$outdir, "/", stringr::str_remove(bam, ".bam"), ".tmp.sam"
)
system(shell_code)

#  Convert sam file to bam file
shell_code <- paste0(
    "samtools view -h -b ",
    paste0(arguments$outdir, "/", stringr::str_remove(bam, ".bam"), ".tmp.sam"), " > ",
    paste0(arguments$outdir, "/", stringr::str_remove(bam, ".bam"), ".tmp.bam")
)
system(shell_code)
