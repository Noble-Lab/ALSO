#!/usr/bin/env Rscript

#  KA

#  Set up working directory ---------------------------------------------------
project <- "2021_kga0_4dn-mouse-cross"
default <- paste0("/Users/kalavattam/Dropbox/UW/projects-etc/", project)
current <- stringr::str_split(getwd(), "/")[[1]][
    length(stringr::str_split(getwd(), "/")[[1]])
]
if(current != project) {
    if(dir.exists(default)) {
        setwd(default)
    } else {
        setwd(
            readline(
                prompt = paste0(
                    "Enter path to and including the project directory, ",
                    project, ":"
                )
            )  #TODO Check again: current == project?
        )
    }
}
rm(project, default, current)


#  Source libraries and functions ---------------------------------------------
source("bin/auxiliary/auxiliary.R")
source("bin/auxiliary/auxiliary_experiment-metadata.R")

packages <- c("parallel", "Rsamtools", "stringr", "tidyverse")
importLibrary(packages)
rm(packages, preload)

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)


`%notin%` <- Negate(`%in%`)


assignFilesToVariables <- function(file, string) {
    split <- file %>% strsplit(., "/")
    name <- split %>% 
        lapply(., `[[`, length(split[[1]])) #%>%
        # unlist() %>%
        # strsplit(., "-") %>%
        # lapply(., `[[`, 1) %>%
        # unlist() %>%
        # paste0(string, .) %>%
        # stringr::str_remove("S1")
    return(name)
}


evaluateOperation <- function(operation = operation) {
    return(eval(parse(text = operation), envir = .GlobalEnv))
}


makeOperation <- function(variable = variable, command = command) {
    operation <- paste0(variable, " <- ", command)
    return(operation)
}


#  Parse arguments ------------------------------------------------------------
script <- "TBD.R"

#  Create a parser
ap <- arg_parser(
    name = script,
    description = "",
    hide.opts = TRUE
)

#  Add command line arguments
ap <- callMetadataArguments()
ap <- add_argument(
    ap,
    short = "-i",
    arg = "--infile",
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
    short = "-o",
    arg = "--outdir",
    type = "character",
    default = NULL,
    help = "directory for saving outfile(s), including path <chr>"
)

#  Parse the command line arguments
directory_base <- "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
directory_data <- "data/kga0"
directory_out <- paste0(directory_base, "/", directory_data)
infile <- paste0(
    directory_out, "/",
    "Disteche_sample_1.dedup.bam"
)
bai <- paste0(
    directory_out, "/",
    "Disteche_sample_1.dedup.bam.bai"
)
cl <- c(
    #  Arguments for analysis
    "--infile", infile,
    "--bai", bai,
    "--outdir", directory_out,

    #  Metadata arguments
    "--username", "kga0",
    "--experiment_date", "2022-0313",
    "--experiment_description", "TBD"
)
arguments <- parse_args(ap, cl)  # RStudio-interactive work
# arguments <- parse_args(ap)  # Command-line calls

rm(directory_base, directory_data, directory_out, infile, bai)


#  Load in .bam information, including mate information -----------------------
file <- c(arguments$infile)
variable <- assignFilesToVariables(file)
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

file <- c(arguments$bai)
index <- assignFilesToVariables(file)
mapply(
    assign, index, file, MoreArgs = list(envir = parent.frame())
)

#  Using Rsamtools, load in standard .bam fields
command <- paste0(
    "Rsamtools::BamFile(",
        variable, ", index = ", index, ", asMates = TRUE",
    ")", " %>% ",
    "Rsamtools::scanBam()", " %>% ",
    "as.data.frame()", " %>% ",
    "tibble::as_tibble()"  #TODO Speed up via parallelization?
)
operation <- makeOperation("full", command)
evaluateOperation(operation)

#  Using Rsamtools, load in .bam AS and MD fields
map_params <- Rsamtools::ScanBamParam(tag = c("AS", "MD"))
command <- paste0(
    variable, " %>% ",
        "Rsamtools::BamFile(",
            variable, ", index = ", index, ", asMates = TRUE",
        ")", " %>% ",
        "Rsamtools::scanBam(param = map_params)", " %>% ",
        "as.data.frame()", " %>% ",
        "tibble::as_tibble()"  #TODO Speed up via parallelization?
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

rm(map_params)

#  Column-bind the standard, AS, and MD fields
command <- paste0(
    "dplyr::bind_cols(", variable, ".full, ", variable, ")"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  Remove unneeded variables
operation <- paste0("rm(", variable, ".full, ", index, ")")
evaluateOperation(operation)

rm(file, index)


#  Split tibbles based on $mate_status ----------------------------------------
split <- c("ambiguous", "mated", "unmated")
for (i in 1:length(split)) {
    variable_split <- paste0(split[i], ".", variable)
    command <- paste0(
        variable, "[", variable, "$mate_status == '", split[i], "', ]"
    )
    operation <- makeOperation(variable_split, command)
    evaluateOperation(operation)
}
rm(i)

#  Remove unneeded initial variables
command <- paste0("rm(", variable, ")")
evaluateOperation(command)

#  Create a vector of names of "split variables"
split <- paste0(split, ".")
variable <- purrr::cross(list(split, variable)) %>%
    purrr::map_chr(paste0, collapse = "")

#  Clean up variables
rm(split, variable_split)
