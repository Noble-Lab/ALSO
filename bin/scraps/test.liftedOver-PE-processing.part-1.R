#!/usr/bin/env Rscript

library(Rsamtools)
library(stringr)
library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)

#  Change to appropriate working directory
project <- "2021_kga0_4dn-mouse-cross"
default <- "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
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
            )
        )
    }
}
rm(project, default, current)


#  Set up functions -----------------------------------------------------------
`%notin%` <- Negate(`%in%`)


createTibbleFromTibbles <- function(vector_string_tibbles) {
    bind_rows(
        eval(parse(text = paste0(vector_string_tibbles)[1])),
        eval(parse(text = paste0(vector_string_tibbles)[2])),
        eval(parse(text = paste0(vector_string_tibbles)[3])),
        eval(parse(text = paste0(vector_string_tibbles)[4])),
        eval(parse(text = paste0(vector_string_tibbles)[5])),
        eval(parse(text = paste0(vector_string_tibbles)[6]))
    )
}


createVectorPlusMinus <- function(vector) {
    #  Takes a numeric vector and returns a vector of numbers ± one
    v1 <- vector
    v2 <- v1 - 1
    v3 <- v1 + 1
    vector_plus_minus <- c(v1, v2, v3) %>% sort() %>% unique()
    return(vector_plus_minus)
}


evaluateOperation <- function(operation = operation) {
    return(eval(parse(text = operation), envir = .GlobalEnv))
}


makeOperation <- function(variable = variable, command = command) {
    operation <- paste(variable, command) #%>% paste(., collapse = "; ")
    return(operation)
}


returnObjectName <- function(object) {
    return(deparse(substitute(object)))
}


testLogicalAlternating <- function(logical) {
    output <- vector(mode = "logical", length = length(logical))
    for (i in 1:length(logical)) {
        output[i] <- (logical[i] == logical[i + 1])
    }
    return(output)
}


testMatesPaired <- function(pos, mpos) {
    output <- vector(mode = "logical", length = length(pos))
    for (i in 1:length(pos)) {
        output[i] <- pos[i] == mpos[i + 1] & mpos[i] == pos[i + 1]
    }
    return(output)
}


testMposInPos <- function(pos, mpos) {
    mpos.in.pos <- mpos %in% pos
    return(mpos.in.pos)
}


testPosInMpos <- function(pos, mpos) {
    pos.in.mpos <- pos %in% mpos
    return(pos.in.mpos)
}


#  Set up work directory, etc. (locations TB∆) --------------------------------
directory_user <- "/Users/kalavattam"
directory_project <- "Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
directory_bam_base <- "Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do"
directory_bam <- paste0(directory_bam_base, "/", "get_unique_fragments/Bonora")
directory_data <- "data/kga0"
path.1 <- paste0(directory_user, "/", directory_project)
path.2 <- paste0(directory_bam)
path.3 <- paste0(directory_user, "/", directory_project, "/", directory_data)

rm(
    directory_bam, directory_bam_base,
    directory_data, directory_project, directory_user
)


#  Script name ----------------------------------------------------------------
script <- "test.liftedOver-PE-processing.part-1.R"


#  Read in sam2pairwise files in order to reconstruct CIGAR sequences ---------
#+
#+ Also, read in pertinent .bam files

#  Switch to location of standard .bam, sam2pairwise .bam, and sam2pairwise
#+ .txt files for KA .bam files
setwd(path.3)

#  Identify standard .bam files
chromosome <- "chrX"

file <- list.files(pattern = paste0(
    "\\", chromosome, ".rmdup.extendedCIGAR.bam$"
))[
    stringr::str_detect(
        list.files(pattern = paste0(
            "\\", chromosome, ".rmdup.extendedCIGAR.bam$"
        )),
        "mm10-CAST-129-Nmasked",
        negate = TRUE
    )
]
file <- file %>% stringr::str_subset("mm10\\.", negate = TRUE)
variable <- file %>%
    strsplit(., "-") %>%
    lapply(., `[[`, 1) %>% unlist() %>%
    paste0("standard.bam.", .)
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

#  Identify standard .bai files
file <- list.files(pattern = paste0(
    "\\", chromosome, ".rmdup.extendedCIGAR.bam.bai$"
))[
    stringr::str_detect(
        list.files(pattern = paste0(
            "\\", chromosome, ".rmdup.extendedCIGAR.bam.bai$"
        )),
        "mm10-CAST-129-Nmasked",
        negate = TRUE
    )
]
index <- file %>%
    strsplit(., "-") %>%
    lapply(., `[[`, 1) %>% unlist() %>%
    paste0("standard.index.", .)
mapply(
    assign, index, file, MoreArgs = list(envir = parent.frame())
)

#  Load in standard .bam fields
command <- paste0(
    "<- ", variable, " %>% ",
        "Rsamtools::BamFile(., index = ", index, ", asMates = TRUE)", " %>% ",
        "Rsamtools::scanBam()", " %>% ",
        "as.data.frame()", " %>% ",
        "tibble::as_tibble()"
)
operation <- makeOperation(paste0(variable, ".full"), command)
evaluateOperation(operation)

#  Load in standard .bam AS and MD fields
map_params <- Rsamtools::ScanBamParam(tag = c("AS", "MD"))
command <- paste0(
    "<- ", variable, " %>% ",
        "Rsamtools::BamFile(., index = ", index, ", asMates = TRUE)", " %>% ",
        "Rsamtools::scanBam(param = map_params)", " %>% ",
        "as.data.frame()", " %>% ",
        "tibble::as_tibble()"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)
rm(map_params)

#  Join the standard and AS fields
command <- paste0(
    "<- dplyr::bind_cols(", variable, ".full, ", variable, ")"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  Remove unneeded variables
operation <- paste0(
    "rm(", index, ", ", variable, ".full", ")"
)
evaluateOperation(operation)

variable.standard.bam <- variable

# -------------------------------------
#  Identify CrossMap .bam files
file <- list.files(pattern = paste0(
    "\\", chromosome, ".rmdup.extendedCIGAR.liftOverMm10.sorted.bam$"
))[
    stringr::str_detect(
        list.files(pattern = paste0(
            "\\", chromosome, ".rmdup.extendedCIGAR.liftOverMm10.sorted.bam$"
        )),
        "mm10-CAST-129-Nmasked",
        negate = TRUE
    )
]
variable <- file %>%
    strsplit(., "-") %>%
    lapply(., `[[`, 1) %>% unlist() %>%
    paste0("cm.bam.", .)
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

#  Identify CrossMap .bai files
file <- list.files(pattern = paste0(
    "\\", chromosome, ".rmdup.liftOverMm10.sorted.bam.bai$"
))[
    stringr::str_detect(
        list.files(pattern = paste0(
            "\\", chromosome, ".rmdup.liftOverMm10.sorted.bam.bai$"
        )),
        "mm10-CAST-129-Nmasked",
        negate = TRUE
    )
]
index <- file %>%
    strsplit(., "-") %>%
    lapply(., `[[`, 1) %>% unlist() %>%
    paste0("s2p.index.", .)
mapply(
    assign, index, file, MoreArgs = list(envir = parent.frame())
)

#  Load in CrossMap .bam fields
command <- paste0(
    "<- ", variable, " %>% ",
        "Rsamtools::BamFile(., index = ", index, ", asMates = TRUE)", " %>% ",
        "Rsamtools::scanBam()", " %>% ",
        "as.data.frame()", " %>% ",
        "tibble::as_tibble()"
)
operation <- makeOperation(paste0(variable, ".full"), command)
evaluateOperation(operation)

#  Load in CrossMap .bam AS and MD fields
map_params <- Rsamtools::ScanBamParam(tag = c("AS", "MD"))
command <- paste0(
    "<- ", variable, " %>% ",
        "Rsamtools::BamFile(., index = ", index, ", asMates = TRUE)", " %>% ",
        "Rsamtools::scanBam(param = map_params)", " %>% ",
        "as.data.frame()", " %>% ",
        "tibble::as_tibble()"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)
rm(map_params)

#  Join the standard and AS fields
command <- paste0(
    "<- dplyr::bind_cols(", variable, ".full, ", variable, ")"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  Remove unneeded variables
operation <- paste0(
    "rm(", index, ", ", variable, ".full", ")"
)
evaluateOperation(operation)

variable.cm.bam <- variable

# -------------------------------------
#  Identify sam2pairwise .txt files for CrossMap-lifted-over .bam files
file <- list.files(pattern = paste0(
    "\\.rmdup.liftOverMm10.sorted.sam2pairwise.munged.txt$"
))[
    stringr::str_detect(
        list.files(pattern = paste0(
            "\\.rmdup.liftOverMm10.sorted.sam2pairwise.munged.txt"
        )),
        "mm10-CAST-129-Nmasked",
        negate = TRUE
    )
]
variable <- file %>%
    strsplit(., "-") %>%
    lapply(., `[[`, 1) %>% unlist() %>%
    paste0("s2p.txt.", .)
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

#  Load in sam2pairwise .txt files for CrossMap-lifted-over .bam files
command <- paste0("<- readr::read_tsv(", variable, ")")
operation <- makeOperation(variable, command)
evaluateOperation(operation)

variable.s2p.txt <- variable

#  Remove variables no longer needed
rm(chromosome, file, index)


#  Sort, munge, and join files ------------------------------------------------
variable <- c(variable.cm.bam, variable.standard.bam)
command <- paste0(
    "<- ", variable, "[order(", variable, "$qname, ", variable, "$seq), ]"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

variable <- variable.s2p.txt
command <- paste0(
    "<- ", variable, "[order(", variable, "$qname, ", variable, "$read_sequence), ]"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  Set up $coordinate, a variable needed for sorting too: qname, rname, pos
variable <- c(variable.cm.bam, variable.standard.bam, variable.s2p.txt)
command <- paste0(
    "<- paste0(",
            variable, "$qname, ", "'_', ",
            variable, "$rname, ", "'_', ",
            variable, "$pos",
    ")"
)
operation <- makeOperation(paste0(variable, "$coordinate"), command)
evaluateOperation(operation)

#  Set up $qpos, a variable needed for sorting
command <- paste0(
    "<- paste0(",
        variable, "$qname, ", "'_', ",
        variable, "$pos",
    ")"
)
operation <- makeOperation(paste0(variable, "$qpos"), command)
evaluateOperation(operation)

#  Set up $qmpos, a variable needed for sorting
command <- paste0(
    "<- paste0(",
        variable, "$qname, ", "'_', ",
        variable, "$mpos",
    ")"
)
operation <- makeOperation(paste0(variable, "$qmpos"), command)
evaluateOperation(operation)


#  Narrow down pertinent fields for column binding
variable <- c(variable.cm.bam, variable.standard.bam)
command <- paste0(
    "<- ", variable, " %>% ",
        "dplyr::select(",
            "coordinate, flag, rname, pos, mapq, cigar, mrnm, mpos, isize, seq",
        ")"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

variable <- variable.s2p.txt
command <- paste0(
    "<- ", variable, " %>% ",
        "dplyr::select(",
            "coordinate, flag, rname, pos, mapq, cigar, mrnm, mpos, isize, read_sequence, matches, reference_sequence",
        ")"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  Upon surveying the files and fields, column-bind the standard.bam.* and
#+ s2p.txt.* tibbles (no need to include information from the cm.bam.* tibbles
#+ because the information is already included in the s2p.txt.* tibbles)

#  First, rename columns
command <- paste0(
    "<- ", variable.standard.bam, " %>% ",
        "dplyr::rename(coordinate.standard = coordinate)", " %>% ",
        "dplyr::rename(flag.standard = flag)", " %>% ",
        "dplyr::rename(rname.standard = rname)", " %>% ",
        "dplyr::rename(pos.standard = pos)", " %>% ",
        "dplyr::rename(mapq.standard = mapq)", " %>% ",
        "dplyr::rename(cigar.standard = cigar)", " %>% ",
        "dplyr::rename(mrnm.standard = mrnm)", " %>% ",
        "dplyr::rename(mpos.standard = mpos)", " %>% ",
        "dplyr::rename(isize.standard = isize)", " %>% ",
        "dplyr::rename(seq.standard = seq)"
)
operation <- makeOperation(variable.standard.bam, command)
evaluateOperation(operation)

command <- paste0(
    "<- ", variable.s2p.txt, " %>% ",
        "dplyr::rename(coordinate.sam2pairwise = coordinate)", " %>% ",
        "dplyr::rename(flag.sam2pairwise = flag)", " %>% ",
        "dplyr::rename(rname.sam2pairwise = rname)", " %>% ",
        "dplyr::rename(pos.sam2pairwise = pos)", " %>% ",
        "dplyr::rename(mapq.sam2pairwise = mapq)", " %>% ",
        "dplyr::rename(cigar.sam2pairwise = cigar)", " %>% ",
        "dplyr::rename(mrnm.sam2pairwise = mrnm)", " %>% ",
        "dplyr::rename(mpos.sam2pairwise = mpos)", " %>% ",
        "dplyr::rename(isize.sam2pairwise = isize)", " %>% ",
        "dplyr::rename(read_sequence.sam2pairwise = read_sequence)", " %>% ",
        "dplyr::rename(matches.sam2pairwise = matches)", " %>% ",
        "dplyr::rename(reference_sequence.sam2pairwise = reference_sequence)"
)
operation <- makeOperation(variable.s2p.txt, command)
evaluateOperation(operation)

#  Bind the columns
variable <- c("bind.129S1", "bind.CAST")
command <- paste0(
    "<- ", "dplyr::bind_cols(",
        variable.standard.bam, ", ", variable.s2p.txt,
    ")"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Remove unneeded variables, clean up the environment ------------------------
operation <- paste0(
    "rm(",
        variable.cm.bam, ", ", variable.s2p.txt, ", ", variable.standard.bam,
    ")"
)
evaluateOperation(operation)

rm(variable, variable.cm.bam, variable.s2p.txt, variable.standard.bam)


#  Save an image, end the script ----------------------------------------------
save.image(stringr::str_remove(script, ".R") %>% paste0(., ".Rdata"))
rm(list = ls())

#  Run 'copy-important-files.sh', then go to 'test.PE-processing.part-7.R'
