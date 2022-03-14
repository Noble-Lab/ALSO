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
        lapply(., `[[`, length(split[[1]])) %>%
        unlist() %>%
        strsplit(., "-") %>%
        lapply(., `[[`, 1) %>%
        unlist() %>%
        paste0(string, .) %>%
        stringr::str_remove("S1")
    return(name)
}


bindRows <- function(tbl) {
    dplyr::bind_rows(
        evaluateOperation(tbl[1]),
        evaluateOperation(tbl[2]),
    )
}


createTibbleFromTibbles <- function(vector_string_tibbles) {
    bind_rows(
        eval(parse(text = paste0(vector_string_tibbles)[1])),
        eval(parse(text = paste0(vector_string_tibbles)[2])),
        eval(parse(text = paste0(vector_string_tibbles)[3])),
        eval(parse(text = paste0(vector_string_tibbles)[4])),
        eval(parse(text = paste0(vector_string_tibbles)[5])),
        eval(parse(text = paste0(vector_string_tibbles)[6])),
        eval(parse(text = paste0(vector_string_tibbles)[7])),
        eval(parse(text = paste0(vector_string_tibbles)[8])),
        eval(parse(text = paste0(vector_string_tibbles)[9]))
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
    operation <- paste0(variable, " <- ", command)
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


#  Parse arguments ------------------------------------------------------------
script <- "script-1.R"

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
    short = "-s1",
    arg = "--strain_1",
    type = "character",
    default = NULL,
    help = "full path and .bam from alignment to strain 1, e.g., 129S1-SvImJ <chr>"
)
ap <- add_argument(
    ap,
    short = "-s2",
    arg = "--strain_2",
    type = "character",
    default = NULL,
    help = "full path and .bam from alignment to strain 2, e.g., CAST-EiJ <chr>"
)
ap <- add_argument(
    ap,
    short = "-st",
    arg = "--strain_test",
    type = "character",
    default = NULL,
    help = "full path and .bam from alignment to strain test, e.g., N-masked mm10 <chr>"
)
ap <- add_argument(
    ap,
    short = "-b1",
    arg = "--bai_1",
    type = "character",
    default = NULL,
    help = "full path and .bai from alignment to strain 1, 129S1-SvImJ <chr>"
)
ap <- add_argument(
    ap,
    short = "-b2",
    arg = "--bai_2",
    type = "character",
    default = NULL,
    help = "full path and .bai from alignment to strain 2, e.g., CAST-EiJ <chr>"
)
ap <- add_argument(
    ap,
    short = "-bt",
    arg = "--bai_test",
    type = "character",
    default = NULL,
    help = "full path and .bai from alignment to strain test, e.g., N-masked mm10 <chr>"
)
ap <- add_argument(
    ap,
    short = "-o",
    arg = "--directory_out",
    type = "character",
    default = NULL,
    help = "full path for saving outfile(s) <chr>"
)

#  Parse the command line arguments
directory_base <- "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
directory_data <- "data/kga0"
directory_out <- paste0(directory_base, "/", directory_data)
s_1 <- paste0(
    directory_out, "/",
    "129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam"
)
s_2 <- paste0(
    directory_out, "/",
    "CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam"
)
s_test <- paste0(
    directory_out, "/",
    "mm10-CAST-129-Nmasked.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam"
)
b_1 <- paste0(
    directory_out, "/",
    "129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam.bai"
)
b_2 <- paste0(
    directory_out, "/",
    "CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam.bai"
)
b_test <- paste0(
    directory_out, "/",
    "mm10-CAST-129-Nmasked.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam.bai"
)
cl <- c(
    #  Arguments for analysis
    "--strain_1", s_1,
    "--strain_2", s_2,
    "--strain_test", s_test,
    "--bai_1", b_1,
    "--bai_2", b_2,
    "--bai_test", b_test,
    "--directory_out", directory_out,

    #  Metadata arguments
    "--username", "kga0",
    "--experiment_date", "2022-0219",
    "--experiment_description", "script-1"
)
arguments <- parse_args(ap, cl)  # RStudio-interactive work
# arguments <- parse_args(ap)  # Command-line calls

rm(
    b_1, b_2, b_test,
    directory_base, directory_data, directory_out,
    s_1, s_2, s_test
)


#  Load in .bam information, including mate information -----------------------
file <- c(arguments$strain_1, arguments$strain_2, arguments$strain_test)
variable <- assignFilesToVariables(file, "dedup.")
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

file <- c(arguments$bai_1, arguments$bai_2, arguments$bai_test)
index <- assignFilesToVariables(file, "index.")
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
operation <- makeOperation(paste0(variable, ".full"), command)
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


#  Drop rows that are not chr1-19, chrX ---------------------------------------
chromosomes <- c(paste0("chr", 1:19), "chrX")
command <- paste0(
    variable, " %>% ", "filter(., rname %in% chromosomes)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Reorder rname factor levels ------------------------------------------------
command <- paste0(
    "forcats::fct_relevel(", variable, "$rname)"
)
operation <- makeOperation(paste0(variable, "$rname"), command)
evaluateOperation(operation)


#  Drop unused rname factor levels --------------------------------------------
command <- paste0(
    variable, "$rname %>% forcats::fct_drop()"
)
operation <- makeOperation(paste0(variable, "$rname"), command)
evaluateOperation(operation)


#  Create and append $pos_end, $mpos_end --------------------------------------
command <- paste0(variable, "$pos + 49")
operation <- makeOperation(paste0(variable, "$pos_end"), command)
evaluateOperation(operation)

command <- paste0(variable, "$mpos + 49")
operation <- makeOperation(paste0(variable, "$mpos_end"), command)
evaluateOperation(operation)

#  Move *_end to appropriate locations
command <- paste0(
    variable, " %>% ",
        "dplyr::relocate(pos_end, .after = pos) %>% ",
        "dplyr::relocate(mpos_end, .after = mpos)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Filter out rows with mapq less than 30 -------------------------------------
command <- paste0(variable, "[", variable, "$mapq >= 30, ]")
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Prepare for tibble-sorting; deduplicate and mate-label the tibbles ---------

#  Set up $criteria, a variable needed for sorting: qname, pos, mpos
command <- paste0(
    "paste0(",
        variable, "$qname, ", "'_', ",
        variable, "$pos, ", "'_', ",
        variable, "$mpos",
    ")"
)
operation <- makeOperation(paste0(variable, "$criteria"), command)
evaluateOperation(operation)

# #  Set up $criteria_flag, a variable needed for sorting: qname, flag, pos, mpos
# command <- paste0(
#     "paste0(",
#         variable, "$qname, ", "'_', ",
#         variable, "$flag, ", "'_', ",
#         variable, "$pos, ", "'_', ",
#         variable, "$mpos",
#     ")"
# )
# operation <- makeOperation(paste0(variable, "$criteria_flag"), command)
# evaluateOperation(operation)

#  Set up $qpos, a variable needed for sorting: qname, pos
command <- paste0(
    "paste0(",
        variable, "$qname, ", "'_', ",
        variable, "$pos",
    ")"
)
operation <- makeOperation(paste0(variable, "$qpos"), command)
evaluateOperation(operation)

#  Set up $qmpos, a variable needed for sorting: qname, mpos
command <- paste0(
    "paste0(",
        variable, "$qname, ", "'_', ",
        variable, "$mpos",
    ")"
)
operation <- makeOperation(paste0(variable, "$qmpos"), command)
evaluateOperation(operation)

#  Set up $isize_abs, i.e., absolute insert-size values, for sorting
command <- paste0(
    variable, "$isize %>% abs() %>% as.numeric()"
)
operation <- makeOperation(paste0(variable, "$isize_abs"), command)
evaluateOperation(operation)

#  To keep the original munged tibbles together as new variable populate the
#+ environment, rename the tibbles associated with 'variable'
command <- paste0(variable)
operation <- makeOperation(paste0("a.", variable), command)
evaluateOperation(operation)

#  Remove unneeded variables
command <- paste0("rm(", variable, ")")
evaluateOperation(command)

#  Adjust values associated 'variable'
variable <- paste0("a.", variable)

#  Rearrange the order of tibble columns for quick, easy surveying
command <- paste0(
    variable, " %>% ",
        "dplyr::relocate(mpos, .after = pos) %>% ",
        "dplyr::relocate(mpos_end, .after = pos_end) %>% ",
        "dplyr::relocate(mrnm, .after = rname) %>% ",
        "dplyr::relocate(isize, .after = mpos_end) %>% ",
        "dplyr::relocate(criteria, .before = groupid) %>% ",
        # "dplyr::relocate(criteria_flag, .before = groupid) %>% ",
        "dplyr::relocate(tag.MD, .after = cigar) %>% ",
        "dplyr::relocate(tag.AS, .after = tag.MD) %>% ",
        "dplyr::relocate(qpos, .after = criteria) %>% ",
        "dplyr::relocate(qmpos, .after = qpos) %>% ",
        "dplyr::relocate(qname, .after = qual) %>% ",
        "dplyr::relocate(isize_abs, .before = isize)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  To survey duplicates, split out entries in which $qpos is the same as
#+ $qmpos, i.e., entries in which each member of the pair maps to the same
#+ location
command <- paste0(
    variable, "[", variable, "$qpos == ", variable, "$qmpos, ]"
)
operation <- makeOperation(paste0("same.", variable), command)
evaluateOperation(operation)

#  Remove unneeded variables
command <- paste0("rm(", "same.", variable, ")")
evaluateOperation(command)

#  Regarding the main tibbles, remove entries in which $qpos is the same as
#+ $qmpos, i.e., entries in which each member of the pair maps to the same
#+ location
command <- paste0(
    variable, "[!(", variable, "$qpos == ", variable, "$qmpos), ]"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  Deduplicate the main tibbles based on $criteria: There should be no more
#+ than one entry with a given $criteria value (combination of $qname,
#+ $pos, and $mpos)
command <- paste0(
    variable, "[!duplicated(", variable, "$criteria), ]"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Is $qpos in $qmpos? Vice versa? Create tables and tibbles of answers -------
command <- paste0(
    "testPosInMpos(",
        "pos = ", variable, "$qpos, ",
        "mpos = ", variable, "$qmpos",
    ") %>% ",
        "table()"
)
operation <- makeOperation(paste0("PIM.", variable), command)
evaluateOperation(operation)

command <- paste0(
    "testMposInPos(",
        "pos = ", variable, "$qpos, ",
        "mpos = ", variable, "$qmpos",
    ") %>% ",
        "table()"
)
operation <- makeOperation(paste0("MIP.", variable), command)
evaluateOperation(operation)

#  Convert the tables to tibbles: PIM
command <- paste0(
    "PIM.", variable, " %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("PIM.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles:" PIM
command <- paste0(
    "tibble::tibble(", "PIM.", variable, ")"
)
operation <- makeOperation(paste0("PIM.", variable), command)
evaluateOperation(operation)

vector_string_tibbles <- paste0("PIM.", variable)
z.check.1.PIM <- createTibbleFromTibbles(vector_string_tibbles)
z.check.1.PIM$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.1.PIM <- z.check.1.PIM %>% dplyr::arrange(rownames)

#  Convert the tables to tibbles: MIP
command <- paste0(
    "MIP.", variable, " %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("MIP.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles:" MIP
command <- paste0(
    "tibble::tibble(", "MIP.", variable, ")"
)
operation <- makeOperation(paste0("MIP.", variable), command)
evaluateOperation(operation)

vector_string_tibbles <- paste0("MIP.", variable)
z.check.1.MIP <- createTibbleFromTibbles(vector_string_tibbles)
z.check.1.MIP$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.1.MIP <- z.check.1.MIP %>% dplyr::arrange(rownames)


#  Does PIM.*$`FALSE` == MIP.*$`FALSE`? ---------------------------------------
command <- paste0(
    "PIM.", variable, "$`FALSE` == MIP.", variable, "$`FALSE`"
)
operation <- makeOperation(paste0("PeqM.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles"
command <- paste0(
    "tibble::tibble(", "PeqM.", variable, ")"
)
operation <- makeOperation(paste0("PeqM.", variable), command)
evaluateOperation(operation)

operation <- paste0("colnames(", "PeqM.", variable, ") ", "<- ", "\"n\"")
evaluateOperation(operation)

vector_string_tibbles <- paste0("PeqM.", variable)
z.check.2.PeqM <- createTibbleFromTibbles(vector_string_tibbles)
z.check.2.PeqM$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.2.PeqM <- z.check.2.PeqM %>% dplyr::arrange(rownames)


#  Does PIM.*$`FALSE` > MIP.*$`FALSE`? ----------------------------------------
command <- paste0(
    "PIM.", variable, "$`FALSE` > MIP.", variable, "$`FALSE`"
)
operation <- makeOperation(paste0("PgtM.", variable), command)
evaluateOperation(operation)

#  Create a tibble
command <- paste0(
    "tibble::tibble(", "PgtM.", variable, ")"
)
operation <- makeOperation(paste0("PgtM.", variable), command)
evaluateOperation(operation)

operation <- paste0("colnames(", "PgtM.", variable, ") ", "<- ", "\"n\"")
evaluateOperation(operation)

vector_string_tibbles <- paste0("PgtM.", variable)
z.check.2.PgtM <- createTibbleFromTibbles(vector_string_tibbles)
z.check.2.PgtM$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.2.PgtM <- z.check.2.PgtM %>% dplyr::arrange(rownames)


#  Is PIM.*$`FALSE` < MIP.*$`FALSE`? ------------------------------------------
command <- paste0(
    "PIM.", variable, "$`FALSE` < MIP.", variable, "$`FALSE`"
)
operation <- makeOperation(paste0("PltM.", variable), command)
evaluateOperation(operation)

#  Create a tibble
command <- paste0(
    "tibble::tibble(", "PltM.", variable, ")"
)
operation <- makeOperation(paste0("PltM.", variable), command)
evaluateOperation(operation)

operation <- paste0("colnames(", "PltM.", variable, ") ", "<- ", "\"n\"")
evaluateOperation(operation)

vector_string_tibbles <- paste0("PltM.", variable)
z.check.2.PltM <- createTibbleFromTibbles(vector_string_tibbles)
z.check.2.PltM$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.2.PltM <- z.check.2.PltM %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0(
    "rm(",
        "PIM.", variable, ", ",
        "MIP.", variable, ", ",
        "PeqM.", variable, ", ",
        "PltM.", variable, ", ",
        "PgtM.", variable,
    ")"
)
evaluateOperation(command)
rm(vector_string_tibbles)


#  Make logical vectors for $qpos in $qmpos and vice versa --------------------

#  PIM
command <- paste0(
    "testPosInMpos(pos = ", variable, "$qpos, mpos = ", variable, "$qmpos)"
)
operation <- makeOperation(paste0("vPIM.", variable), command)
evaluateOperation(operation)

# #  Check: NAs in vPIM vector?
# command <- paste0("is.na(", "vPIM.", variable, ")")
# operation <- makeOperation(paste0("isNA.vPIM.", variable), command)
# evaluateOperation(operation)
# 
# command <- paste0(
#     paste0("isNA.vPIM.", variable, " %>% ", "table()")
# )
# operation <- makeOperation(paste0("isNA.vPIM.", variable), command)
# evaluateOperation(operation)  # No NAs
# 
# command <- paste0("rm(", "isNA.vPIM.", variable, ")")
# evaluateOperation(command)

#  MIP
command <- paste0(
    "testMposInPos(pos = ", variable, "$qpos, mpos = ", variable, "$qmpos)"
)
operation <- makeOperation(paste0("vMIP.", variable), command)
evaluateOperation(operation)

# #  Check: NAs in vMIP vector?
# command <- paste0("is.na(", "vMIP.", variable, ")")
# operation <- makeOperation(paste0("isNA.vMIP.", variable), command)
# evaluateOperation(operation)
# 
# command <- paste0(
#     paste0("isNA.vMIP.", variable, " %>% ", "table()")
# )
# operation <- makeOperation(paste0("isNA.vMIP.", variable), command)
# evaluateOperation(operation)  # No NAs
# 
# command <- paste0("rm(", "isNA.vMIP.", variable, ")")
# evaluateOperation(command)

#NOTE
#  Based on PeqM.*, PgtM.*, and PltM.* results, use the following vectors:
#+ All vPIM.* are fine, apparently

#  Remove unneeded variables
command <- paste0("rm(", "vMIP.", variable, ")")
evaluateOperation(command)


#  Subset tibbles based on logical vectors for $qpos in $qmpos ----------------

#NOTE vPIM: logical vectors for $qpos %in% $qmpos
command <- paste0(variable, "[", "vPIM.", variable, ", ]")
operation <- makeOperation(paste0("sub.", variable), command)
evaluateOperation(operation)

#  Remove unneeded vPIM vectors
command <- paste0("rm(", "vPIM.", variable, ")")
evaluateOperation(command)


#  For the subsetted tibbles, save qpos %in$ qmpos, vice versa in vectors -----
variable <- paste0("sub.", variable)
command <- paste0(  # Table qpos %in$ qmpos
    "testPosInMpos(",
        "pos = ", variable, "$qpos, ",
        "mpos = ", variable, "$qmpos",
    ") %>% ",
        "table()"
)
operation <- makeOperation(paste0("PIM.", variable), command)
evaluateOperation(operation)

command <- paste0(  # Table qmpos %in$ qpos
    "testMposInPos(",
        "pos = ", variable, "$qpos, ",
        "mpos = ", variable, "$qmpos",
    ") %>% ",
        "table()"
)
operation <- makeOperation(paste0("MIP.", variable), command)
evaluateOperation(operation)

#  Convert the tables to tibbles
command <- paste0(
    "PIM.", variable, " %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("PIM.", variable), command)
evaluateOperation(operation)

command <- paste0(
    "MIP.", variable, " %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("MIP.", variable), command)
evaluateOperation(operation)

#  Collect pertinent readouts for the subsetted tibbles in one tibble; i.e.,
#+ create a "tibble of tibbles"

#  PIM
vector_string_tibbles <- paste0("PIM.", variable)
z.check.sub.0.PIM <- createTibbleFromTibbles(vector_string_tibbles)
z.check.sub.0.PIM$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.sub.0.PIM <- z.check.sub.0.PIM %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0("rm(", "PIM.", variable, ")")
evaluateOperation(command)
rm(vector_string_tibbles)

#  MIP
vector_string_tibbles <- paste0("MIP.", variable)
z.check.sub.0.MIP <- createTibbleFromTibbles(vector_string_tibbles)
z.check.sub.0.MIP$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.sub.0.MIP <- z.check.sub.0.MIP %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0("rm(", "MIP.", variable, ")")
evaluateOperation(command)
rm(vector_string_tibbles)


#  Order the tibbles based on $isize_abs, $qpos, and $qmpos -------------------
command <- paste0(
    variable, "[",
        "order(",
            variable, "$isize_abs, ",
            variable, "$qpos, ",
            variable, "$qmpos",
        "), ",
    "]"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Clean up global environment, save environment image (part 1/2) -------------
operation <- paste0(
    "rm(", stringr::str_remove(variable, "sub."), ")"
)
evaluateOperation(operation)

save.image(paste0(
        arguments$directory_out, "/",
        stringr::str_remove(script, ".R"), "_part-1", ".Rdata"
))


#  Add row numbers to sub.* tibbles -------------------------------------------
operation <- paste0(
    variable, "$row_n", " <- ", "seq.int(nrow(", variable, "))"
)
evaluateOperation(operation)

command <- paste0(
    variable, " %>% ", "dplyr::relocate(row_n, .before = flag)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Check on mate pairing; generate vectors and tables for mate pairing --------
command <- paste0(
    "testMatesPaired(pos = ", variable, "$qpos, mpos = ", variable, "$qmpos)"
)  # VMP: "vector mates paired"
operation <- makeOperation(paste0("VMP.", variable), command)
evaluateOperation(operation)

#  Check: NAs in VMP vector?
command <- paste0("is.na(", "VMP.", variable, ")")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)

command <- paste0("isNA.VMP.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)  # Yes, 1 NA

command <- paste0("rm(", "isNA.VMP.", variable, ")")
evaluateOperation(command)

#  Convert NAs
#  For "ambiguous" and "mated," assign FALSE to NAs
operation <- paste0(
    "VMP.", variable[!stringr::str_detect(variable, "unmated")], "[",
        "is.na(", "VMP.", variable[!stringr::str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  For "unmated," assign FALSE to NAs
operation <- paste0(
    "VMP.", variable[stringr::str_detect(variable, "unmated")], "[",
        "is.na(", "VMP.", variable[stringr::str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  Check: NAs in VMP vector? Should be gone now
command <- paste0("is.na(", "VMP.", variable, ")")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)

command <- paste0("isNA.VMP.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)  # No NAs

command <- paste0("rm(", "isNA.VMP.", variable, ")")
evaluateOperation(command)

#  Save vectors as tibbles (`TRUE`, `FALSE` columns)
command <- paste0(
    "VMP.", variable, " %>% ",
        "table() %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("tbl.VMP.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles"
vector_string_tibbles <- paste0("tbl.VMP.", variable)
z.check.sub.1.VMP <- createTibbleFromTibbles(vector_string_tibbles)
z.check.sub.1.VMP$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.sub.1.VMP <- z.check.sub.1.VMP %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0("rm(", "tbl.VMP.", variable, ")")
evaluateOperation(command)
rm(vector_string_tibbles)


#  In the mate-pairing test variable, are `TRUE`, `FALSE` *not* alternating? --
#+ (They should be alternating.)

#  That is, `TRUE` means rows *are not* alternating; `FALSE` means they *are*
#+ alternating

#  Generate vectors and tables
command <- paste0(
    "testLogicalAlternating(logical = ", "VMP.", variable, ")"
)  # VA: "vector alternating"
operation <- makeOperation(paste0("VA.", variable), command)
evaluateOperation(operation)

#  Check: NAs in VA vector?
command <- paste0("is.na(", "VA.", variable, ")")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)

command <- paste0("isNA.VA.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)  # Yes, 1 NA

command <- paste0("rm(", "isNA.VA.", variable, ")")
evaluateOperation(command)

#  Convert NAs
#  For "ambiguous" and "mated," assign FALSE to NAs
operation <- paste0(
    "VA.", variable[!stringr::str_detect(variable, "unmated")], "[",
        "is.na(", "VA.", variable[!stringr::str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  For "unmated," assign TRUE to NAs
operation <- paste0(
    "VA.", variable[stringr::str_detect(variable, "unmated")], "[",
        "is.na(", "VA.", variable[stringr::str_detect(variable, "unmated")], ")",
    "]", " <- ", "TRUE"
)
evaluateOperation(operation)

#  Check: NAs in VA vector? Should be gone now
command <- paste0("is.na(", "VA.", variable, ")")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)

command <- paste0("isNA.VA.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)  # No NAs

command <- paste0("rm(", "isNA.VA.", variable, ")")
evaluateOperation(command)

#  Save vectors as tibbles (`TRUE`, `FALSE` columns)
command <- paste0(
    "VA.", variable, " %>% ",
        "table() %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("tbl.VA.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles"
vector_string_tibbles <- paste0("tbl.VA.", variable)
z.check.sub.2.VA <- createTibbleFromTibbles(vector_string_tibbles)
z.check.sub.2.VA$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.sub.2.VA <- z.check.sub.2.VA %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0("rm(", "tbl.VA.", variable, ")")
evaluateOperation(command)
rm(vector_string_tibbles)


#  Create tbl columns denoting lack of alternation, i.e., $discrepancy --------
command <- paste0("VA.", variable)
operation <- makeOperation(paste0(variable, "$discrepancy"), command)
evaluateOperation(operation)

#  In $discrepancy, NA should be FALSE for "mated" and "ambiguous," TRUE for
#+ "unmated;" check via tables
command <- paste0(
    "is.na(", variable, "$discrepancy", ")", " %>% ", "table()" 
)
operation <- makeOperation(paste0("isNA.", variable), command)
evaluateOperation(operation)  # No NAs present

#  Remove the temporary tables
command <- paste0("rm(", "isNA.", variable, ")")
evaluateOperation(command)


#  Make tibble indices (numeric vectors) for discrepancies ± 1 ----------------
command <- paste0(
    "createVectorPlusMinus(",
        "vector = which(", variable, "$discrepancy == TRUE)",
    ")"
)  # VD: "vector discrepancy"
operation <- makeOperation(paste0("VD.", variable), command)
evaluateOperation(operation)


#  Subset tibbles for only discrepancies ± 1 ----------------------------------
command <- paste0(
    variable, "[", "VD.", variable, ", ]"
)
operation <- makeOperation(paste0("discrepancy.", variable), command)
evaluateOperation(operation)

#  Sort the "discrepancy tibbles"
command <- paste0(
    "discrepancy.", variable, "[",
        "order(",
            "discrepancy.", variable, "$isize_abs, ",
            "discrepancy.", variable, "$qpos, ",
            "discrepancy.", variable, "$pos, ",  # Addition 2/21/22
            "discrepancy.", variable, "$qmpos",
        "), ",
    "]"
)
operation <- makeOperation(paste0("discrepancy.", variable), command)
evaluateOperation(operation)


#  From the sub tibbles, filter out duplicate $qpos, $qmpos -------------------

#QUESTION Is this step handled properly?
command <- paste0(
    variable, "[",
        "!duplicated(", variable, "$qpos) & ",
        "!duplicated(", variable, "$qmpos), ",
    "]"
)
operation <- makeOperation(
    paste0("filter.", variable) %>% stringr::str_remove("sub."),
    command
)
evaluateOperation(operation)

#  Remove unneeded variables
command <- paste0(
    "rm(",
        "VA.", variable, ", ",
        "VD.", variable, ", ",
        "VMP.", variable,
    ")"
)
evaluateOperation(command)


#  Are mates paired after filtering out duplicate $qpos, $qmpos? --------------
variable <- paste0("filter.", variable) %>% stringr::str_remove("sub.")

#  Add row numbers to filter.* tibbles
operation <- paste0(
    variable, "$row_n", " <- ", "seq.int(nrow(", variable, "))"
)
evaluateOperation(operation)

#  Generate vectors and tibbles 
command <- paste0(
    "testMatesPaired(",
        "pos = ", variable, "$qpos, ",
        "mpos = ", variable, "$qmpos",
    ")"
)
operation <- makeOperation(paste0("VMP.", variable), command)
evaluateOperation(operation)

#  Check: NAs in VMP vector?
command <- paste0("is.na(", "VMP.", variable, ")")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)

command <- paste0("isNA.VMP.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)  # Yes, 1 NA

command <- paste0("rm(", "isNA.VMP.", variable, ")")
evaluateOperation(command)

#  Convert NAs
#  For "ambiguous" and "mated," assign FALSE to NAss
operation <- paste0(
    "VMP.", variable[!stringr::str_detect(variable, "unmated")], "[",
        "is.na(", "VMP.", variable[!stringr::str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  For "unmated," assign FALSE to NAs
operation <- paste0(
    "VMP.", variable[stringr::str_detect(variable, "unmated")], "[",
        "is.na(", "VMP.", variable[stringr::str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  Check: NAs in VMP vector?
command <- paste0("is.na(", "VMP.", variable, ")")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)

command <- paste0("isNA.VMP.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)  # No NAs

command <- paste0("rm(", "isNA.VMP.", variable, ")")
evaluateOperation(command)

command <- paste0(
    "VMP.", variable, " %>% ",
        "table() %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("tbl.VMP.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles"
vector_string_tibbles <- paste0("tbl.VMP.", variable)
z.check.filter.1.VMP <- createTibbleFromTibbles(vector_string_tibbles)
z.check.filter.1.VMP$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.filter.1.VMP <- z.check.filter.1.VMP %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0("rm(", "tbl.VMP.", variable, ")")
evaluateOperation(command)
rm(vector_string_tibbles)


#  In the mate-pairing test variable, are `TRUE`, `FALSE` *not* alternating? --
#+ (They should be alternating.)

#  That is, `TRUE` means rows *are not* alternating; `FALSE` means they *are*
#+ alternating

#  Create vectors and tibbles
command <- paste0(
    "testLogicalAlternating(logical = ", "VMP.", variable, ")"
)  # VA: "vector alternating"
operation <- makeOperation(paste0("VA.", variable), command)
evaluateOperation(operation)

#  Check: NAs in VA vector?
command <- paste0("is.na(", "VA.", variable, ")")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)

command <- paste0("isNA.VA.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)  # Yes, 1 NA

command <- paste0("rm(", "isNA.VA.", variable, ")")
evaluateOperation(command)

#  Convert NAs
#  For "ambiguous" and "mated," assign FALSE to NAs
operation <- paste0(
    "VA.", variable[!stringr::str_detect(variable, "unmated")], "[",
        "is.na(", "VA.", variable[!stringr::str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  For "unmated," assign TRUE to NAs
operation <- paste0(
    "VA.", variable[stringr::str_detect(variable, "unmated")], "[",
        "is.na(", "VA.", variable[stringr::str_detect(variable, "unmated")], ")",
    "]", " <- ", "TRUE"
)
evaluateOperation(operation)

#  Check: NAs in VA vector?
command <- paste0("is.na(", "VA.", variable, ")")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)

command <- paste0("isNA.VA.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)  # No NAs

command <- paste0("rm(", "isNA.VA.", variable, ")")
evaluateOperation(command)

command <- paste0(
    "VA.", variable, " %>% ",
        "table() %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("tbl.VA.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles"
vector_string_tibbles <- paste0("tbl.VA.", variable)
z.check.filter.2.VA <- createTibbleFromTibbles(vector_string_tibbles)
z.check.filter.2.VA$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.filter.2.VA <- z.check.filter.2.VA %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0("rm(", "tbl.VA.", variable, ")")
evaluateOperation(command)
rm(vector_string_tibbles)


#  Create tbl columns denoting lack of alternation, i.e., "discrepancies" -----
command <- paste0("VA.", variable)
operation <- makeOperation(paste0(variable, "$discrepancy"), command)
evaluateOperation(operation)

#  In $discrepancy, should be no NAs? Check via tables
command <- paste0(
    "is.na(", variable, "$discrepancy", ")", " %>% ", "table()" 
)
operation <- makeOperation(paste0("isNA.", variable), command)
evaluateOperation(operation)  # No NAs present

#  Remove the temporary tables
command <- paste0("rm(", "isNA.", variable, ")")
evaluateOperation(command)


#  Make numeric vectors (i.e., tibble indices) for discrepancies ± 1 ----------
command <- paste0(
    "createVectorPlusMinus(",
        "vector = which(", variable, "$discrepancy == TRUE)",
    ")"
)  # VD: "vector discrepancy"
operation <- makeOperation(paste0("VD.", variable), command)
evaluateOperation(operation)


#  Subset tibbles for only discrepancies ± 1 ----------------------------------
command <- paste0(
    variable, "[", "VD.", variable, ", ]"
)
operation <- makeOperation(paste0("discrepancy.", variable), command)
evaluateOperation(operation)

#  Sort the "discrepancy tibbles"
command <- paste0(
    "discrepancy.", variable, "[",
        "order(",
            "discrepancy.", variable, "$isize_abs, ",
            "discrepancy.", variable, "$qname, ",
            "discrepancy.", variable, "$mpos, ",  # Addition 2/21/22
            "discrepancy.", variable, "$pos",
        "), ",
    "]"
)
operation <- makeOperation(paste0("discrepancy.", variable), command)
evaluateOperation(operation)


#TODO Some automated way to do the following
#  Manually correct filter.a.ambiguous.dedup.CAST -----------------------------
# View(filter.a.ambiguous.dedup.CAST[paste0("584", c(19, 20, 21, 22, 23, 24)), ])
# # A tibble: 6 x 26
#   row_n  flag rname mrnm  strand       pos      mpos   pos_end  mpos_end
#   <int> <int> <fct> <fct> <fct>      <int>     <int>     <dbl>     <dbl>
# 1 58419    99 chrX  chrX  +       40482871  40483010  40482920  40483059
# 2 58420   147 chrX  chrX  -       40483010  40482871  40483059  40482920
# 3 58421   163 chrX  chrX  +        1211507   1211646   1211556   1211695
# 4 58422    99 chrX  chrX  +      121154285 121154424 121154334 121154473
# 5 58423   147 chrX  chrX  -      121154424 121154285 121154473 121154334
# 6 58424    83 chrX  chrX  -        1211646   1211507   1211695   1211556

filter.a.ambiguous.dedup.CAST[paste0("584", c(19, 20, 21, 22, 23, 24)), ] <-
    discrepancy.filter.a.ambiguous.dedup.CAST

#  Check
filter.a.ambiguous.dedup.CAST[c(58419:58429), ]
filter.a.ambiguous.dedup.CAST$row_n[c(58419:58429)]
# [1] 58419 58420 58421 58424 58422 58423 58425 58426 58427 58428 58429

#  Re-label $row_n
filter.a.ambiguous.dedup.CAST$row_n <- seq.int(
    nrow(filter.a.ambiguous.dedup.CAST)
)
filter.a.ambiguous.dedup.CAST[c(58419:58429), ]
filter.a.ambiguous.dedup.CAST$row_n[c(58419:58429)]
# [1] 58419 58420 58421 58422 58423 58424 58425 58426 58427 58428 58429


#  Manually correct filter.a.mated.dedup.129 ----------------------------------
# discrepancy.filter.a.mated.dedup.129
# # A tibble: 6 x 26
#   row_n  flag rname mrnm  strand       pos      mpos   pos_end  mpos_end
#   <int> <int> <fct> <fct> <fct>      <int>     <int>     <dbl>     <dbl>
# 1 34379   163 chrX  chrX  +      100268668 100268780 100268717 100268829
# 2 34380    83 chrX  chrX  -      100268780 100268668 100268829 100268717
# 3 34381    99 chrX  chrX  +      130098649 130098761 130098698 130098810
# 4 34382    99 chrX  chrX  +      130098654 130098766 130098703 130098815
# 5 34383   147 chrX  chrX  -      130098761 130098649 130098810 130098698
# 6 34384   147 chrX  chrX  -      130098766 130098654 130098815 130098703

desired_order <- paste0("343", c(79, 80, 81, 83, 82, 84))
discrepancy.filter.a.mated.dedup.129.adjust <-
    discrepancy.filter.a.mated.dedup.129 %>%
        dplyr::mutate(row_n = factor(row_n, levels = desired_order)) %>%
        dplyr::arrange(row_n)
discrepancy.filter.a.mated.dedup.129.adjust$row_n

discrepancy.filter.a.mated.dedup.129.adjust$row_n <-
    discrepancy.filter.a.mated.dedup.129.adjust$row_n %>%
        as.character() %>%
        as.integer()

filter.a.mated.dedup.129[discrepancy.filter.a.mated.dedup.129$row_n, ] <-
    discrepancy.filter.a.mated.dedup.129.adjust

discrepancy.filter.a.mated.dedup.129 <-
    discrepancy.filter.a.mated.dedup.129.adjust

#  Check
filter.a.mated.dedup.129[c(34379:34389), ]
filter.a.mated.dedup.129$row_n[c(34379:34389)]
# [1] 34379 34380 34381 34383 34382 34384 34385 34386 34387 34388 34389

#  Re-label $row_n
filter.a.mated.dedup.129$row_n <- seq.int(nrow(filter.a.mated.dedup.129))

#  Check
filter.a.mated.dedup.129[c(34379:34389), ]
filter.a.mated.dedup.129$row_n[c(34379:34389)]
# [1] 34379 34380 34381 34382 34383 34384 34385 34386 34387 34388 34389

#  Remove unneeded variables
rm(desired_order, discrepancy.filter.a.mated.dedup.129.adjust)


#  Clean up unneeded variables before setting up new ones ---------------------
command <- paste0(
    "rm(",
        "VMP.", variable, ", ",
        "VA.", variable, ", ",
        "VD.", variable,
    ")"
)
evaluateOperation(command)


#  Set up new variables -------------------------------------------------------
command <- paste0(variable)
operation <- makeOperation(
    paste0("adjust.", stringr::str_remove(variable, "filter.")),
    command
)
evaluateOperation(operation)

variable <- paste0("adjust.", stringr::str_remove(variable, "filter."))


#  Generate vectors and tables for mate-pairing readouts ----------------------
command <- paste0(
    "testMatesPaired(pos = ", variable, "$qpos, mpos = ", variable, "$qmpos)"
)  # VMP: "vector mates paired"
operation <- makeOperation(paste0("VMP.", variable), command)
evaluateOperation(operation)

#  Check: NAs in VMP vector?
command <- paste0("is.na(", "VMP.", variable, ")")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)

command <- paste0("isNA.VMP.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)  # Yes, 1 NA

command <- paste0("rm(", "isNA.VMP.", variable, ")")
evaluateOperation(command)

#  Convert NAs
#  For "ambiguous" and "mated," assign FALSE to NAs
operation <- paste0(
    "VMP.", variable[!stringr::str_detect(variable, "unmated")], "[",
        "is.na(", "VMP.", variable[!stringr::str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  For "unmated," assign FALSE to NAs
operation <- paste0(
    "VMP.", variable[stringr::str_detect(variable, "unmated")], "[",
        "is.na(", "VMP.", variable[stringr::str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  Check: NAs in VMP vector? Should be gone now
command <- paste0("is.na(", "VMP.", variable, ")")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)

command <- paste0("isNA.VMP.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)  # No NAs

command <- paste0("rm(", "isNA.VMP.", variable, ")")
evaluateOperation(command)

#  Save vectors as tibbles (`TRUE`, `FALSE` columns)
command <- paste0(
    "VMP.", variable, " %>% ",
        "table() %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("tbl.VMP.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles"
vector_string_tibbles <- paste0("tbl.VMP.", variable)
z.check.adjust.1.VMP <- createTibbleFromTibbles(vector_string_tibbles)
z.check.adjust.1.VMP$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.adjust.1.VMP <- z.check.adjust.1.VMP %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0("rm(", "tbl.VMP.", variable, ")")
evaluateOperation(command)
rm(vector_string_tibbles)


#  In the mate-pairing test variable, are `TRUE`, `FALSE` *not* alternating? --
#+ (They should be alternating.)

#  Generate vectors and tables
command <- paste0(
    "testLogicalAlternating(logical = ", "VMP.", variable, ")"
)  # VA: "vector alternating"
operation <- makeOperation(paste0("VA.", variable), command)
evaluateOperation(operation)

#  Check: NAs in VA vector?
command <- paste0("is.na(", "VA.", variable, ")")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)

command <- paste0("isNA.VA.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)  # Yes, 1 NA

command <- paste0("rm(", "isNA.VA.", variable, ")")
evaluateOperation(command)

#  Convert NAs
#  For "ambiguous" and "mated," assign FALSE to NAs
operation <- paste0(
    "VA.", variable[!stringr::str_detect(variable, "unmated")], "[",
        "is.na(", "VA.", variable[!stringr::str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  For "unmated," assign TRUE to NAs
operation <- paste0(
    "VA.", variable[stringr::str_detect(variable, "unmated")], "[",
        "is.na(", "VA.", variable[stringr::str_detect(variable, "unmated")], ")",
    "]", " <- ", "TRUE"
)
evaluateOperation(operation)

#  Check: NAs in VA vector? Should be gone now
command <- paste0("is.na(", "VA.", variable, ")")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)

command <- paste0("isNA.VA.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)  # No NAs

command <- paste0("rm(", "isNA.VA.", variable, ")")
evaluateOperation(command)

#  Save vectors as tibbles (`TRUE`, `FALSE` columns)
command <- paste0(
    "VA.", variable, " %>% ",
        "table() %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("tbl.VA.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles"
vector_string_tibbles <- paste0("tbl.VA.", variable)
z.check.adjust.2.VA <- createTibbleFromTibbles(vector_string_tibbles)
z.check.adjust.2.VA$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.adjust.2.VA <- z.check.adjust.2.VA %>% dplyr::arrange(rownames)


#  Remove unneeded variables --------------------------------------------------
variable_non_adjust <- variable %>% stringr::str_remove("adjust.")
variable_non_a <- variable_non_adjust %>% stringr::str_remove("a.")
operation <- paste0(
    "rm(",
        "VMP.", variable, ", ",
        "VA.", variable, ", ",
        "tbl.VA.", variable, ", ",
        "same.", variable_non_adjust, ", ",
        "sub.", variable_non_adjust, ", ",
        "discrepancy.sub.", variable_non_adjust, ", ",
        "discrepancy.filter.", variable_non_adjust, ", ",
        "filter.", variable_non_adjust, ", ",
        variable_non_adjust, ", ",
        vector_string_tibbles, ", ",
        variable_non_a, ", ",
        variable_non_adjust,
    ")"
)
evaluateOperation(operation)

rm(variable_non_a, variable_non_adjust, vector_string_tibbles)

command <- paste0(
    "rm(",
        paste0("z.check.", c(
            "1.MIP", "1.PIM",
            "2.PeqM", "2.PgtM", "2.PltM",
            "sub.0.MIP", "sub.0.PIM", "sub.1.VMP", "sub.2.VA",
            "filter.1.VMP", "filter.2.VA",
            "adjust.1.VMP", "adjust.2.VA"
        )),
    ")"
)
evaluateOperation(command)


#  Generate $groupid_2 for munged "mated" and "ambiguous" tibbles -------------

#  Work with only "mated" and "ambiguous" tibbles, not "unmated" tibbles
variable <- variable[!stringr::str_detect(variable, "unmated")]

for (i in 1:length(variable)) {
    df <- eval(parse(text = variable[i]))
    odd <- seq(1, nrow(df), 2)
    even <- (seq(1, nrow(df), 2)) + 1
    
    odd <- df[odd, ]
    even <- df[even, ]
    
    odd$tibble <- "1"
    even$tibble <- "2"
    
    odd <- odd %>% dplyr::relocate(tibble, .before = row_n)
    even <- even %>% dplyr::relocate(tibble, .before = row_n)
    
    command <- paste0(
        "odd %>%",
            "dplyr::mutate(groupid_2 = row_number()) %>%",
            "dplyr::bind_rows(even %>% mutate(groupid_2 = row_number())) %>%",
            "dplyr::arrange(groupid_2, tibble) %>%",
            "dplyr::relocate(groupid_2, .after = groupid) %>%",
            "dplyr::select(-tibble)"
    )
    operation <- makeOperation(variable[i], command)
    evaluateOperation(operation)
    
    n_occur <- data.frame(table(eval(parse(text = variable[i]))$groupid_2))
    print(n_occur[n_occur$Freq > 2, ])
}

rm(df, odd, even, n_occur)


#  Combine "mated" and "ambiguous" tibbles ------------------------------------
tbl.129 <- stringr::str_subset(variable, "129") %>% bindRows()
tbl.CAST <- stringr::str_subset(variable, "CAST") %>% bindRows()
tbl.mm10 <- stringr::str_subset(variable, "mm10") %>% bindRows()

#  Remove unneeded data from environment
rm(bindRows)


#  Generate groupid.* variables for combination mated/ambiguous tibbles -------
variable_tbl <- paste0("tbl.", c("129", "CAST", "mm10"))

for (i in 1:length(variable_tbl)) {
    #  Assign tibble of interest to variable 'df'
    df <- eval(parse(text = variable_tbl[i]))
    
    #  Create tibbles from odd and even rows
    odd.seq <- seq(1, nrow(df), 2)
    even.seq <- (seq(1, nrow(df), 2)) + 1
    
    odd <- df[odd.seq, ]
    even <- df[even.seq, ]
    
    #  Odd tibble is tibble #1, even tibble is tibble #2
    odd$tibble <- "1"
    even$tibble <- "2"
    
    odd <- odd %>% dplyr::relocate(tibble, .before = row_n)
    even <- even %>% dplyr::relocate(tibble, .before = row_n)
    
    #  Interleave the odd/even rows, then arrange them by group ID and tibble
    #+ number
    command <- paste0(
        "odd %>%",
            "dplyr::mutate(groupid_3 = row_number()) %>% ",
            "dplyr::bind_rows(even %>% mutate(groupid_3 = row_number())) %>% ",
            "dplyr::arrange(groupid_3, tibble) %>% ",
            "dplyr::relocate(groupid_3, .after = groupid_2) %>% ",
            "dplyr::select(-tibble)"
    )
    operation <- makeOperation(variable_tbl[i], command)
    evaluateOperation(operation)
    
    #  Check to make sure that there are no more than two entries per group ID
    n_occur <- data.frame(table(eval(parse(text = variable_tbl[i]))$groupid_3))
    print(n_occur[n_occur$Freq > 2, ])
    rm(n_occur)
    
    #  Sort the tibble of interest by pos while maintaining proper mate pairs
    df <- eval(parse(text = variable_tbl[i]))
    
    odd <- df[odd.seq, ]
    even <- df[even.seq, ]
    
    #  Create a tibble with two mates per row (joined by "groupid_3"); then,
    #+ sort by pos.x (the odd-row pos value)
    df <- full_join(odd, even, by = "groupid_3") %>%
        dplyr::arrange(pos.x) %>%  # Sort by pos.x
        dplyr::rename(groupid_3.x = groupid_3)  # new_name = old_name
    df$groupid_3.y <- df$groupid_3.x
    df <- df %>% dplyr::relocate(groupid_3.y, .after = groupid_2.y)
    
    #  Rename column names to denote odd/even status
    colnames(df) <- str_replace_all(colnames(df), "\\.x", "\\.odd")
    colnames(df) <- str_replace_all(colnames(df), "\\.y", "\\.even")
    
    #  Split the tibble by odd or even status
    odd <- df[stringr::str_subset(colnames(df), "\\.odd")]
    even <- df[stringr::str_subset(colnames(df), "\\.even")]
    
    #  Strip suffixes from column names
    colnames(odd) <- str_replace_all(colnames(odd), "\\.odd", "")
    colnames(even) <- str_replace_all(colnames(even), "\\.even", "")
    
    #  Odd tibble is tibble #1, even tibble is tibble #2
    odd$tibble <- "1"
    even$tibble <- "2"
    
    #  Again, interleave the odd/even rows, then arrange them by group ID and
    #+ tibble number
    command <- paste0(
        "odd %>%",
            "dplyr::mutate(groupid_4 = row_number()) %>% ",
            "dplyr::bind_rows(even %>% mutate(groupid_4 = row_number())) %>% ",
            "dplyr::arrange(groupid_4, tibble) %>% ",
            "dplyr::relocate(groupid_4, .after = groupid_3) %>% ",
            "dplyr::select(-tibble)"
    )
    operation <- makeOperation(variable_tbl[i], command)
    evaluateOperation(operation)
    
    #  Again, check to make sure that there are no more than two entries per
    #+ group ID
    n_occur <- data.frame(table(eval(parse(text = variable_tbl[i]))$groupid_4))
    print(n_occur[n_occur$Freq > 2, ])
    rm(n_occur)
}

rm(df, i, even, even.seq, odd, odd.seq)

#  Remove unneeded columns from the tbl.* tibbles
command <- paste0(
    variable_tbl, " %>% ",
        "dplyr::select(-groupid, -groupid_2, -groupid_3, -row_n)", " %>% ",
        "dplyr::rename(groupid = groupid_4)"
)
operation <- makeOperation(variable_tbl, command)
evaluateOperation(operation)


#  Clean up global environment, save environment image ------------------------
operation <- paste0("rm(", variable, ")")
evaluateOperation(operation)

variable <- variable_tbl
rm(chromosomes, variable_tbl)
rm(
    adjust.a.unmated.dedup.129,
    adjust.a.unmated.dedup.CAST,
    adjust.a.unmated.dedup.mm10
)

save.image(paste0(
        arguments$directory_out, "/",
        stringr::str_remove(script, ".R"), "_part-2", ".Rdata"
))


#  Script is completed; clean up environment ----------------------------------
rm(list = ls())
