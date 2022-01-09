#!/usr/bin/env Rscript

library(parallel)
library(Rsamtools)
library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)


# -----------------------------------------------------------------------------
#  Temporary: Set up work directory (location TB∆)
"/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora" %>% setwd()


#  Set up functions -----------------------------------------------------------
`%notin%` <- Negate(`%in%`)


createVectorPlusMinus <- function(vector) {
    #  Takes a numeric vector and returns a vector of numbers ± one
    v1 <- vector
    v2 <- v1 - 1
    v3 <- v1 + 1
    vector_plus_minus <- c(v1, v2, v3) %>% sort() %>% unique()
    return(vector_plus_minus)
}


makeOperation <- function(variable, command) {
    operation <- paste(variable, command) %>% paste(., collapse = "; ")
    return(operation)
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


#  Assign variables necessary for subsequent commands -------------------------
chromosome <- "chrX"

file <- list.files(pattern = paste0("\\", chromosome, ".rmdup.bam$"))
file <- file[c(1, 2, 3)]
variable <- file %>%
    strsplit(., "-") %>%
    lapply(., `[[`, 1) %>% unlist() %>%
    paste0("dedup.", .)
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

file <- list.files(pattern = paste0("\\", chromosome, ".rmdup.bam.bai$"))
file <- file[c(1, 2, 3)]
index <- file %>%
    strsplit(., "-") %>%
    lapply(., `[[`, 1) %>% unlist() %>%
    paste0("index.", .)
mapply(
    assign, index, file, MoreArgs = list(envir = parent.frame())
)


#  Assign .bam information ----------------------------------------------------
#  Load in standard .bam fields
command <- paste0(
    "<- ", variable,
    " %>% Rsamtools::BamFile(., index = ", index, ", asMates = TRUE)",
    " %>% Rsamtools::scanBam()",
    " %>% as.data.frame()",
    " %>% tibble::as_tibble()"
)
operation <- makeOperation(paste0(variable, ".full"), command)
eval(parse(text = operation))

#  Load in .bam AS field
map_params <- Rsamtools::ScanBamParam(tag = "AS")
command <- paste0(
    "<- ", variable,
    " %>% Rsamtools::BamFile(., index = ", index, ", asMates = TRUE)",
    " %>% Rsamtools::scanBam(param = map_params)",
    " %>% as.data.frame()",
    " %>% tibble::as_tibble()"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))
rm(map_params)

#  Join the standard and AS fields
command <- paste0(
    "<- dplyr::bind_cols(", variable, ".full, ", variable, ")"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))

#  Remove variables no longer needed
command <- paste0("rm(", variable, ".full)")
eval(parse(text = command))


#  Check .bam flag and mate_status attributes ---------------------------------
#  Check $flag
variable_flag <- paste0("flag.", variable)
command <- paste0("<- ", variable, "$flag %>% table()")
operation <- makeOperation(variable_flag, command)
eval(parse(text = operation))

#  Remove variables no longer needed
command <- paste0("rm(", variable_flag, ")")
eval(parse(text = command))
rm(variable_flag)

#  Check $mate_status
variable_status <- paste0("status.", variable)
command <- paste0("<- ", variable, "$mate_status %>% table()")
operation <- makeOperation(variable_status, command)
eval(parse(text = operation))

#  Print out mate statuses
for (i in 1:length(variable_status)) {
    cat(variable_status[i], "\n")
    cat("mated ambiguous unmated\n")
    command <- paste0(variable_status[i])
    eval(parse(text = command)) %>% cat()
    cat("\n\n")
}

#  Remove variables no longer needed
# command <- paste0("rm(", variable_status, ")")
# eval(parse(text = command))
rm(i, variable_status)


#  Split tibbles based on $mate_status ----------------------------------------
split <- c("ambiguous", "mated", "unmated")
for (i in 1:length(split)) {
    variable_split <- paste0(split[i], ".", variable)
    command <- paste0(
        "<- ", variable, "[", variable, "$mate_status == '", split[i], "', ]"
    )
    operation <- makeOperation(variable_split, command)
    eval(parse(text = operation))
}

#  Create a vector of names of "split variables"
split <- paste0(split, ".")
variable <- purrr::cross(list(split, variable)) %>% purrr::map_chr(paste0, collapse = "")

#  Clean up variables
rm(i, split, variable_split)

#  For each "split variable", print out duplicated status based on $pos
# for (i in 1:length(variable)) {
#     paste0(variable[i], "$pos") %>% cat("\n")
#     cat("FALSE TRUE\n")
#     command <- paste0(variable[i], "$pos %>% duplicated() %>% table()")
#     eval(parse(text = command)) %>% cat()
#     cat("\n\n")
# }

variable <- variable.bak
variable.bak <- variable

#  Reorder rname factor levels ------------------------------------------------
command <- paste0("<- forcats::fct_relevel(", variable, "$rname, chromosome)")
operation <- makeOperation(paste0(variable, "$rname"), command)
eval(parse(text = operation))


#  Drop unused rname factor levels --------------------------------------------
command <- paste0("<- ", variable, "$rname %>% forcats::fct_drop()")
operation <- makeOperation(paste0(variable, "$rname"), command)
eval(parse(text = operation))


#  Drop rows that are not chr1-19, chrX ---------------------------------------
chromosomes <- c(paste0("chr", 1:19), "chrX")
command_pipe <- paste0("<- ", variable, " %>% ")
command <- paste0(command_pipe, "filter(., rname %in% chromosomes)")
operation <- makeOperation(variable, command)
eval(parse(text = operation))


#  Create and append pos_end, mpos_end columns --------------------------------
command <- paste0("<- ", variable, "$pos + 49")
operation <- makeOperation(paste0(variable, "$pos_end"), command)
eval(parse(text = operation))

command <- paste0("<- ", variable, "$mpos + 49")
operation <- makeOperation(paste0(variable, "$mpos_end"), command)
eval(parse(text = operation))

command <- paste0(
    command_pipe,
    "dplyr::relocate(pos_end, .after = pos) %>% ",
    "dplyr::relocate(mpos_end, .after = mpos)"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))


#  Filter out rows with mapq less than 30 -------------------------------------
command <- paste0("<- ", variable, "[", variable, "$mapq >= 30, ]")
operation <- makeOperation(variable, command)
eval(parse(text = operation))


#  Prior to sorting, deduplicate and mate-label the tibbles -------------------
#  Set up $criteria
command <- paste0(
    "<- paste0(",
    variable, "$qname, ", "'_', ",
    variable, "$flag, ", "'_', ",
    variable, "$pos, ", "'_', ",
    variable, "$mpos",
    ")"
)
operation <- makeOperation(paste0(variable, "$criteria"), command)
eval(parse(text = operation))

#  Set up $qpos
command <- paste0(
    "<- paste0(",
    variable, "$qname, ", "'_', ",
    variable, "$pos",
    ")"
)
operation <- makeOperation(paste0(variable, "$qpos"), command)
eval(parse(text = operation))

#  Set up $qmpos
command <- paste0(
    "<- paste0(",
    variable, "$qname, ", "'_', ",
    variable, "$mpos",
    ")"
)
operation <- makeOperation(paste0(variable, "$qmpos"), command)
eval(parse(text = operation))

#  Set up $isize_abs, i.e., absolute insert-size values
command <- paste0(
    "<- ", variable, "$isize %>% abs() %>% as.numeric()"
)
operation <- makeOperation(paste0(variable, "$isize_abs"), command)
eval(parse(text = operation))

# Create back-ups of variables
command <- paste0("<- ", variable)
operation <- makeOperation(paste0("bak.", variable), command)
eval(parse(text = operation))

#  Load back-ups of variables
command <- paste0("<- bak.", variable)
operation <- makeOperation(variable, command)
eval(parse(text = operation))

#  Rearrange the order of tibble columns for quick, easy surveying
command <- paste0(
    "<- ", variable, " %>% ",
    "dplyr::relocate(mpos, .after = pos) %>% ",
    "dplyr::relocate(mpos_end, .after = pos_end) %>% ",
    "dplyr::relocate(mrnm, .after = rname) %>% ",
    "dplyr::relocate(isize, .after = mpos_end) %>% ",
    "dplyr::relocate(criteria, .before = groupid) %>% ",
    "dplyr::relocate(AS, .after = cigar) %>% ",
    "dplyr::relocate(qpos, .after = criteria) %>% ",
    "dplyr::relocate(qmpos, .after = qpos) %>% ",
    "dplyr::relocate(qname, .after = qual) %>% ",
    "dplyr::relocate(isize_abs, .before = isize)"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))

#  Convert $groupid to factor  #MAYBE
command <- paste0("<- ", variable, "$groupid %>% as_factor()")
operation <- makeOperation(paste0(variable, "$groupid"), command)
eval(parse(text = operation))

#  Calculate the numbers of entries per $groupid level
# command <- paste0(
#     "<- ", variable, " %>% ",
#     "dplyr::group_by(groupid) %>% ",
#     "summarize(no_rows = length(groupid))"
# )
# operation <- makeOperation(paste0("groupid.", variable), command)
# eval(parse(text = operation))

#  Clean up groupid.* variables
command <- paste0("rm(groupid.", variable, ")")
eval(parse(text = command))

#  Split out entries in which $qpos is the same as $qmpos, i.e., entries in
#+ which each member of the pair maps to the same location
command <- paste0(
    "<- ", variable, "[", variable, "$qpos == ", variable, "$qmpos, ]"
)
operation <- makeOperation(paste0("same.", variable), command)
eval(parse(text = operation))

#  Regarding the main tibbles, remove entries in which $qpos is the same as
#+ $qmpos, i.e., entries in which each member of the pair maps to the same
#+ location
command <- paste0(
    "<- ", variable, "[!(", variable, "$qpos == ", variable, "$qmpos), ]"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))

#HERE - 1: OK
#  Deduplicate the main tibbles based on $criteria: There should be no more
#+ than one entry with a given $criteria (combination of $qname, $flag, $pos,
#+ and $mpos)
command <- paste0(
    "<- ", variable, "[!duplicated(", variable, "$criteria), ]"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))
#HERE - 2: OK

# Check
# nrow(bak.ambiguous.dedup.129S1) - nrow(ambiguous.dedup.129S1)


#  For each tibble, score the numbers of entries after removing "same" entires
#+ and deduplication
command <- paste0(
    "<- nrow(bak.", variable, ") - nrow(", variable, ")"
)
operation <- makeOperation(paste0("change.", variable), command)
eval(parse(text = operation))

#  Does the length() of change.mated.* equal the nrow() of same.mated.*?
command <- paste0(
    "<- all(nrow(same.", variable, ") == ", "change.", variable, ")"
)
operation <- makeOperation(paste0("equal.", variable), command)
eval(parse(text = operation))

#  Clean up unneeded variables
command <- paste0(
    "rm(", "change.", variable, ", equal.", variable, ")"
)
eval(parse(text = command))


# -------------------------------------
#  Is $qpos in $qmpos? Vice versa?
command <- paste0(
    "<- testPosInMpos(pos = ", variable, "$qpos, mpos = ", variable, "$qmpos) %>% ",
    "table()"
)
operation <- makeOperation(paste0("PIM.", variable), command)
eval(parse(text = operation))

command <- paste0(
    "<- testMposInPos(pos = ", variable, "$qpos, mpos = ", variable, "$qmpos) %>% ",
    "table()"
)
operation <- makeOperation(paste0("MIP.", variable), command)
eval(parse(text = operation))

command <- paste0(
    "<- PIM.", variable, " %>% ",
    "t() %>% ",
    "cbind() %>% ",
    "as_tibble()"
)
operation <- makeOperation(paste0("PIM.", variable), command)
eval(parse(text = operation))

command <- paste0(
    "<- MIP.", variable, " %>% ",
    "t() %>% ",
    "cbind() %>% ",
    "as_tibble()"
)
operation <- makeOperation(paste0("MIP.", variable), command)
eval(parse(text = operation))

command <- paste0(
    "<- PIM.", variable, "$`FALSE` == MIP.", variable, "$`FALSE`"
)
operation <- makeOperation(paste0("PeqM.", variable), command)
eval(parse(text = operation))

command <- paste0(
    "<- PIM.", variable, "$`FALSE` > MIP.", variable, "$`FALSE`"
)
operation <- makeOperation(paste0("PgtM.", variable), command)
eval(parse(text = operation))

command <- paste0(
    "<- PIM.", variable, "$`FALSE` < MIP.", variable, "$`FALSE`"
)
operation <- makeOperation(paste0("PltM.", variable), command)
eval(parse(text = operation))
#HERE - 3

# -------------------------------------
#  Make logical vectors for $qpos in $qmpos and vice versa
command <- paste0(
    "<- testPosInMpos(pos = ", variable, "$qpos, mpos = ", variable, "$qmpos)"
)
operation <- makeOperation(paste0("vPIM.", variable), command)
eval(parse(text = operation))

command <- paste0(
    "<- testMposInPos(pos = ", variable, "$qpos, mpos = ", variable, "$qmpos)"
)
operation <- makeOperation(paste0("vMIP.", variable), command)
eval(parse(text = operation))
#HERE - 4

#  Based on PeqM.*, PgtM.*, and PltM.* results, use these vectors
#+ 
#+ All vPIM.* are fine, apparently


#  Remove unneeded variables
command <- paste0(
    "rm(", "PIM.", variable,
    ", MIP.", variable,
    ", PeqM.", variable,
    ", PgtM.", variable,
    ", PltM.", variable,
    ", vMIP.", variable,
    ")"
)
eval(parse(text = command))


#  Subset tibbles based on logical vectors for $qpos in $qmpos ----------------
command <- paste0(
    "<- ", variable, "[", "vPIM.", variable, ", ]"
)
operation <- makeOperation(paste0("test.", variable), command)
eval(parse(text = operation))
#HERE - 5: OK

pos <- paste0(test.ambiguous.dedup.129S1$qname, "_", test.ambiguous.dedup.129S1$pos)
mpos <- paste0(test.ambiguous.dedup.129S1$qname, "_", test.ambiguous.dedup.129S1$mpos)

pos.in.mpos <- pos %in% mpos
mpos.in.pos <- mpos %in% pos

pos.in.mpos %>% table()
# .
#   TRUE
# 108138
mpos.in.pos %>% table()
# .
#   TRUE
# 108138
#HERE - 6: OK
#HERE - 6.5: OK


#  Order the tibbles based on $isize_abs, $qpos, and $qmpos -------------------
variable <- paste0("test.", variable)
command <- paste0(
    "<- ", variable, "[order(",
    variable, "$isize_abs, ",
    variable, "$qpos, ",
    variable, "$qmpos), ",
    "]"
)
# operation <- makeOperation(paste0("test_2.", variable), command)
operation <- makeOperation(variable, command)
eval(parse(text = operation))
#HERE - 7: OK

vector_mates_paired <- testMatesPaired(
    pos = test.ambiguous.dedup.129S1$qpos,
    mpos = test.ambiguous.dedup.129S1$qmpos
)
vector_mates_paired %>% table()
# .
# FALSE  TRUE 
# 54070 54067
#HERE - 8: OK

vector_alternating <- testLogicalAlternating(logical = vector_mates_paired)
vector_alternating %>% table()
# .
# FALSE   TRUE 
# 108132      4
#HERE - 9: OK

test.ambiguous.dedup.129S1$discrepancy <- vector_alternating
vector_discrepancy <- createVectorPlusMinus(
    vector = which(test.ambiguous.dedup.129S1$discrepancy == TRUE)
)
#HERE - 10: OK

discrepancy.ambiguous.dedup.129S1 <- test.ambiguous.dedup.129S1[vector_discrepancy, ]
discrepancy.ambiguous.dedup.129S1
#HERE - 11: OK

test.ambiguous.dedup.129S1 <- test.ambiguous.dedup.129S1[
    !duplicated(test.ambiguous.dedup.129S1$pos) &
    !duplicated(test.ambiguous.dedup.129S1$mpos), 
]
#HERE - 12: OK

vector_mates_paired <- testMatesPaired(
    pos = test.ambiguous.dedup.129S1$qpos,
    mpos = test.ambiguous.dedup.129S1$qmpos
)
vector_mates_paired %>% table()
# .
# FALSE  TRUE 
# 53520 53521
#HERE - 13: OK

vector_alternating <- testLogicalAlternating(logical = vector_mates_paired)
vector_alternating %>% table()
# .
# FALSE
# 107040
#HERE - 14: OK
