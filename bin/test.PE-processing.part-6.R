#!/usr/bin/env Rscript

library(BSgenome)
library(ggplot2)
library(Mmusculus129S1SvImJ)
library(MmusculusCASTEiJ)
library(Mmusculus129Inserted)
library(MmusculusCASTInserted)
library(MmusculusCAST129Inserted)
library(Mmusculus129Nmasked)
library(MmusculusCASTNmasked)
library(MmusculusCAST129Nmasked)
library(magrittr)
library(pheatmap)
library(Rsamtools)
library(scales)
library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)


#  Set up functions -----------------------------------------------------------
`%notin%` <- Negate(`%in%`)


convertPercent <- function(x) {
    if(is.numeric(x)) {
        ifelse(is.na(x), x, paste0(round(x * 100L, 2), "%"))
    } else x
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


loadRDS <- readRDS


loadTibbleFromRDS <- function(variable, file) {
    command <- paste0(variable, " <- loadRDS(file = \"", file, "\")") %>% as.list()
    return(eval(parse(text = command), envir = globalenv()))
}


makeOperation <- function(variable, command) {
    operation <- paste(variable, command)
    return(operation)
}


mungeMatesIntoOneRow <- function(variable_in, sort_by) {
    #  Input should be of class "character"
    if (class(variable_in) != "character") {
        stop("\"variable_in\" should be class \"character\"")
    } else {
        variable_out <- eval(parse(text = variable_in))
    }
    
    if (class(sort_by) != "character") {
        stop("\"sort_by\" should be class \"character\"")
    }

    odd.seq <- seq(1, nrow(variable_out), 2)
    even.seq <- (seq(1, nrow(variable_out), 2)) + 1
    
    odd <- variable_out[odd.seq, ]
    even <- variable_out[even.seq, ]
    
    odd <- odd %>%
        dplyr::mutate(groupid = row_number()) %>%
        dplyr::relocate(groupid, .after = assignment)
    even <- even %>%
        dplyr::mutate(groupid = row_number()) %>%
        dplyr::relocate(groupid, .after = assignment)
    
    variable_out <- full_join(odd, even, by = "groupid") %>%
        dplyr::arrange(sort_by) %>%  # Sort by pos.x
        dplyr::rename(groupid.x = groupid)  # new_name = old_name
    variable_out$groupid.y <- variable_out$groupid.x
    variable_out <- variable_out %>% dplyr::relocate(groupid.y, .after = assignment.y)
    
    colnames(variable_out) <- str_replace_all(
        colnames(variable_out), "\\.x", "\\.odd"
    )
    colnames(variable_out) <- str_replace_all(
        colnames(variable_out), "\\.y", "\\.even"
    )
    
    return(variable_out)
}


mungeMatesIntoTwoRows <- function(
    variable_in, sort_by_rname, sort_by_pos, sort_by_mpos, position_after
) {
    #  All function arguments should be of class "character"
    variable_out <- eval(parse(text = variable_in))
    
    #  Sort the uniline tibbles by rname, pos, mpos
    variable_out <- variable_out %>%
        dplyr::arrange(sort_by_rname, sort_by_pos, sort_by_mpos)
    
    #  Split the uniline tibbles by odd or even status
    odd <- variable_out[
        stringr::str_subset(colnames(variable_out), "\\.odd")
    ]
    even <- variable_out[
        stringr::str_subset(colnames(variable_out), "\\.even")
    ]
    
    #  Strip suffixes from column names
    colnames(odd) <- str_replace_all(colnames(odd), "\\.odd", "")
    colnames(even) <- str_replace_all(colnames(even), "\\.even", "")
    
    #  Odd tibble is tibble #1, even tibble is tibble #2
    odd$tibble <- "1"
    even$tibble <- "2"
    
    #  Again, interleave the odd/even rows, then arrange them by group ID and
    #+ tibble number; save results to initial "tbl.*" variables
    variable_out <- odd %>% 
        dplyr::mutate(groupid = row_number()) %>%
        dplyr::bind_rows(even %>% mutate(groupid = row_number())) %>%
        dplyr::arrange(groupid, tibble) %>%
        dplyr::relocate(
            groupid, .after = tidyselect::all_of(position_after)
        ) %>%
        dplyr::select(-tibble)
    
    return(variable_out)
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
directory_bam <- "results/kga0/2022-0207_segregated_reads.thresh_SNP_1.thresh_Q_30"
directory_data <- "data/kga0"
path.1 <- paste0(directory_user, "/", directory_project)
path.2 <- paste0(directory_bam)
path.3 <- paste0(directory_data)

rm(directory_bam, directory_data, directory_project, directory_user)

setwd(paste0(path.1))


#  Script name ----------------------------------------------------------------
script <- "test.PE-processing.part-6.R"

#  Make sure to have run 'test.PE-processing.part-4.R',
#+ 'run_GB-assignment-pipeline.2022-0207.sh',
#+ 'generate-sam2pairwise-files_4dn-mouse-cross_GB.sh', and
#+ 'copy-important-files.sh' before running the rest of this script

#  Set up what to query from .bam files ---------------------------------------
map_info <- c(
    "qname", "flag", "rname", "pos", "mapq", "cigar",
    "mrnm", "mpos", "isize", "seq", "qual"
)
map_params <- Rsamtools::ScanBamParam(what = map_info)


#  Munge the KA-generated GB .bam files ---------------------------------------
setwd(paste0(path.1, "/", path.2))

# chromosome <- "chr1"
chromosome <- "chrX"

file <- list.files(pattern = "\\.bam")
variable_GB <- c(
    "GB.alt.CAST",
    "GB.ambig",
    "GB.contra",
    "GB.ref.129S1"
)

#  Check that the correct files will be assigned to the correct variables
mapply(
    assign, variable_GB, file, MoreArgs = list(envir = parent.frame())
)


#  Assign .bam information as list --------------------------------------------
command <- paste0(
    "<- ", variable_GB, " %>% ",
        "Rsamtools::scanBam(., param = map_params)"
)
operation <- makeOperation(variable_GB, command)
evaluateOperation(operation)

rm(map_info, map_params)


#  Convert .bam information from list to dataframe to tibble ------------------
command <- paste0(
    "<- ", variable_GB, " %>% ",
        "as.data.frame() %>% ",
        "as_tibble()"
)
operation <- makeOperation(variable_GB, command)
evaluateOperation(operation)


#  Reorder rname factor levels ------------------------------------------------
chromosomes <- c(paste0("chr", c(1:19)), "chrX", "chrY", "chrM")
command <- paste0(
    "<- forcats::fct_relevel(", variable_GB, "$rname, chromosomes)"
)
operation <- makeOperation(paste0(variable_GB, "$rname"), command)
evaluateOperation(operation)


#  Drop unused rname factor levels --------------------------------------------
command <- paste0(
    "<- ", variable_GB, "$rname %>% forcats::fct_drop()"
)
operation <- makeOperation(paste0(variable_GB, "$rname"), command)
evaluateOperation(operation)


#  Drop rows that are not chr1-19, chrX ---------------------------------------
chromosomes <- c(paste0("chr", 1:19), "chrX")
command <- paste0(
    "<- ", variable_GB, " %>% ", "filter(., rname %in% chromosomes)"
)
operation <- makeOperation(variable_GB, command)
evaluateOperation(operation)


#  Create and append pos_end, mpos_end columns --------------------------------
command <- paste0("<- ", variable_GB, "$pos + 49")
operation <- makeOperation(paste0(variable_GB, "$pos_end"), command)
evaluateOperation(operation)

command <- paste0("<- ", variable_GB, "$mpos + 49")
operation <- makeOperation(paste0(variable_GB, "$mpos_end"), command)
evaluateOperation(operation)

command <- paste0(
    "<- ", variable_GB, " %>% ",
        "dplyr::relocate(pos_end, .after = pos) %>% ",
        "dplyr::relocate(mpos_end, .after = mpos)"
)
operation <- makeOperation(variable_GB, command)
evaluateOperation(operation)


#  Create a temporary tags such as "coordinate" to distinguish entries --------
command <- paste0(
    "<- ", variable_GB, " %>% ",
        "tidyr::unite(",
            "coordinate, ",
            "c(\"qname\", \"rname\", \"pos\"), ",
            # "c(\"qname\", \"rname\", \"pos\", \"pos_end\"), ",
            "sep = \"_\", ",
            "remove = FALSE",
        ")", " %>% ",
            "dplyr::relocate(c(qname, coordinate), .after = qual)"
)
operation <- makeOperation(variable_GB, command)
evaluateOperation(operation)

#  Set up $qpos, a variable needed for sorting
command <- paste0(
    "<- paste0(",
        variable_GB, "$qname, ", "'_', ",
        variable_GB, "$pos",
    ")"
)
operation <- makeOperation(paste0(variable_GB, "$qpos"), command)
evaluateOperation(operation)

#  Set up $qmpos, a variable needed for sorting
command <- paste0(
    "<- paste0(",
        variable_GB, "$qname, ", "'_', ",
        variable_GB, "$mpos",
    ")"
)
operation <- makeOperation(paste0(variable_GB, "$qmpos"), command)
evaluateOperation(operation)

#  Set up $isize_abs, i.e., absolute insert-size values, for sorting
command <- paste0(
    "<- ", variable_GB, "$isize", " %>% ",
        "abs()", " %>% ",
        "as.numeric()"
)
operation <- makeOperation(paste0(variable_GB, "$isize_abs"), command)
evaluateOperation(operation)


#  Explicitly label the assignments from Giancarlo's pipeline
GB.alt.CAST$assignment <- "GB.CAST"
GB.ref.129S1$assignment <- "GB.129S1"
GB.ambig$assignment <- "GB.Ambiguous"
GB.contra$assignment <- "GB.Contra"


#  Reorganize the tibbles for easy reading
command <- paste0(
    "<- ", variable_GB, " %>% ",
        "dplyr::relocate(mpos, .after = pos)", " %>% ",
        "dplyr::relocate(c(qpos, qmpos), .after = qname)", " %>% ",
        "dplyr::relocate(isize_abs, .after = isize)", " %>% ",
        "dplyr::relocate(assignment, .after = isize_abs)"
)
operation <- makeOperation(variable_GB, command)
evaluateOperation(operation)


#  Check on mate pairing; generate vectors and tables for mate pairing --------
command <- paste0(
    "<- testMatesPaired(",
        "pos = ", variable_GB, "$qpos, mpos = ", variable_GB, "$qmpos",
    ")"
)  # VMP: "vector mates paired"
operation <- makeOperation(paste0("VMP.", variable_GB), command)
evaluateOperation(operation)

#  Convert NAs: assign FALSE to NAs
operation <- paste0(
    "VMP.", variable_GB[!str_detect(variable_GB, "unmated")], "[",
        "is.na(", "VMP.", variable_GB[!str_detect(variable_GB, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

# #  Check: number of FALSE should equal number of TRUE
# VMP.GB.alt.CAST %>% table()  # Yes
# VMP.GB.ambig %>% table()  # Yes
# VMP.GB.contra %>% table()  # Yes
# VMP.GB.ref.129S1 %>% table()  # Yes

#  Are `TRUE`, `FALSE` *not* alternating? (They should be)
#+
#+ Generate vectors and tables
command <- paste0(
    "<- testLogicalAlternating(logical = ", "VMP.", variable_GB, ")"
)  # VA: "vector alternating"
operation <- makeOperation(paste0("VA.", variable_GB), command)
evaluateOperation(operation)

#  Convert NAs: assign FALSE to NAs
operation <- paste0(
    "VA.", variable_GB[!str_detect(variable_GB, "unmated")], "[",
        "is.na(", "VA.", variable_GB[!str_detect(variable_GB, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

# #  Check (want to see FALSE and not TRUE here)
# VA.GB.alt.CAST %>% table()
# VA.GB.ambig %>% table()
# VA.GB.contra %>% table()
# VA.GB.ref.129S1 %>% table()

#  Create tbl columns denoting lack of alternation, i.e., $discrepancy
command <- paste0("<- ", "VA.", variable_GB)
operation <- makeOperation(paste0(variable_GB, "$discrepancy"), command)
evaluateOperation(operation)

#  Clean up
operation <- paste0(
    "rm(", paste0("VMP.", variable_GB), ", ", paste0("VA.", variable_GB), ")"
)
evaluateOperation(operation)

#  Concatenate even row to preceding odd row ----------------------------------
#+
#+ Get the mate pairs into one row each instead of two separate rows:
#+ "uniline" variables
variable_uniline <- paste0("uniline.", variable_GB)
mapply(
    assign, variable_uniline, "", MoreArgs = list(envir = parent.frame())
)
uniline.GB.alt.CAST <- mungeMatesIntoOneRow(variable_GB[1], "pos.x")
uniline.GB.ambig <- mungeMatesIntoOneRow(variable_GB[2], "pos.x")
uniline.GB.contra <- mungeMatesIntoOneRow(variable_GB[3], "pos.x")
uniline.GB.ref.129S1 <- mungeMatesIntoOneRow(variable_GB[4], "pos.x")


#  Join the related tibbles together ------------------------------------------
uniline.joint.GB <- dplyr::bind_rows(
    uniline.GB.alt.CAST,
    uniline.GB.ref.129S1,
    uniline.GB.ambig,
    uniline.GB.contra
)
colnames(uniline.joint.GB) <- colnames(uniline.joint.GB) %>%
    stringr::str_replace("\\.odd", "\\.GB.odd")
colnames(uniline.joint.GB) <- colnames(uniline.joint.GB) %>%
    stringr::str_replace("\\.even", "\\.GB.even")

uniline.joint.GB <- uniline.joint.GB %>%
    dplyr::arrange(rname.GB.odd, pos.GB.odd, mpos.GB.odd)


#  Deconcatenate uniline.joint.GB ---------------------------------------------
joint.GB <- mungeMatesIntoTwoRows(
    "uniline.joint.GB",
    "rname.GB.odd",
    "pos.GB.odd",
    "mpos.GB.odd",
    "assignment.GB"
)
colnames(joint.GB)[14] <- "groupid.initial.GB"
colnames(joint.GB)[13] <- "groupid.GB"


#  Clean up, save environment, etc. -------------------------------------------
setwd(paste0(path.1, "/", path.3))

operation <- paste0(
    "rm(", variable_GB, ", ", variable_uniline, ")"
)
evaluateOperation(operation)

rm(
    chromosome, chromosomes, command, file,
    operation, variable_GB, variable_uniline
)

#  Save the image
save.image(stringr::str_remove(script, ".R") %>% paste0(., ".Rdata"))


#  Script is completed; clean up environment ----------------------------------
rm(list = ls())
