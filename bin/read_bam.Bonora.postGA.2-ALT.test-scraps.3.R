#!/usr/bin/env Rscript

library(parallel)
library(Rsamtools)
library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)


#  Set up work directory (locations TB∆) --------------------------------------
directory_user <- "/Users/kalavattam"
directory_base <- "Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do"
directory_work <- "get_unique_fragments/Bonora"
setwd(
    paste0(directory_user, "/", directory_base, "/", directory_work)
)
rm(directory_user, directory_base, directory_work)


#  Set up functions -----------------------------------------------------------
`%notin%` <- Negate(`%in%`)


evaluateOperation <- function(operation = operation) {
    return(eval(parse(text = operation), envir = .GlobalEnv))
}


makeOperation <- function(variable, command) {
    operation <- paste(variable, command) #%>% paste(., collapse = "; ")
    return(operation)
}


#  Assign variables necessary for subsequent commands -------------------------
# chromosome <- "chr1"
chromosome <- "chrX"
file <- list.files(pattern = paste0("\\", chromosome, ".bam$")) %>%
    stringr::str_subset("^dedup\\.", negate = TRUE) %>%
    stringr::str_subset("^mm10\\.", negate = TRUE)
variable <- file %>%
    stringr::str_split(., "-") %>%
    lapply(., `[[`, 1) %>%
    unlist() %>%
    paste0("dedup.", .)
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)


#  Assign .bam information as list --------------------------------------------

#  What to query from .bam files
map_info <- c(
    "qname", "flag", "rname", "pos", "mapq", "cigar", 
    "mrnm", "mpos", "isize", "seq", "qual"
)
tag_info <- "AS"
map_params <- Rsamtools::ScanBamParam(what = map_info, tag = tag_info)

command <- paste0(
    "<- ", variable, " %>% ",
        "Rsamtools::scanBam(., param = map_params)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

rm(map_info, tag_info, map_params)


#  Convert .bam information from list to dataframe to tibble ------------------
command <- paste0("<- ", variable, " %>% ", "as.data.frame() %>% as_tibble()")
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Removal of Picard MarkDuplicates flags -------------------------------------
dedup.129S1$flag %>% table()
duplicate_flags <- c(1107, 1123, 1171, 1187)
command <- paste0(
    "<- ", variable, " %>% ", "dplyr::filter(flag %notin% duplicate_flags)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


# #  Reorder rname factor levels ------------------------------------------------
# command <- paste0("<- forcats::fct_relevel(", variable, "$rname, chromosome)")
# operation <- makeOperation(paste0(variable, "$rname"), command)
# evaluateOperation(operation)
# 
# 
# #  Drop unused rname factor levels --------------------------------------------
# command <- paste0("<- ", variable, "$rname %>% forcats::fct_drop()")
# operation <- makeOperation(paste0(variable, "$rname"), command)
# evaluateOperation(operation)
# 
# 
# #  Drop rows that are not chr1-19, chrX ---------------------------------------
# chromosomes <- c(paste0("chr", 1:19), "chrX")
# command <- paste0(
#     "<- ", variable, " %>% ", "filter(., rname %in% chromosomes)"
#     )
# operation <- makeOperation(variable, command)
# evaluateOperation(operation)
# 
# 
# #  Create and append pos_end, mpos_end columns --------------------------------
# command <- paste0("<- ", variable, "$pos + 49")
# operation <- makeOperation(paste0(variable, "$pos_end"), command)
# evaluateOperation(operation)
# 
# command <- paste0("<- ", variable, "$mpos + 49")
# operation <- makeOperation(paste0(variable, "$mpos_end"), command)
# evaluateOperation(operation)
# 
# command <- paste0(
#     "<- ", variable, " %>% ",
#         "dplyr::relocate(pos_end, .after = pos) %>% ",
#         "dplyr::relocate(mpos_end, .after = mpos)"
# )
# operation <- makeOperation(variable, command)
# evaluateOperation(operation)


#  Order columns by coordinate, tag.AS, and mapq ------------------------------
command <- paste0(
    "<- ", variable, "[",
        "order(",
            variable, "$qname, ",
            variable, "$rname, ",
            variable, "$pos, ",
            variable, "$mpos, -",
            variable, "$AS, -",
            variable, "$mapq",
        "), ",
    "]"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


# #  Based on coordinate value, identify if a given row is a duplicate ----------
# command <- paste0(
#     "<- ave(",
#         variable, "$coordinate, ",
#         variable, "$coordinate, ",
#         "FUN = length",
#     ") > 1L"
# )
# operation <- makeOperation(paste0(variable, "$duplicated"), command)
# evaluateOperation(operation)


#  Create a "temporary tag" for use in distinguishing entries
#+ 
#+ This is necessary to filter out duplicates and join liftOver data to initial
#+ tibbles; it's named "old_coordinate", i.e., the coordinate prior to liftOver
command <- paste0(
    "<- ", variable, " %>% ",
        "tidyr::unite(",
            "old_coordinate, ",
            "c(\"qname\", \"rname\", \"pos\"), ",
            "sep = \"_\", ",
            "remove = FALSE",
        ")"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Filter out all rows without unique coordinates -----------------------------
command <- paste0(
    "<- ", variable, " %>% ",
        "dplyr::distinct(old_coordinate, .keep_all = TRUE)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)
# dedup.129S1 %>% nrow()  # [1] 494694 (chr1); [1] 274025 (chrX, w/o rm MarkDuplicates); [1] 243239 (chrX, w/rm MarkDuplicates)
# dedup.CAST %>% nrow()  # [1] 492644 (chr1); [1] 281850 (chrX, w/o rm MarkDuplicates); [1] 249783 (chrX, w/rm MarkDuplicates)
# dedup.mm10 %>% nrow()  # [1] 491416 (chr1); [1] 283022 (chrX, w/o rm MarkDuplicates); [1] 252171 (chrX, w/rm MarkDuplicates)


#  Remove the "temporary tag" -------------------------------------------------
command <- paste0(
    "<- ", variable, " %>% ",
        "dplyr::select(-old_coordinate, -AS)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Write out .sam files -------------------------------------------------------
for (i in 1:length(variable)) {
    readr::write_tsv(
        x = eval(parse(text = variable[i])),
        file = paste0(variable[i], ".", chromosome, ".sam"),
        col_names = FALSE
    )
}


#  Add headers to the .sam files ----------------------------------------------
for (i in 1:length(variable)) {
shell_code <- paste0(
"
samtools view -H '", file[i], "' >> header.txt

cat ", variable[i], ".", chromosome, ".sam >> header.txt

mv header.txt ", variable[i], ".", chromosome, ".sam
"
)

system(shell_code)
}


#  Convert the .sam files to .bam files ---------------------------------------
for (i in 1:length(variable)) {
shell_code <- paste0(
"
samtools view -S -b ", variable[i], ".", chromosome, ".sam > ", variable[i], ".", chromosome, ".bam
"
)

system(shell_code)
}


#  Sort the .bam files --------------------------------------------------------
for (i in 1:length(variable)) {
shell_code <- paste0(
"
samtools sort ", variable[i], ".", chromosome, ".bam -o sorted.bam

mv sorted.bam ", variable[i], ".", chromosome, ".bam
"
)

system(shell_code)
}


#  Index the .bam files -------------------------------------------------------
for (i in 1:length(variable)) {
    shell_code <- paste0(
        "samtools index ", variable[i], ".", chromosome, ".bam"
    )
    
    system(shell_code)
}


#  Clean up: Remove the .sam files --------------------------------------------
for (i in 1:length(variable)) {
    shell_code <- paste0("rm ", variable[i], ".", chromosome, ".sam")

    system(shell_code)
}

rm(list = ls())
