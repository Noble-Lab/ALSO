#!/usr/bin/env Rscript

library(ggplot2)
library(stringr)
library(tidyverse)
library(scales)

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)


#  Set up work directory (locations TB∆) --------------------------------------
directory_user <- "/Users/kalavattam"
directory_base <- "Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do"
directory_work <- "get_unique_fragments/Bonora"

setwd(
    paste0(directory_user, "/", directory_base, "/", directory_work)
)
rm(directory_user, directory_base, directory_work)


#  Set up functions -----------------------------------------------------------
returnObjectName <- function(object) {
    return(deparse(substitute(object)))
}


#  Having run part 3 first, load part-3 .Rdata into environment ---------------
load("test.PE-processing.part-3.Rdata")
#  Doing so loads in appropriate functions, variables, etc.

#  Clean up environment
# rm(
#     uniline.a.ambiguous.dedup.129S1,
#     uniline.a.ambiguous.dedup.CAST,
#     uniline.a.ambiguous.dedup.mm10,
#     uniline.a.mated.dedup.129S1,
#     uniline.a.mated.dedup.CAST,
#     uniline.a.mated.dedup.mm10
# )

rm(
    adjust.a.ambiguous.dedup.129S1,
    adjust.a.ambiguous.dedup.CAST,
    adjust.a.ambiguous.dedup.mm10,
    adjust.a.mated.dedup.129S1,
    adjust.a.mated.dedup.CAST,
    adjust.a.mated.dedup.mm10,
    adjust.a.unmated.dedup.129S1,
    adjust.a.unmated.dedup.CAST,
    adjust.a.unmated.dedup.mm10
)


#  Script name ----------------------------------------------------------------
script <- "test.PE-processing.part-4.R"


#  Create $coordinate columns from pertinent liftOver values ------------------
variable <- paste0("tbl.", c("129S1", "CAST", "mm10"))
command <- paste0(
    "<- ", stringr::str_subset(variable, "mm10", negate = TRUE), " %>% ",
        "tidyr::unite(",
            "coordinate, ",
            "c(\"qname\", \"lO_rname\", \"lO_pos\"), ",
            "sep = \"_\", ",
            "remove = FALSE",
        ")", " %>% ",
            "dplyr::relocate(", "coordinate, ", ".before = qname", ")"
)
operation <- makeOperation(
    stringr::str_subset(variable, "mm10", negate = TRUE), command
)
evaluateOperation(operation)

#  Create $coordinate columns from pertinent values
command <- paste0(
    "<- ", stringr::str_subset(variable, "mm10", negate = FALSE), " %>% ",
        "tidyr::unite(",
            "coordinate, ",
            "c(\"qname\", \"rname\", \"pos\"), ",
            "sep = \"_\", ",
            "remove = FALSE",
        ")", " %>% ",
            "dplyr::relocate(", "coordinate, ", ".before = qname", ")"
)
operation <- makeOperation(
    stringr::str_subset(variable, "mm10", negate = FALSE), command
)
evaluateOperation(operation)


# #  Check on duplicated() for $coordinate --------------------------------------
# duplicated(tbl.129S1$coordinate) %>% table()
# duplicated(tbl.CAST$coordinate) %>% table()
# duplicated(tbl.mm10$coordinate) %>% table()
# 
# tbl.129S1[is.na(tbl.129S1$lO_rname), ] %>% nrow()
# test.129S1 <- tbl.129S1[!is.na(tbl.129S1$lO_rname), ]
# 
# tbl.CAST[is.na(tbl.CAST$lO_rname), ] %>% nrow()
# test.CAST <- tbl.CAST[!is.na(tbl.CAST$lO_rname), ]
# 
# duplicated(test.129S1$coordinate) %>% table()
# duplicated(test.CAST$coordinate) %>% table()
# 
# rm(test.129S1, test.CAST)


# #  Munging 129S1 --------------------------------------------------------------
# df <- tbl.129S1  # 225856
# df[is.na(df$flag), ]  # 0 x 31
# 
# df <- df[df$lO_reason == "liftOver successful", ]  # 223826
# df[is.na(df$flag), ]  # 0 x 31
# 
# df <- df[df$lO_rname == "chrX", ]  # 223452
# df[is.na(df$flag), ]  # 0 x 31
# 
# df <- df[!is.na(df$lO_mrnm), ]  # 222864
# df[is.na(df$flag), ]  # 0 x 31
# 
# df <- df[df$lO_mrnm == "chrX", ]  # 222860
# df[is.na(df$flag), ]  # 0 x 31
# 
# df[is.na(df$lO_pos), ]  # 0 x 31
# df[is.na(df$lO_mpos), ]  # 0 x 31
# df[is.na(df$lO_pos_end), ]  # 0 x 31
# df[is.na(df$lO_mpos_end), ]  # 0 x 31
# df[is.na(df$lO_reason), ]  # 0 x 31
# df[is.na(df$lO_rname), ]  # 0 x 31
# df[is.na(df$lO_mrnm), ]  # 0 x 31
# 
# df %>% sapply(., function(x) sum(is.na(x)))  # 0 for everything
# tbl.129S1 %>% sapply(., function(x) sum(is.na(x)))  # 2030 for lO_*


# #  Munging CAST ---------------------------------------------------------------
# df <- tbl.CAST  # 230790
# df[is.na(df$flag), ]  # 0 x 31
# 
# df <- df[df$lO_reason == "liftOver successful", ]  # 223401
# df[is.na(df$flag), ]  # 0 x 31
# 
# df <- df[df$lO_rname == "chrX", ]  # 219778
# df[is.na(df$flag), ]  # 0 x 31
# 
# df <- df[!is.na(df$lO_mrnm), ]  # 217530
# df[is.na(df$flag), ]  # 0 x 31
# 
# df <- df[df$lO_mrnm == "chrX", ]  # 217510
# df[is.na(df$flag), ]  # 0 x 31
# 
# df[is.na(df$lO_pos), ]  # 0 x 31
# df[is.na(df$lO_mpos), ]  # 0 x 31
# df[is.na(df$lO_pos_end), ]  # 0 x 31
# df[is.na(df$lO_mpos_end), ]  # 0 x 31
# df[is.na(df$lO_reason), ]  # 0 x 31
# df[is.na(df$lO_rname), ]  # 0 x 31
# df[is.na(df$lO_mrnm), ]  # 0 x 31
# 
# df %>% sapply(., function(x) sum(is.na(x)))  # 0 for everything
# tbl.CAST %>% sapply(., function(x) sum(is.na(x)))  # 7389 for lO_*


#  Systematic removal of NAs from tibbles subjected to liftOver ---------------
variable_lO <- variable %>% stringr::str_subset("mm10", negate = TRUE)

for (i in 1:length(variable_lO)) {
    #  Check NAs
    print(eval(parse(text = variable_lO[i])) %>%
              sapply(., function(x) sum(is.na(x))))
    
    #  Assign tibble to df, then munge
    df <- eval(parse(text = variable_lO[i]))
    df <- df[df$lO_reason == "liftOver successful", ]
    df <- df[df$lO_rname == "chrX", ]
    df <- df[!is.na(df$lO_mrnm), ]
    df <- df[df$lO_mrnm == "chrX", ]
    
    #  Assign df back to initial tibble
    command <- paste0("<- df")
    operation <- makeOperation(variable_lO[i], command)
    evaluateOperation(operation)
    
    #  Check NAs
    print(eval(parse(text = variable_lO[i])) %>%
              sapply(., function(x) sum(is.na(x))))
}
rm(df, i)


# #  Check to ensure mates are still paired -------------------------------------
# n_occur.129S1 <- tbl.129S1$groupid %>% table() %>% data.frame()  # 111430
# n_occur.129S1[n_occur.129S1$Freq == 2, ] %>% nrow()  # 111430
# n_occur.129S1[n_occur.129S1$Freq > 2, ] %>% nrow()  # 0
# n_occur.129S1[n_occur.129S1$Freq < 2, ] %>% nrow()  # 0
# 
# n_occur.CAST <- tbl.CAST$groupid %>% table() %>% data.frame()  # 108755
# n_occur.CAST[n_occur.CAST$Freq == 2, ] %>% nrow()  # 108755
# n_occur.CAST[n_occur.CAST$Freq > 2, ] %>% nrow()  # 0
# n_occur.CAST[n_occur.CAST$Freq < 2, ] %>% nrow()  # 0
# 
# n_occur.mm10 <- tbl.mm10$groupid %>% table() %>% data.frame()  # 117342
# n_occur.mm10[n_occur.mm10$Freq == 2, ] %>% nrow()  # 117342
# n_occur.mm10[n_occur.mm10$Freq > 2, ] %>% nrow()  # 0
# n_occur.mm10[n_occur.mm10$Freq < 2, ] %>% nrow()  # 0
# 
# rm(list = ls(pattern = "^n_occur"))


#  Check for duplicated $coordinate entries -----------------------------------

#  If present, then remove duplicates and remove any rows with less than two
#+ groupid entries

# #  Checks galore
# duplicated(tbl.129S1$coordinate) %>% table()
# duplicated(tbl.CAST$coordinate) %>% table()
# duplicated(tbl.mm10$coordinate) %>% table()
# 
# # -------------------------------------
# test.129S1 <- tbl.129S1[!duplicated(tbl.129S1$coordinate), ]  # 222707
# test.CAST <- tbl.CAST[!duplicated(tbl.CAST$coordinate), ]  # 217363
# test.mm10 <- tbl.mm10[!duplicated(tbl.mm10$coordinate), ]  # 234511
# 
# n_occur.129S1 <- test.129S1$groupid %>% table() %>% data.frame()  # 111430
# n_occur.CAST <- test.CAST$groupid %>% table() %>% data.frame()  # 108755
# n_occur.mm10 <- test.mm10$groupid %>% table() %>% data.frame()  # 117342
# 
# n_occur.129S1[n_occur.129S1$Freq == 2, ] %>% nrow()  # 111277
# n_occur.129S1[n_occur.129S1$Freq > 2, ] %>% nrow()  # 0
# n_occur.129S1[n_occur.129S1$Freq < 2, ] %>% nrow()  # 153
# 
# n_occur.CAST[n_occur.CAST$Freq == 2, ] %>% nrow()  # 108608
# n_occur.CAST[n_occur.CAST$Freq > 2, ] %>% nrow()  # 0
# n_occur.CAST[n_occur.CAST$Freq < 2, ] %>% nrow()  # 147
# 
# n_occur.mm10[n_occur.mm10$Freq == 2, ] %>% nrow()  # 117169
# n_occur.mm10[n_occur.mm10$Freq > 2, ] %>% nrow()  # 0
# n_occur.mm10[n_occur.mm10$Freq < 2, ] %>% nrow()  # 173
# 
# test.129S1[test.129S1$groupid %in% n_occur.129S1$.[n_occur.129S1$Freq < 2], ]  # 153 x 31
# test.CAST[test.CAST$groupid %in% n_occur.CAST$.[n_occur.CAST$Freq < 2], ]  # 147 x 31
# test.mm10[test.mm10$groupid %in% n_occur.mm10$.[n_occur.mm10$Freq < 2], ]  # 173 x 24
# 
# # -------------------------------------
# test.129S1 <- test.129S1[!test.129S1$groupid %in% n_occur.129S1$.[n_occur.129S1$Freq < 2], ]  # 222554 x 31
# test.CAST <- test.CAST[!test.CAST$groupid %in% n_occur.CAST$.[n_occur.CAST$Freq < 2], ]  # 217216 x 31
# test.mm10 <- test.mm10[!test.mm10$groupid %in% n_occur.mm10$.[n_occur.mm10$Freq < 2], ]  # 234338 x 24
# 
# n_occur.129S1 <- test.129S1$groupid %>% table() %>% data.frame()  # 111277
# n_occur.CAST <- test.CAST$groupid %>% table() %>% data.frame()  # 108608
# n_occur.mm10 <- test.mm10$groupid %>% table() %>% data.frame()  # 117169
# 
# n_occur.129S1[n_occur.129S1$Freq == 2, ] %>% nrow()  # 111277
# n_occur.129S1[n_occur.129S1$Freq > 2, ] %>% nrow()  # 0
# n_occur.129S1[n_occur.129S1$Freq < 2, ] %>% nrow()  # 0
# 
# n_occur.CAST[n_occur.CAST$Freq == 2, ] %>% nrow()  # 108608
# n_occur.CAST[n_occur.CAST$Freq > 2, ] %>% nrow()  # 0
# n_occur.CAST[n_occur.CAST$Freq < 2, ] %>% nrow()  # 0
# 
# n_occur.mm10[n_occur.mm10$Freq == 2, ] %>% nrow()  # 117169
# n_occur.mm10[n_occur.mm10$Freq > 2, ] %>% nrow()  # 0
# n_occur.mm10[n_occur.mm10$Freq < 2, ] %>% nrow()  # 0

#  Do the actual work
tbl.129S1 <- tbl.129S1[!duplicated(tbl.129S1$coordinate), ]  # 222707
tbl.CAST <- tbl.CAST[!duplicated(tbl.CAST$coordinate), ]  # 217363
tbl.mm10 <- tbl.mm10[!duplicated(tbl.mm10$coordinate), ]  # 234511

n_occur.129S1 <- tbl.129S1$groupid %>% table() %>% data.frame()  # 111430
n_occur.CAST <- tbl.CAST$groupid %>% table() %>% data.frame()  # 108755
n_occur.mm10 <- tbl.mm10$groupid %>% table() %>% data.frame()  # 117342

tbl.129S1 <- tbl.129S1[!tbl.129S1$groupid %in% n_occur.129S1$.[n_occur.129S1$Freq < 2], ]  # 222554 x 31
tbl.CAST <- tbl.CAST[!tbl.CAST$groupid %in% n_occur.CAST$.[n_occur.CAST$Freq < 2], ]  # 217216 x 31
tbl.mm10 <- tbl.mm10[!tbl.mm10$groupid %in% n_occur.mm10$.[n_occur.mm10$Freq < 2], ]  # 234338 x 24

n_occur.129S1 <- tbl.129S1$groupid %>% table() %>% data.frame()  # 111277
n_occur.CAST <- tbl.CAST$groupid %>% table() %>% data.frame()  # 108608
n_occur.mm10 <- tbl.mm10$groupid %>% table() %>% data.frame()  # 117169

rm(list = ls(pattern = "^n_occur."))


#  Concatenate even row to preceding odd row in  m.* variables ----------------
#+
#+ Get the mate pairs into one row each instead of two separate rows:
#+ "uniline" variables

#  129S1, CAST
variable_uniline <- paste0("uniline.", variable_lO)
mapply(
    assign, variable_uniline, "", MoreArgs = list(envir = parent.frame())
)
for (i in 1:length(variable_lO)) {
    df <- eval(parse(text = variable_lO[i]))

    odd.seq <- seq(1, nrow(df), 2)
    even.seq <- (seq(1, nrow(df), 2)) + 1
    
    odd <- df[odd.seq, ]
    even <- df[even.seq, ]
    
    df <- full_join(odd, even, by = "groupid") %>%
        dplyr::arrange(lO_pos.x) %>%  # Sort by pos.x
        dplyr::rename(groupid.x = groupid)  # new_name = old_name
    df$groupid.y <- df$groupid.x
    df <- df %>% dplyr::relocate(groupid.y, .after = qmpos.y)
    
    colnames(df) <- str_replace_all(colnames(df), "\\.x", "\\.odd")
    colnames(df) <- str_replace_all(colnames(df), "\\.y", "\\.even")
    
    command <- paste0("<- ", "df")
    operation <- makeOperation(variable_uniline[i], command)
    evaluateOperation(operation)
}

#  mm10
variable_uniline <- paste0("uniline.", stringr::str_subset(variable, "mm10"))
mapply(
    assign, variable_uniline, "", MoreArgs = list(envir = parent.frame())
)
for (i in 1:length(stringr::str_subset(variable, "mm10"))) {
    df <- eval(parse(text = stringr::str_subset(variable, "mm10")[i]))
    
    odd.seq <- seq(1, nrow(df), 2)
    even.seq <- (seq(1, nrow(df), 2)) + 1
    
    odd <- df[odd.seq, ]
    even <- df[even.seq, ]
    
    df <- full_join(odd, even, by = "groupid") %>%
        dplyr::arrange(pos.x) %>%  # Sort by pos.x
        dplyr::rename(groupid.x = groupid)  # new_name = old_name
    df$groupid.y <- df$groupid.x
    df <- df %>% dplyr::relocate(groupid.y, .after = qmpos.y)
    
    colnames(df) <- str_replace_all(colnames(df), "\\.x", "\\.odd")
    colnames(df) <- str_replace_all(colnames(df), "\\.y", "\\.even")
    
    command <- paste0("<- ", "df")
    operation <- makeOperation(variable_uniline[i], command)
    evaluateOperation(operation)
}

rm(df, i, even, even.seq, odd, odd.seq)


#  Output paired, sorted .sam files -------------------------------------------

#  Sort and deconcatenate uniline tibbles for 129, CAST
variable_uniline <- paste0("uniline.", variable) %>%
    stringr::str_subset("mm10", negate = TRUE)

for (i in 1:length(variable_uniline)) {
    df <- eval(parse(text = variable_uniline[i]))
    
    #  Sort the uniline tibbles by rname, pos, mpos
    df <- df %>% dplyr::arrange(lO_rname.odd, lO_pos.odd, lO_mpos.odd)
    
    #  Split the uniline tibbles by odd or even status
    odd <- df[stringr::str_subset(colnames(df), "\\.odd")]
    even <- df[stringr::str_subset(colnames(df), "\\.even")]
    
    #  Strip suffixes from column names
    colnames(odd) <- str_replace_all(colnames(odd), "\\.odd", "")
    colnames(even) <- str_replace_all(colnames(even), "\\.even", "")
    
    #  Odd tibble is tibble #1, even tibble is tibble #2
    odd$tibble <- "1"
    even$tibble <- "2"
    
    #  Again, interleave the odd/even rows, then arrange them by group ID and
    #+ tibble number; save results to initial "tbl.*" variables
    command <- paste0(
        "<- odd %>%",
            "dplyr::mutate(groupid = row_number()) %>% ",
            "dplyr::bind_rows(even %>% mutate(groupid = row_number())) %>% ",
            "dplyr::arrange(groupid, tibble) %>% ",
            "dplyr::relocate(groupid, .after = qmpos) %>% ",
            "dplyr::select(-tibble)"
    )
    operation <- makeOperation(variable[i], command)
    evaluateOperation(operation)
    
    #  Again, check to make sure that there are no more than two entries per
    #+ group ID
    n_occur <- data.frame(table(eval(parse(text = variable[i]))$groupid))
    print(n_occur[n_occur$Freq > 2, ])
    
    #  Clean up
    rm(df, n_occur, even, odd)
}
rm(i)

#  Sort and deconcatenate uniline tibbles for mm10
variable_uniline <- paste0("uniline.", variable) %>%
    stringr::str_subset("mm10", negate = FALSE)

for (i in 1:length(variable_uniline)) {
    df <- eval(parse(text = variable_uniline[i]))
    
    #  Sort the uniline tibbles by rname, pos, mpos
    df <- df %>% dplyr::arrange(rname.odd, pos.odd, mpos.odd)
    
    #  Split the uniline tibbles by odd or even status
    odd <- df[stringr::str_subset(colnames(df), "\\.odd")]
    even <- df[stringr::str_subset(colnames(df), "\\.even")]
    
    #  Strip suffixes from column names
    colnames(odd) <- str_replace_all(colnames(odd), "\\.odd", "")
    colnames(even) <- str_replace_all(colnames(even), "\\.even", "")
    
    #  Odd tibble is tibble #1, even tibble is tibble #2
    odd$tibble <- "1"
    even$tibble <- "2"
    
    #  Again, interleave the odd/even rows, then arrange them by group ID and
    #+ tibble number; save results to initial "tbl.*" variables
    command <- paste0(
        "<- odd %>%",
        "dplyr::mutate(groupid = row_number()) %>% ",
        "dplyr::bind_rows(even %>% mutate(groupid = row_number())) %>% ",
        "dplyr::arrange(groupid, tibble) %>% ",
        "dplyr::relocate(groupid, .after = qmpos) %>% ",
        "dplyr::select(-tibble)"
    )
    operation <- makeOperation(variable[i], command)
    evaluateOperation(operation)
    
    #  Again, check to make sure that there are no more than two entries per
    #+ group ID
    n_occur <- data.frame(table(eval(parse(text = variable[i]))$groupid))
    print(n_occur[n_occur$Freq > 2, ])
    
    #  Clean up
    rm(df, n_occur, even, odd)
}
rm(i)

#  Create sam.*.tbl for 1291, CAST
variable_sam <- paste0(
    "sam.",
    variable %>%
        stringr::str_subset("mm10", negate = TRUE) %>%
        stringr::str_remove("tbl.")
)
command <- paste0(
    "<- ", variable, " %>% ",
        "dplyr::select(",
            "qname, flag, rname, pos, mapq, cigar, ",
            "mrnm, mpos, isize, seq, qual",
        ")"
)
operation <- makeOperation(variable_sam, command)
evaluateOperation(operation)

#  Create sam.*.tbl for mm10
variable_sam <- paste0(
    "sam.",
    variable %>%
        stringr::str_subset("mm10", negate = FALSE) %>%
        stringr::str_remove("tbl.")
)
command <- paste0(
    "<- ", variable, " %>% ",
    "dplyr::select(",
        "qname, flag, rname, pos, mapq, cigar, ",
        "mrnm, mpos, isize, seq, qual",
    ")"
)
operation <- makeOperation(variable_sam, command)
evaluateOperation(operation)


#  Write out .sam files -------------------------------------------------------
chromosome <- "chrX"
variable_sam <- paste0("sam.", variable %>% stringr::str_remove("tbl."))

file_sam <- vector(mode = "character", length = length(variable_sam))
for (i in 1:length(variable_sam)) {
    file_sam[i] <- paste0(
        "processed", ".",
        "mate-paired", ".",
        stringr::str_remove(variable_sam[i], "sam."), ".",
        chromosome, ".",
        "sam"
    )
}
rm(i)


#  Write out .sam files -------------------------------------------------------
for (i in 1:length(variable_sam)) {
    readr::write_tsv(
        x = eval(parse(text = variable_sam[i])),
        file = file_sam[i],
        col_names = FALSE
    )
}


#  Add headers to the .sam files ----------------------------------------------
file <- list.files(pattern = paste0("\\", chromosome, ".rmdup.bam$"))
file <- file %>% stringr::str_subset("mm10\\.", negate = TRUE)
variable_bam_initial <- file %>%
    strsplit(., "-") %>%
    lapply(., `[[`, 1) %>% unlist() %>%
    paste0("dedup.", .)
mapply(
    assign, variable_bam_initial, file, MoreArgs = list(envir = parent.frame())
)

for (i in 1:length(variable_sam)) {
shell_code <- paste0(
"
samtools view -H '", eval(parse(text = variable_bam_initial[i])), "' >> header.txt

cat ", file_sam[i], " >> header.txt

mv header.txt ", file_sam[i]
)

system(shell_code)
}
rm(i)


#  Convert the .sam files to .bam files ---------------------------------------
for (i in 1:length(variable_sam)) {
shell_code <- paste0(
"
samtools view -S -b ", file_sam[i], " > ", stringr::str_replace(file_sam[i], ".sam", ".bam")
)

system(shell_code)
}


#  Sort the .bam files --------------------------------------------------------
for (i in 1:length(variable_sam)) {
shell_code <- paste0(
"
samtools sort ", stringr::str_replace(file_sam[i], ".sam", ".bam"), " -o sorted.bam

mv sorted.bam ", stringr::str_replace(file_sam[i], ".sam", ".sorted.bam")
)

system(shell_code)
}


#  Index the .bam files -------------------------------------------------------
for (i in 1:length(variable_sam)) {
    shell_code <- paste0(
        "samtools index ", stringr::str_replace(file_sam[i], ".sam", ".sorted.bam")
    )
    
    system(shell_code)
}
rm(i)

#  Save the .sam files as .rds files ------------------------------------------
operation <- paste0(
    "saveRDS(",
        variable_sam, ", ",
        "file = ", "\"", stringr::str_replace(file_sam, ".sam", ".bam.rds"), "\"",
    ")"
)
evaluateOperation(operation)

#  Clean up environment, save environment image -------------------------------
rm(
    chromosome,
    dedup.129S1, dedup.CAST, dedup.mm10,
    file, file_sam,
    sam.129S1, sam.CAST, sam.mm10,
    shell_code,
    tbl.129S1, tbl.CAST, tbl.mm10,
    variable, variable_bam_initial, variable_lO, variable_sam, variable_uniline
)

save.image(stringr::str_remove(script, ".R") %>% paste0(., ".Rdata"))


#  Script is completed; clean up environment ----------------------------------
rm(list = ls())


# #  Clean up: Remove the .sam files --------------------------------------------
# for (i in 1:length(variable_sam)) {
#     shell_code <- paste0("rm ", file_sam[i])
#     
#     system(shell_code)
# }
