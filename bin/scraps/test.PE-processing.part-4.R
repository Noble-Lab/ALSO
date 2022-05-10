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
mm10omit <- function(variable) {
    stringr::str_subset(variable, "mm10", negate = TRUE)
}


mm10only <- function(variable) {
    stringr::str_subset(variable, "mm10", negate = FALSE)
}


returnObjectName <- function(object) {
    return(deparse(substitute(object)))
}


#  Having run part 3 first, load part-3 .Rdata into environment ---------------
load("test.PE-processing.part-3.Rdata")
#  Doing so loads in appropriate functions, variables, etc.

#TODO Can remove these variables earlier: figure out when
rm(
    adjust.a.ambiguous.dedup.129,
    adjust.a.ambiguous.dedup.CAST,
    adjust.a.ambiguous.dedup.mm10,
    adjust.a.mated.dedup.129,
    adjust.a.mated.dedup.CAST,
    adjust.a.mated.dedup.mm10,
    adjust.a.unmated.dedup.129,
    adjust.a.unmated.dedup.CAST,
    adjust.a.unmated.dedup.mm10
)


#  Script name ----------------------------------------------------------------
script <- "test.PE-processing.part-4.R"


#TODO 1/2 Should "coordinate" be called "criteria" to be consistent with naming
#TODO 2/3 in "script-1.R" (but the change is not made in
#TODO 3/3 test.PE-processing.part-1.R)
#  Create $coordinate columns from pertinent liftOver values ------------------
variable <- paste0("tbl.", c("129", "CAST", "mm10"))

#  coordinate pre-lO: 129, CAST, mm10
command <- paste0(
    variable, " %>% ",
        "tidyr::unite(",
            "coordinate, ",
            "c(\"qname\", \"rname\", \"pos\"), ",
            "sep = \"_\", ",
            "remove = FALSE",
        ")", " %>% ",
            "dplyr::relocate(", "coordinate, ", ".before = criteria", ")"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  coordinate post-lO: 129, CAST
command <- paste0(
    mm10omit(variable), " %>% ",
        "tidyr::unite(",
            "lO_coordinate, ",
            "c(\"qname\", \"lO_rname\", \"lO_pos\"), ",
            "sep = \"_\", ",
            "remove = FALSE",
        ")", " %>% ",
            "dplyr::relocate(", "lO_coordinate, ", ".after = coordinate", ")"
)
operation <- makeOperation(mm10omit(variable), command)
evaluateOperation(operation)

#  coordinate post-lO: mm10
command <- paste0(
    mm10only(variable), " %>% ",
        "tidyr::unite(",
            "lO_coordinate, ",
            "c(\"qname\", \"rname\", \"pos\"), ",
            "sep = \"_\", ",
            "remove = FALSE",
        ")", " %>% ",
            "dplyr::relocate(", "lO_coordinate, ", ".after = coordinate", ")"
)
operation <- makeOperation(mm10only(variable), command)
evaluateOperation(operation)


<<<<<<< HEAD:bin/scraps/test.PE-processing.part-4.R
=======
#TODO Comment out the following section (see #TODO above)
>>>>>>> main:bin/test.PE-processing.part-4.R
#  Create $criteria columns from pertinent liftOver values ------------------
command <- paste0(
    mm10omit(variable), " %>% ",
        "tidyr::unite(",
            "lO_criteria, ",
            "c(\"qname\", \"flag\", \"lO_pos\", \"lO_mpos\"), ",
            "sep = \"_\", ",
            "remove = FALSE",
        ")", " %>% ",
            "dplyr::relocate(", "lO_criteria, ", ".after = criteria", ")"
)
operation <- makeOperation(mm10omit(variable), command)
evaluateOperation(operation)

#  mm10
command <- paste0(
    mm10only(variable), " %>% ",
        "tidyr::unite(",
            "lO_criteria, ",
            "c(\"qname\", \"flag\", \"pos\", \"mpos\"), ",
            "sep = \"_\", ",
            "remove = FALSE",
        ")", " %>% ",
            "dplyr::relocate(", "lO_criteria, ", ".after = criteria", ")"
)
operation <- makeOperation(mm10only(variable), command)
evaluateOperation(operation)


# #  Check on duplicated() for $lO_criteria -------------------------------------
# duplicated(tbl.129$lO_criteria) %>% table()  # 225194 662
# duplicated(tbl.CAST$lO_criteria) %>% table()  # 227518 3272
# duplicated(tbl.mm10$lO_criteria) %>% table()  # 234684 0
# 
# tbl.129[is.na(tbl.129$lO_rname), ] %>% nrow()  # 2030
# test.129 <- tbl.129[!is.na(tbl.129$lO_rname), ]  # 223826
# 
# tbl.CAST[is.na(tbl.CAST$lO_rname), ] %>% nrow()  # 7389
# test.CAST <- tbl.CAST[!is.na(tbl.CAST$lO_rname), ]  # 223401
# 
# duplicated(test.129$lO_criteria) %>% table()  # 223822 4
# duplicated(test.CAST$lO_criteria) %>% table()  # 223390 11
# 
# rm(test.129, test.CAST)


# #  Munging 129 --------------------------------------------------------------
# df <- tbl.129  # 225856
# df[is.na(df$flag), ]  # 0 x 34
# 
# df <- df[df$lO_reason == "liftOver successful", ]  # 223826
# df[is.na(df$flag), ]  # 0 x 34
# 
# df <- df[df$lO_rname == "chrX", ]  # 223452
# df[is.na(df$flag), ]  # 0 x 34
# 
# df <- df[!is.na(df$lO_mrnm), ]  # 222864
# df[is.na(df$flag), ]  # 0 x 34
# 
# df <- df[df$lO_mrnm == "chrX", ]  # 222860
# df[is.na(df$flag), ]  # 0 x 34
# 
# df[is.na(df$lO_pos), ]  # 0 x 34
# df[is.na(df$lO_mpos), ]  # 0 x 34
# df[is.na(df$lO_pos_end), ]  # 0 x 34
# df[is.na(df$lO_mpos_end), ]  # 0 x 34
# df[is.na(df$lO_reason), ]  # 0 x 34
# df[is.na(df$lO_rname), ]  # 0 x 34
# df[is.na(df$lO_mrnm), ]  # 0 x 34
# 
# df %>% sapply(., function(x) sum(is.na(x)))  # 0 for everything
# tbl.129 %>% sapply(., function(x) sum(is.na(x)))  # 2030 for lO_*
# rm(df)


# #  Munging CAST ---------------------------------------------------------------
# df <- tbl.CAST  # 230790
# df[is.na(df$flag), ]  # 0 x 34
# 
# df <- df[df$lO_reason == "liftOver successful", ]  # 223401
# df[is.na(df$flag), ]  # 0 x 34
# 
# df <- df[df$lO_rname == "chrX", ]  # 219778
# df[is.na(df$flag), ]  # 0 x 34
# 
# df <- df[!is.na(df$lO_mrnm), ]  # 217530
# df[is.na(df$flag), ]  # 0 x 34
# 
# df <- df[df$lO_mrnm == "chrX", ]  # 217510
# df[is.na(df$flag), ]  # 0 x 34
# 
# df[is.na(df$lO_pos), ]  # 0 x 34
# df[is.na(df$lO_mpos), ]  # 0 x 34
# df[is.na(df$lO_pos_end), ]  # 0 x 34
# df[is.na(df$lO_mpos_end), ]  # 0 x 34
# df[is.na(df$lO_reason), ]  # 0 x 34
# df[is.na(df$lO_rname), ]  # 0 x 34
# df[is.na(df$lO_mrnm), ]  # 0 x 34
# 
# df %>% sapply(., function(x) sum(is.na(x)))  # 0 for everything
# tbl.CAST %>% sapply(., function(x) sum(is.na(x)))  # 7389 for lO_*
# rm(df)


# #  Munging mm10 ---------------------------------------------------------------
# df <- tbl.mm10  # 234684
# df[is.na(df$flag), ]  # 0 x 27
# 
# df <- df[df$rname == "chrX", ]  # 234684
# df[is.na(df$flag), ]  # 0 x 27
# 
# df <- df[!is.na(df$mrnm), ]  # 234684
# df[is.na(df$flag), ]  # 0 x 27
# 
# df <- df[df$mrnm == "chrX", ]  # 234684
# df[is.na(df$flag), ]  # 0 x 27
# 
# df[is.na(df$pos), ]  # 0 x 27
# df[is.na(df$mpos), ]  # 0 x 27
# df[is.na(df$pos_end), ]  # 0 x 27
# df[is.na(df$mpos_end), ]  # 0 x 27
# df[is.na(df$rname), ]  # 0 x 27
# df[is.na(df$mrnm), ]  # 0 x 27
# 
# df %>% sapply(., function(x) sum(is.na(x)))  # 0 for everything
# tbl.mm10 %>% sapply(., function(x) sum(is.na(x)))  # 0 for everything
# rm(df)


#  Systematic removal of NAs from tibbles subjected to liftOver ---------------
variable_lO <- mm10omit(variable)

for (i in 1:length(variable_lO)) {
    #  Check NAs
    cat(paste0(
        "\n", "Working with ", variable_lO[i], ": ",
        "Table before removing NAs lO_* NAs..."
    ))
    print(eval(parse(text = variable_lO[i])) %>%
              sapply(., function(x) sum(is.na(x))))
    
    #  Assign tibble to df, then munge
    df <- eval(parse(text = variable_lO[i]))
    df <- df[df$lO_reason == "liftOver successful", ]
    df <- df[df$lO_rname == "chrX", ]
    df <- df[!is.na(df$lO_mrnm), ]
    df <- df[df$lO_mrnm == "chrX", ]
    
    #  Assign df back to initial tibble
    command <- paste0("df")
    operation <- makeOperation(variable_lO[i], command)
    evaluateOperation(operation)
    
    #  Check NAs
    cat(paste0(
        "\n", "Working with ", variable_lO[i], ": ",
        "Table after removing NAs lO_* NAs..."
    ))
    print(eval(parse(text = variable_lO[i])) %>%
              sapply(., function(x) sum(is.na(x))))
}
rm(df, i)


# #  Check to ensure mates are still paired -------------------------------------
# n_occur.129 <- tbl.129$groupid %>% table() %>% data.frame()  # 111430
# n_occur.129[n_occur.129$Freq == 2, ] %>% nrow()  # 111430
# n_occur.129[n_occur.129$Freq > 2, ] %>% nrow()  # 0
# n_occur.129[n_occur.129$Freq < 2, ] %>% nrow()  # 0
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


###############################################################################
#  Check for duplicated $lO_criteria entries -------------------------------------

#  If present, then remove duplicates and remove any rows with less than two
#+ groupid entries

# #  Checks galore
# duplicated(tbl.129$lO_criteria) %>% table(useNA = "ifany")  # FALSE 222860 TRUE 0
# duplicated(tbl.CAST$lO_criteria) %>% table(useNA = "ifany")  # FALSE 217510 TRUE 0
# duplicated(tbl.mm10$lO_criteria) %>% table(useNA = "ifany")  # FALSE 234684 TRUE 0
# 
# dup.129 <- tbl.129[duplicated(tbl.129$lO_criteria), ]  # 0
# dup.CAST <- tbl.CAST[duplicated(tbl.CAST$lO_criteria), ]  # 0
# dup.mm10 <- tbl.mm10[duplicated(tbl.mm10$lO_criteria), ]  # 0
# 
# dup.129 <- tbl.129 %>% dplyr::group_by(lO_criteria) %>% dplyr::filter(n() > 1)  # 0
# dup.CAST <- tbl.CAST %>% dplyr::group_by(lO_criteria) %>% dplyr::filter(n() > 1)  # 0
# dup.mm10 <- tbl.mm10 %>% dplyr::group_by(lO_criteria) %>% dplyr::filter(n() > 1)  # 0
# 
# dup.129$frequency <- table(dup.129$lO_criteria)[dup.129$lO_criteria]
# dup.CAST$frequency <- table(dup.CAST$lO_criteria)[dup.CAST$lO_criteria]
# dup.mm10$frequency <- table(dup.mm10$lO_criteria)[dup.mm10$lO_criteria]
# 
# test.129 <- tbl.129[!duplicated(tbl.129$lO_criteria), ]  # 222860
# test.CAST <- tbl.CAST[!duplicated(tbl.CAST$lO_criteria), ]  # 217510
# test.mm10 <- tbl.mm10[!duplicated(tbl.mm10$lO_criteria), ]  # 234684
# 
# n_occur.129 <- test.129$groupid %>% table() %>% data.frame()  # 111430
# n_occur.CAST <- test.CAST$groupid %>% table() %>% data.frame()  # 108755
# n_occur.mm10 <- test.mm10$groupid %>% table() %>% data.frame()  # 117342
# 
# n_occur.129[n_occur.129$Freq == 2, ] %>% nrow()  # 111430
# n_occur.129[n_occur.129$Freq > 2, ] %>% nrow()  # 0
# n_occur.129[n_occur.129$Freq < 2, ] %>% nrow()  # 0
# 
# n_occur.CAST[n_occur.CAST$Freq == 2, ] %>% nrow()  # 108755
# n_occur.CAST[n_occur.CAST$Freq > 2, ] %>% nrow()  # 0
# n_occur.CAST[n_occur.CAST$Freq < 2, ] %>% nrow()  # 0
# 
# n_occur.mm10[n_occur.mm10$Freq == 2, ] %>% nrow()  # 117342
# n_occur.mm10[n_occur.mm10$Freq > 2, ] %>% nrow()  # 0
# n_occur.mm10[n_occur.mm10$Freq < 2, ] %>% nrow()  # 0
# 
# test.129[test.129$groupid %in% n_occur.129$.[n_occur.129$Freq < 2], ]  # 0 x 34
# test.CAST[test.CAST$groupid %in% n_occur.CAST$.[n_occur.CAST$Freq < 2], ]  # 0 x 34
# test.mm10[test.mm10$groupid %in% n_occur.mm10$.[n_occur.mm10$Freq < 2], ]  # 0 x 27
# 
# test.129 <- test.129[!test.129$groupid %in% n_occur.129$.[n_occur.129$Freq < 2], ]  # 222860 x 34
# test.CAST <- test.CAST[!test.CAST$groupid %in% n_occur.CAST$.[n_occur.CAST$Freq < 2], ]  # 217510 x 34
# test.mm10 <- test.mm10[!test.mm10$groupid %in% n_occur.mm10$.[n_occur.mm10$Freq < 2], ]  # 234684 x 27
# 
# n_occur.129 <- test.129$groupid %>% table() %>% data.frame()  # 111430
# n_occur.CAST <- test.CAST$groupid %>% table() %>% data.frame()  # 108755
# n_occur.mm10 <- test.mm10$groupid %>% table() %>% data.frame()  # 117342
# 
# n_occur.129[n_occur.129$Freq == 2, ] %>% nrow()  # 111430
# n_occur.129[n_occur.129$Freq > 2, ] %>% nrow()  # 0
# n_occur.129[n_occur.129$Freq < 2, ] %>% nrow()  # 0
# 
# n_occur.CAST[n_occur.CAST$Freq == 2, ] %>% nrow()  # 108755
# n_occur.CAST[n_occur.CAST$Freq > 2, ] %>% nrow()  # 0
# n_occur.CAST[n_occur.CAST$Freq < 2, ] %>% nrow()  # 0
# 
# n_occur.mm10[n_occur.mm10$Freq == 2, ] %>% nrow()  # 117342
# n_occur.mm10[n_occur.mm10$Freq > 2, ] %>% nrow()  # 0
# n_occur.mm10[n_occur.mm10$Freq < 2, ] %>% nrow()  # 0

tbl.129 <- tbl.129[!duplicated(tbl.129$lO_criteria), ]  # 222860
tbl.CAST <- tbl.CAST[!duplicated(tbl.CAST$lO_criteria), ]  # 217510
tbl.mm10 <- tbl.mm10[!duplicated(tbl.mm10$lO_criteria), ]  # 234684

n_occur.129 <- tbl.129$groupid %>% table() %>% data.frame()  # 111430
n_occur.CAST <- tbl.CAST$groupid %>% table() %>% data.frame()  # 108755
n_occur.mm10 <- tbl.mm10$groupid %>% table() %>% data.frame()  # 117342

tbl.129 <- tbl.129[!tbl.129$groupid %in% n_occur.129$.[n_occur.129$Freq < 2], ]  # 222860 x 34
tbl.CAST <- tbl.CAST[!tbl.CAST$groupid %in% n_occur.CAST$.[n_occur.CAST$Freq < 2], ]  # 217510 x 34
tbl.mm10 <- tbl.mm10[!tbl.mm10$groupid %in% n_occur.mm10$.[n_occur.mm10$Freq < 2], ]  # 234684 x 27

n_occur.129 <- tbl.129$groupid %>% table() %>% data.frame()  # 111430
n_occur.CAST <- tbl.CAST$groupid %>% table() %>% data.frame()  # 108755
n_occur.mm10 <- tbl.mm10$groupid %>% table() %>% data.frame()  # 117342

rm(list = ls(pattern = "^n_occur."))


#  Concatenate even row to preceding odd row in  m.* variables ----------------
#+
#+ Get the mate pairs into one row each instead of two separate rows:
#+ "uniline" variables

#  129, CAST
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
    
    command <- paste0("df")
    operation <- makeOperation(variable_uniline[i], command)
    evaluateOperation(operation)
}

#  mm10
variable_uniline <- paste0("uniline.", stringr::str_subset(variable, "mm10"))
mapply(
    assign, variable_uniline, "", MoreArgs = list(envir = parent.frame())
)
for (i in 1:length(mm10only(variable))) {
    df <- eval(parse(text = mm10only(variable))[i])
    
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
    
    command <- paste0("df")
    operation <- makeOperation(variable_uniline[i], command)
    evaluateOperation(operation)
}

rm(df, i, even, even.seq, odd, odd.seq)


#  Output paired, sorted .sam files -------------------------------------------

#  Sort and deconcatenate uniline tibbles for 129, CAST
variable_uniline <- paste0("uniline.", variable) %>% mm10omit()

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
        "odd %>%",
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
variable_uniline <- paste0("uniline.", variable) %>% mm10only()

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
        "odd %>%",
            "dplyr::mutate(groupid = row_number()) %>% ",
            "dplyr::bind_rows(even %>% mutate(groupid = row_number())) %>% ",
            "dplyr::arrange(groupid, tibble) %>% ",
            "dplyr::relocate(groupid, .after = qmpos) %>% ",
            "dplyr::select(-tibble)"
    )
    operation <- makeOperation(mm10only(variable)[i], command)
    evaluateOperation(operation)
    
    #  Again, check to make sure that there are no more than two entries per
    #+ group ID
    n_occur <- data.frame(table(eval(parse(text = variable[i]))$groupid))
    print(n_occur[n_occur$Freq > 2, ])
    
    #  Clean up
    rm(df, n_occur, even, odd)
}
rm(i)

# #  Clean up unneeded variable
# operation <- paste0(
#     "rm(", "uniline.", variable, ")"
# )
# evaluateOperation(operation)

#  Add "proper" labels for MD and AS
command <- paste0(
    "paste0(\"MD:Z:\", ", variable, "$tag.MD)"
)
operation <- makeOperation(paste0(variable, "$tag.MD.proper"), command)
evaluateOperation(operation)

#  Add "proper" labels for AS
command <- paste0(
    "paste0(\"AS:i:\", ", variable, "$tag.AS)"
)
operation <- makeOperation(paste0(variable, "$tag.AS.proper"), command)
evaluateOperation(operation)

# #  Create sam.* tibbles for 1291, CAST: non-liftOver columns
# variable_sam <- paste0(
#     "sam.",
#     variable %>%
#         mm10omit %>%
#         stringr::str_remove("tbl.") #%>%
#         #paste0(., ".unlifted")
# )
# command <- paste0(
#     mm10omit(variable), " %>% ",
#         "dplyr::select(",
#             "qname, flag, rname, pos, mapq, cigar, ",
#             "mrnm, mpos, isize, seq, qual, tag.MD.proper, tag.AS.proper",
#         ")"
# )
# operation <- makeOperation(variable_sam, command)
# evaluateOperation(operation)

#  Create sam.* tibbles for 1291, CAST: liftOver columns
variable_sam <- paste0(
    "sam.",
    variable %>%
        mm10omit %>%
        stringr::str_remove("tbl.") #%>%
        # paste0(., ".lifted")
)
command <- paste0(
    mm10omit(variable), " %>% ",
        "dplyr::select(",
            "qname, flag, lO_rname, lO_pos, mapq, cigar, ",
            "lO_mrnm, lO_mpos, isize, seq, qual, tag.MD.proper, tag.AS.proper",
        ")"
)
operation <- makeOperation(variable_sam, command)
evaluateOperation(operation)

#  Create sam.*.tbl for mm10
variable_sam <- paste0(
    "sam.", variable %>%
        mm10only %>%
        stringr::str_remove("tbl.")
)
command <- paste0(
    mm10only(variable), " %>% ",
    "dplyr::select(",
        "qname, flag, rname, pos, mapq, cigar, ",
        "mrnm, mpos, isize, seq, qual, tag.MD.proper, tag.AS.proper",
    ")"
)
operation <- makeOperation(variable_sam, command)
evaluateOperation(operation)


#  Write out .sam files -------------------------------------------------------
chromosome <- "chrX"
# variable_sam <- paste0(
#     "sam.",
#     c("129.lifted", "129.unlifted", "CAST.lifted", "CAST.unlifted", "mm10")
# )
variable_sam <- paste0("sam.", c("129", "CAST", "mm10"))

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


#  Remove the (unneeded) .sam files -------------------------------------------
for (i in 1:length(variable_sam)) {
    shell_code <- paste0("rm ", file_sam[i])

    system(shell_code)
}
rm(i)


# #  Save the .sam files as .rds files ------------------------------------------
# operation <- paste0(
#     "saveRDS(",
#         variable_sam, ", ",
#         "file = ", "\"", stringr::str_replace(file_sam, ".sam", ".bam.rds"), "\"",
#     ")"
# )
# evaluateOperation(operation)


# #  Check on NAs between pairs -------------------------------------------------
# uniline.tbl.129$pos.odd %>% is.na() %>% table(useNA = "ifany")  # 111430
# uniline.tbl.129$pos.even %>% is.na() %>% table(useNA = "ifany")  # 111430
# uniline.tbl.CAST$pos.odd %>% is.na() %>% table(useNA = "ifany")  # 108755
# uniline.tbl.CAST$pos.even %>% is.na() %>% table(useNA = "ifany")  # 108755
# uniline.tbl.mm10$pos.odd %>% is.na() %>% table(useNA = "ifany")  # 117342
# uniline.tbl.mm10$pos.even %>% is.na() %>% table(useNA = "ifany")  # 117342


#  Clean up environment, save environment image -------------------------------
rm(
    chromosome,
    dedup.129, dedup.CAST, dedup.mm10,
    file, file_sam,
    sam.129.lifted, sam.129.unlifted, sam.CAST.lifted, sam.CAST.unlifted,
    sam.129, sam.CAST, sam.mm10,
    shell_code,
    tbl.129, tbl.CAST, tbl.mm10,
    variable, variable_bam_initial, variable_lO, variable_sam, variable_uniline
)

save.image(stringr::str_remove(script, ".R") %>% paste0(., ".Rdata"))


#  Script is completed; clean up environment ----------------------------------
rm(list = ls())

#  Go to part 5
