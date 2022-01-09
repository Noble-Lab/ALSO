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
    #  Returns a vector of rows for vector entries ± one entry
    v1 <- vector
    v2 <- v1 - 1
    v3 <- v1 + 1
    vector_plus_minus <- c(v1, v2, v3) %>% sort() %>% unique()
    return(vector_plus_minus)
}


createVectorDiscrepancy <- function(discrepancy) {
    #  Returns a vector of rows for discrepancy entries ± one entry
    v1 <- which(discrepancy == TRUE)
    v2 <- v1 - 1
    v3 <- v1 + 1
    discrepancy_plus_minus <- c(v1, v2, v3) %>% sort() %>% unique()
    return(discrepancy_plus_minus)
}


makeOperation <- function(variable, command) {
    operation <- paste(variable, command) %>% paste(., collapse = "; ")
    return(operation)
}


testLogicalAlternating <- function(logical) {
    output <- vector(mode = "logical", length = length(output.1))
    for (i in 1:length(output.1)) {
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


#  Set up what to query from .bam files ---------------------------------------
# map_info <- c(
#     "qname", "flag", "rname", "pos", "mapq", "cigar", "mrnm", "mpos", "isize",
#     "seq", "qual"
# )
# flag_info <- Rsamtools::scanBamFlag(isProperPair = TRUE)
tag_info <- "AS"
map_params <- Rsamtools::ScanBamParam(
    # what = map_info,
    # flag = flag_info,
    tag = tag_info
)


#  Assign variables necessary for subsequent commands -------------------------
# chromosome <- "chr1"
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
command <- paste0(
    "<- ", variable,
    " %>% Rsamtools::BamFile(., index = ", index, ", asMates = TRUE)",
    " %>% Rsamtools::scanBam()",
    " %>% as.data.frame()",
    " %>% tibble::as_tibble()"
)
operation <- makeOperation(paste0(variable, ".full"), command)
eval(parse(text = operation))

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


command <- paste0(
    "<- dplyr::bind_cols(", variable, ".full, ", variable, ")"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))


command <- paste0("rm(", variable, ".full)")
eval(parse(text = command))


#  Tests ----------------------------------------------------------------------
# dedup.129S1$flag %>% table()
# dedup.129S1$mate_status %>% table()

o.ambiguous <- dedup.129S1[dedup.129S1$mate_status == "ambiguous", ]
o.mated <- dedup.129S1[dedup.129S1$mate_status == "mated", ]
o.unmated <- dedup.129S1[dedup.129S1$mate_status == "unmated", ]

# o.ambiguous$pos %>% duplicated() %>% table()
# # .
# #  FALSE   TRUE
# # 111462 148792
# 
# o.mated$pos %>% duplicated() %>% table()
# # .
# #  FALSE   TRUE
# # 126542   3426
# 
# o.unmated$pos %>% duplicated() %>% table()
# # .
# # FALSE  TRUE
# #  1499  1537


#  Reorder rname factor levels ------------------------------------------------
variable <- paste0("o.", c("ambiguous", "mated", "unmated"))
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
    command_pipe, "dplyr::relocate(pos_end, .after = pos) %>% dplyr::relocate(mpos_end, .after = mpos)"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))


#  Filter out rows with mapq less than 30 -------------------------------------
command <- paste0("<- ", variable, "[", variable, "$mapq >= 30, ]")
operation <- makeOperation(variable, command)
eval(parse(text = operation))

#  Test: Deduplicating prior to sorting, etc. ---------------------------------
# o.ambiguous$pos_mpos <- paste0(o.ambiguous$pos, "_", o.ambiguous$mpos)
o.ambiguous$criteria <- paste0(
    o.ambiguous$qname, "_",
    o.ambiguous$flag, "_",
    o.ambiguous$pos, "_",
    o.ambiguous$mpos
)

# o.ambiguous$criteria %>% head(n = 100)
# o.ambiguous$criteria %>% duplicated() %>% head(n = 100)

# Create a back-up for the test variable
o.ambiguous <- o.ambiguous.pre
o.ambiguous.pre <- o.ambiguous

# Get mpos close to pos
o.ambiguous <- o.ambiguous %>% dplyr::relocate(mpos, .after = pos)

#MAYBE
o.ambiguous$groupid <- o.ambiguous$groupid %>% as_factor()
# length_groupid <- o.ambiguous %>%
#     dplyr::group_by(groupid) %>%
#     summarize(no_rows = length(groupid))

o.ambiguous$qpos <- paste0(o.ambiguous$qname, "_", o.ambiguous$pos)
o.ambiguous$qmpos <- paste0(o.ambiguous$qname, "_", o.ambiguous$mpos)
o.ambiguous.same <- o.ambiguous[o.ambiguous$qpos == o.ambiguous$qmpos, ]
o.ambiguous <- o.ambiguous[!(o.ambiguous$qpos == o.ambiguous$qmpos), ]

test <- o.ambiguous[!duplicated(o.ambiguous$criteria), ]
# nrow(o.ambiguous) - nrow(test)  # [1] 140199

# -----------------
pos <- paste0(test$qname, "_", test$pos)
mpos <- paste0(test$qname, "_", test$mpos)

pos.in.mpos <- pos %in% mpos
mpos.in.pos <- mpos %in% pos

pos.in.mpos %>% table()
mpos.in.pos %>% table()

# -------
test.2 <- test[pos.in.mpos, ]

pos <- paste0(test.2$qname, "_", test.2$pos)
mpos <- paste0(test.2$qname, "_", test.2$mpos)

pos.in.mpos <- pos %in% mpos
mpos.in.pos <- mpos %in% pos

pos.in.mpos %>% table()
mpos.in.pos %>% table()

test.2$isize_abs <- test.2$isize %>% abs()
test.2 <- test.2 %>% dplyr::relocate(isize_abs, .after = mpos)

# -------
test.3 <- test.2[order(test.2$isize_abs, test.2$qpos, test.2$qmpos), ]

test.3$strand

pos <- test.3$qpos
mpos <- test.3$qmpos

output.1 <- vector(mode = "logical", length = length(pos))
for (i in 1:length(pos)) {
    output.1[i] <- pos[i] == mpos[i + 1] & mpos[i] == pos[i + 1]
}
output.1 %>% table()  # 54070 - 54067 = 3
# .
# FALSE  TRUE 
# 54070 54067

output.2 <- vector(mode = "logical", length = length(output.1))
for (i in 1:length(output.1)) {
    output.2[i] <- (output.1[i] == output.1[i + 1])
}
output.2 %>% table()
# .
# FALSE   TRUE
# 108132     4

test.3$discrepancy <- output.2
v1 <- which(test.3$discrepancy == TRUE)
v2 <- v1 - 1
v3 <- v1 + 1
discrepancy_plus_minus <- c(v1, v2, v3) %>% sort() %>% unique()
rm(v1, v2, v3)


# -----------------
test.4 <- test.3[discrepancy_plus_minus, ]
test.4 <- test.4[order(test.4$isize_abs, test.4$qpos, test.4$qmpos, test.4$strand), ]


# -----------------
# test.5 <- test.4[!duplicated(test.4$pos) & !duplicated(test.4$mpos), ]
test.5 <- test.3[!duplicated(test.3$pos) & !duplicated(test.3$mpos), ]

pos <- test.5$qpos
mpos <- test.5$qmpos

output.1 <- vector(mode = "logical", length = length(pos))
for (i in 1:length(pos)) {
    output.1[i] <- pos[i] == mpos[i + 1] & mpos[i] == pos[i + 1]
}
output.1 %>% table()  # 54081 - 54056 = 25
# .
# FALSE TRUE 
#     1    2

output.2 <- vector(mode = "logical", length = length(output.1))
for (i in 1:length(output.1)) {
    output.2[i] <- (output.1[i] == output.1[i + 1])
}
output.2 %>% table()
# .
# FALSE 
#     2
rm(pos, mpos, output.1, output.2, i)


