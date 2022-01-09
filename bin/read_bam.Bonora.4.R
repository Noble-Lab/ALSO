#!/usr/bin/env Rscript

#  Load libraries
library(dplyr)
library(forcats)
library(ggplot2)
library(magrittr)
library(parallel)
library(purrr)
library(Rsamtools)
library(scales)

options(pillar.sigfig = 8, scipen = 10000)


# -----------------------------------------------------------------------------
#  Temporary: Set up work directory (location TB∆)
"/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora/segregatedReads.SNPTHRESH1.Q30" %>% setwd()

#  Files are from
#+ /net/noble/vol2/home/gbonora/proj/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouseDiff_Nmasked/segregatedReads.SNPTHRESH1.Q30
#+ 
#+ For more information, see the following script:
#+ /net/noble/vol2/home/gbonora/proj/2019_sciATAC_analysis/results/gbonora/20191105_sciATAC_mouseDiff_Nmasked/20191105_sciATAC_mouseDiff_Nmasked_allelicSegregation_workflow.sh

# -----------------------------------------------------------------------------
`%notin%` <- Negate(`%in%`)

makeOperation <- function(variable, command) {
    operation <- paste(variable, command) %>% paste(., collapse = "; ")
    return(operation)
}

#  What to query from .bam files
map_info <- c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar", "mrnm", "mpos", "isize", "seq", "qual")
tag_info <- c("AS", "XS", "NM")
map_params <- Rsamtools::ScanBamParam(what = map_info, tag = tag_info)

#  Assign variables necessary for subsequent commands
files <- list.files(pattern = "\\.chr1.bam$")
n <- length(files)
variables <- paste0("bam", 1:n)
variables_rname <- paste0("bam", 1:n, "$rname")
command_pipe <- paste0("<- bam", 1:n, " %>% ")
chromosomes <- c(paste0("chr", 1:19), "chrX")

mapply(
    assign, variables, files, MoreArgs = list(envir = parent.frame())
)

#  Assign .bam information as list
GB.alt.CAST <- bam1 %>% Rsamtools::scanBam(., param = map_params)
GB.ref.129 <- bam4 %>% Rsamtools::scanBam(., param = map_params)
GB.ambig <- bam2 %>% Rsamtools::scanBam(., param = map_params)
GB.contra <- bam3 %>% Rsamtools::scanBam(., param = map_params)

#  Convert .bam information from list to dataframe to tibble
variables <- c(
    "GB.alt.CAST",
    "GB.ref.129",
    "GB.ambig",
    "GB.contra"
)
command_pipe <- paste0("<- ", variables, " %>% ")

command <- paste0(command_pipe, "as.data.frame() %>% as_tibble()")
operation <- makeOperation(variables, command)
eval(parse(text = operation))
rm("bam1", "bam2", "bam3", "bam4")

#  Reorder rname factor levels
variables_rname <- paste0(variables, "$rname")
command <- paste0("<- forcats::fct_relevel(", variables_rname, ", chromosomes)")
operation <- makeOperation(variables_rname, command)
eval(parse(text = operation))

#  Drop unused rname factor levels
command <- paste0("<- ", variables_rname, " %>% forcats::fct_drop()")
operation <- makeOperation(variables_rname, command)
eval(parse(text = operation))

#  Drop rows that are not chr1-19, chrX
command <- paste0(command_pipe, "filter(., rname %in% chromosomes)")
operation <- makeOperation(variables, command)
eval(parse(text = operation))

GB.alt.CAST %>% head()
GB.ref.129 %>% head()
GB.ambig %>% head()
GB.contra %>% head()


# -----------------------------------------------------------------------------
#  Create a column of concatenated barcode_fragment strings; arrange data
#+ frames by rname, pos, and AS; then, filter out all rows without non-unique
#+ combinations of qname (barcode) and seq (read sequence) values
b_f.GB.alt.CAST <- GB.alt.CAST %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    ) %>%
    dplyr::arrange(rname, pos)

b_f.GB.ref.129 <- GB.ref.129 %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    ) %>%
    dplyr::arrange(rname, pos)

b_f.GB.ambig <- GB.ambig %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    ) %>%
    dplyr::arrange(rname, pos)

b_f.GB.contra <- GB.contra %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    ) %>%
    dplyr::arrange(rname, pos)


# b_f.GB.alt.CAST.distinct <- GB.alt.CAST %>%
#     tidyr::unite(
#         b_f,
#         c("qname", "seq"),
#         sep = "_",
#         remove = FALSE
#     ) %>%
#     dplyr::arrange(rname, pos) %>%
#     dplyr::distinct(qname, seq, .keep_all = TRUE)
# 
# b_f.GB.ref.129.distinct <- GB.ref.129 %>%
#     tidyr::unite(
#         b_f,
#         c("qname", "seq"),
#         sep = "_",
#         remove = FALSE
#     ) %>%
#     dplyr::arrange(rname, pos) %>%
#     dplyr::distinct(qname, seq, .keep_all = TRUE)
# 
# b_f.GB.ambig.distinct <- GB.ambig %>%
#     tidyr::unite(
#         b_f,
#         c("qname", "seq"),
#         sep = "_",
#         remove = FALSE
#     ) %>%
#     dplyr::arrange(rname, pos) %>%
#     dplyr::distinct(qname, seq, .keep_all = TRUE)
# 
# b_f.GB.contra.distinct <- GB.contra %>%
#     tidyr::unite(
#         b_f,
#         c("qname", "seq"),
#         sep = "_",
#         remove = FALSE
#     ) %>%
#     dplyr::arrange(rname, pos) %>%
#     dplyr::distinct(qname, seq, .keep_all = TRUE)


# -----------------------------------------------------------------------------
# source("/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/bin/read_bam.Bonora.3.R")
#  For now, do the above manually

# -----------------
assignment.CAST[assignment.CAST$b_f %in% assignment.CAST$b_f, ]  # 122,033
assignment.CAST[assignment.CAST$b_f %notin% assignment.CAST$b_f, ]  # 0

assignment.CAST[assignment.CAST$b_f %in% b_f.GB.alt.CAST$b_f, ]  # 68,582
assignment.CAST[assignment.CAST$b_f %notin% b_f.GB.alt.CAST$b_f, ]  # 53,451

assignment.CAST[assignment.CAST$b_f %in% b_f.GB.ref.129$b_f, ]  # 3,245
assignment.CAST[assignment.CAST$b_f %notin% b_f.GB.ref.129$b_f, ]  # 118,788

assignment.CAST[assignment.CAST$b_f %in% b_f.GB.ambig$b_f, ]  # 16,196
assignment.CAST[assignment.CAST$b_f %notin% b_f.GB.ambig$b_f, ]  # 105,837

assignment.CAST[assignment.CAST$b_f %in% b_f.GB.contra$b_f, ]  # 117
assignment.CAST[assignment.CAST$b_f %notin% b_f.GB.contra$b_f, ]  # 121,916

# -----------------
assignment.129[assignment.129$b_f %in% assignment.129$b_f, ]  # 121,076
assignment.129[assignment.129$b_f %notin% assignment.129$b_f, ]  # 0

assignment.129[assignment.129$b_f %in% b_f.GB.ref.129$b_f, ]  # 68,712
assignment.129[assignment.129$b_f %notin% b_f.GB.ref.129$b_f, ]  # 52,364

assignment.129[assignment.129$b_f %in% b_f.GB.alt.CAST$b_f, ]  # 4,385
assignment.129[assignment.129$b_f %notin% b_f.GB.alt.CAST$b_f, ]  # 116,691

assignment.129[assignment.129$b_f %in% b_f.GB.ambig$b_f, ]  # 21,720
assignment.129[assignment.129$b_f %notin% b_f.GB.ambig$b_f, ]  # 99,356

assignment.129[assignment.129$b_f %in% b_f.GB.contra$b_f, ]  # 98
assignment.129[assignment.129$b_f %notin% b_f.GB.contra$b_f, ]  # 120,978

# -----------------
#  Logical "and"
assignment.CAST[assignment.CAST$b_f %in% b_f.GB.alt.CAST$b_f &
                assignment.CAST$b_f %in% b_f.GB.ref.129$b_f, ]  # 180
assignment.CAST[assignment.CAST$b_f %in% b_f.GB.alt.CAST$b_f &
                assignment.CAST$b_f %in% b_f.GB.ambig$b_f, ]  # 1,290

assignment.129[assignment.129$b_f %in% b_f.GB.ref.129$b_f &
               assignment.129$b_f %in% b_f.GB.alt.CAST$b_f, ]  # 190
assignment.129[assignment.129$b_f %in% b_f.GB.ref.129$b_f &
               assignment.129$b_f %in% b_f.GB.ambig$b_f, ]  # 1,508

assignment.CAST[assignment.CAST$b_f %notin% b_f.GB.alt.CAST$b_f &
                    assignment.CAST$b_f %notin% b_f.GB.ref.129$b_f, ]  # 50,386
assignment.CAST[assignment.CAST$b_f %notin% b_f.GB.alt.CAST$b_f &
                    assignment.CAST$b_f %notin% b_f.GB.ambig$b_f, ]  # 38,545

assignment.129[assignment.129$b_f %notin% b_f.GB.ref.129$b_f &
                   assignment.129$b_f %notin% b_f.GB.alt.CAST$b_f, ]  # 48,169
assignment.129[assignment.129$b_f %notin% b_f.GB.ref.129$b_f &
                   assignment.129$b_f %notin% b_f.GB.ambig$b_f, ]  # 32,152

#  Logical "or"
assignment.CAST[assignment.CAST$b_f %in% c(b_f.GB.alt.CAST$b_f, b_f.GB.ref.129$b_f), ]  # 71,647
assignment.CAST[assignment.CAST$b_f %notin% c(b_f.GB.alt.CAST$b_f, b_f.GB.ref.129$b_f), ]  # 50,386

assignment.129[assignment.129$b_f %in% c(b_f.GB.alt.CAST$b_f, b_f.GB.ref.129$b_f), ]  # 72,907
assignment.129[assignment.129$b_f %notin% c(b_f.GB.alt.CAST$b_f, b_f.GB.ref.129$b_f), ]  # 48,169

assignment.CAST[assignment.CAST$b_f %in% b_f.GB.alt.CAST$b_f |
                    assignment.CAST$b_f %in% b_f.GB.ref.129$b_f, ]  # 71,647
assignment.CAST[assignment.CAST$b_f %notin% b_f.GB.alt.CAST$b_f |
                    assignment.CAST$b_f %notin% b_f.GB.ref.129$b_f, ]  # 121,853
assignment.CAST[!assignment.CAST$b_f %in% b_f.GB.alt.CAST$b_f |
                    !assignment.CAST$b_f %in% b_f.GB.ref.129$b_f, ]  # 121,853

# -----------------
b_f.GB.alt.CAST[b_f.GB.alt.CAST$b_f %in% b_f.GB.alt.CAST$b_f, ]  # 2,599,984
b_f.GB.alt.CAST[b_f.GB.alt.CAST$b_f %notin% b_f.GB.alt.CAST$b_f, ]  # 0

b_f.GB.alt.CAST[b_f.GB.alt.CAST$b_f %in% b_f.GB.ref.129$b_f, ]  # 27,660
b_f.GB.alt.CAST[b_f.GB.alt.CAST$b_f %notin% b_f.GB.ref.129$b_f, ]  # 2,572,324

b_f.GB.alt.CAST[b_f.GB.alt.CAST$b_f %in% b_f.GB.ambig$b_f, ]  # 265,572
b_f.GB.alt.CAST[b_f.GB.alt.CAST$b_f %notin% b_f.GB.ambig$b_f, ]  # 2,334,402

b_f.GB.alt.CAST[b_f.GB.alt.CAST$b_f %in% b_f.GB.contra$b_f, ]  # 628
b_f.GB.alt.CAST[b_f.GB.alt.CAST$b_f %notin% b_f.GB.contra$b_f, ]  # 2,599,356

# -----------------
b_f.GB.ref.129[b_f.GB.ref.129$b_f %in% b_f.GB.alt.CAST$b_f, ]  # 27,621
b_f.GB.ref.129[b_f.GB.ref.129$b_f %notin% b_f.GB.alt.CAST$b_f, ]  # 2,628,270

b_f.GB.ref.129[b_f.GB.ref.129$b_f %in% b_f.GB.ref.129$b_f, ]  # 2,655,891
b_f.GB.ref.129[b_f.GB.ref.129$b_f %notin% b_f.GB.ref.129$b_f, ]  # 0

b_f.GB.ref.129[b_f.GB.ref.129$b_f %in% b_f.GB.ambig$b_f, ]  # 271,372
b_f.GB.ref.129[b_f.GB.ref.129$b_f %notin% b_f.GB.ambig$b_f, ]  # 2,384,519

b_f.GB.ref.129[b_f.GB.ref.129$b_f %in% b_f.GB.contra$b_f, ]  # 668
b_f.GB.ref.129[b_f.GB.ref.129$b_f %notin% b_f.GB.contra$b_f, ]  # 2,655,223

#  However, GB assigned the reads, qname_seq combinations appear to have been
#+ assigned to more than one type...


# -----------------------------------------------------------------------------
x <- c(10,11,12)
y <- c(6,7,8,9)

# x <- c(1,5,10,20,20,24,45)
# y <- c(1,5,10,20)

z <- c(5,10)
# z <- c(5,45)

z %in% x
z %in% y
z %in% c(x, y)  # Logical "or," i.e., "z in x" OR "z in y"
z %in% x | z %in% y  # Logical "or," i.e., "z in x" OR "z in y"
z %in% x & z %in% y  # Logical "and," i.e., "z in x" AND "z in y"

z %notin% x
z %notin% y
z %notin% c(x, y)  # Negation of logical "or," i.e., "z not in x" AND "z not in y"
z %notin% x | z %notin% y  # Logical "or," i.e., "z not in x" AND "z not in y"
z %notin% x & z %notin% y  # Logical "and," i.e., "z not in x" OR "z not in y"

rm(x, y, z)
