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

#  Important resources
#+ - stats.stackexchange.com/questions/10838/produce-a-list-of-variable-name-in-a-for-loop-then-assign-values-to-them
#+ - stackoverflow.com/questions/47962255/change-factor-levels-to-custom-order-of-a-column
#+ - stackoverflow.com/questions/51272510/how-to-count-unique-rows-in-a-data-frame

options(pillar.sigfig = 8, scipen = 10000)


# -----------------------------------------------------------------------------
#  Temporary: Set up work directory (location TB∆)
"/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Disteche" %>% setwd()
# getwd()


# -----------------------------------------------------------------------------
makeOperation <- function(variable, command) {
    operation <- paste(variable, command) %>% paste(., collapse = "; ")
    return(operation)
}

#  What to query from .bam files
map_info <- c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar", "mrnm", "mpos", "isize", "seq", "qual")
tag_info <- c("AS", "XS", "NM")
map_params <- Rsamtools::ScanBamParam(what = map_info, tag = tag_info)

#  Assign variables necessary for subsequent commands
# files <- list.files(pattern = "\\.dedup.MarkDuplicates.bam$")
files <- list.files(pattern = "\\.sort.chr1.bam$")
n <- files %>% length()
variables <- paste0("bam", 1:n)
variables_rname <- paste0("bam", 1:n, "$rname")
command_pipe <- paste0("<- bam", 1:n, " %>% ")
chromosomes <- c(paste0("chr", 1:19), "chrX")

#  Read in .bam filenames
# command <- paste0("<- files[", 1:n, "]")
# operation <- makeOperation(variables, command)
# eval(parse(text = operation))

mapply(
    assign, variables, files, MoreArgs = list(envir = parent.frame())
)
# rm("bam1", "bam2", "bam3", "bam4")

# cores <- parallel::detectCores() - 1 %>% as.integer()
# parallel::mcmapply(
#     assign, variables, files, MoreArgs = list(envir = parent.frame())
# )

#  Assign .bam information as list
dedup.CAST <- bam1 %>% Rsamtools::scanBam(., param = map_params)
dedup.mm10_CAST_Nmasked <- bam2 %>% Rsamtools::scanBam(., param = map_params)
dedup.mm10 <- bam3 %>% Rsamtools::scanBam(., param = map_params)

#  (A previous call; have deleted the predup files to save space...)
# dedup.CAST <- bam1 %>% Rsamtools::scanBam(., param = map_params)
# predup.CAST <- bam2 %>% Rsamtools::scanBam(., param = map_params)
# dedup.mm10_CAST_Nmasked <- bam3 %>% Rsamtools::scanBam(., param = map_params)
# predup.mm10_CAST_Nmasked <- bam4 %>% Rsamtools::scanBam(., param = map_params)
# dedup.mm10 <- bam5 %>% Rsamtools::scanBam(., param = map_params)
# predup.mm10 <- bam6 %>% Rsamtools::scanBam(., param = map_params)

# command <- paste0(command_pipe, "Rsamtools::scanBam(., param = map_params)")
# operation <- makeOperation(variables, command)
# eval(parse(text = operation))
# parallel::mcparallel(eval(parse(text = operation)))  #TODO

#TODO Get this working
# mapply(
#     assign, variables, eval(parse(text = operation)), MoreArgs = list(envir = parent.frame())
# )

#TODO Get this working
# cores <- parallel::detectCores() - 1 %>% as.integer()
# parallel::mcmapply(
#     assign,
#     variables,
#     command,
#     MoreArgs = list(envir = parent.frame()),
#     mc.cores = 4
# )

#  Convert .bam information from list to dataframe to tibble
variables <- c(
    "dedup.CAST",
    # "predup.CAST",
    "dedup.mm10_CAST_Nmasked",
    # "predup.mm10_CAST_Nmasked",
    "dedup.mm10"#,
    # "predup.mm10"
)
command_pipe <- paste0("<- ", variables, " %>% ")

command <- paste0(command_pipe, "as.data.frame() %>% as_tibble()")
operation <- makeOperation(variables, command)
eval(parse(text = operation))

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

# #  Sort rows by qname, rname, pos, mpos, and qual (very slow)
# command <- paste0(command_pipe, "arrange(., qname, rname, pos, mpos, qual)")
# operation <- makeOperation(variables, command)
# eval(parse(text = operation))

dedup.CAST %>% head()
dedup.mm10 %>% head()
dedup.mm10_CAST_Nmasked %>% head()

#  (A previous call; have deleted the predup files to save space...)
# predup.CAST %>% head()
# dedup.CAST %>% head()
# predup.mm10 %>% head()
# dedup.mm10 %>% head()
# predup.mm10_CAST_Nmasked %>% head()
# dedup.mm10_CAST_Nmasked %>% head()

# #  How many rows are there with unique rname, pos, mpos combinations (very, very slow)?
# variables_count_unique <- paste0("bam", 1:n, "_count_unique")
# command <- paste0("<- subset(", variables, ", select = c(qname, rname, pos, mpos, qual)) %>% unique() %>% nrow()")
# operation <- makeOperation(variables_count_unique, command)
# eval(parse(text = operation))
# 
# variables_difference <- paste0("bam", 1:n, "_difference")
# command <- paste0("<- nrow(", variables, ") - ", variables_count_unique)
# operation <- makeOperation(variables_difference, command)
# eval(parse(text = operation))
# 
# #  What are the duplicate rows (very, very slow)?
# variables_what <- paste0("bam", 1:n, "_what_is_duplicated")
# command <- paste0(command_pipe, "dplyr::group_by(qname, rname, pos, mpos, qual) %>% dplyr::count() %>% arrange(desc(n))")
# operation <- makeOperation(variables_what, command)
# eval(parse(text = operation))


# #  Numbers of unique qname (barcodes)
# predup.CAST.qname.unique <- predup.CAST$qname %>% unique() %>% length()  # [1] 4083661
# length(predup.CAST$qname) - predup.CAST.qname.unique  # [1] 4083697
# 
# dedup.CAST.qname.unique <- dedup.CAST$qname %>% unique() %>% length()  # [1] 2869453
# length(dedup.CAST$qname) - dedup.CAST.qname.unique  # [1] 2869481
# 
# predup.mm10.qname.unique <- predup.mm10$qname %>% unique() %>% length()  # [1] 4311150
# length(predup.mm10$qname) - predup.mm10.qname.unique  # [1] 4311184
# 
# dedup.mm10.qname.unique <- dedup.mm10$qname %>% unique() %>% length()  # [1] 3033308
# length(dedup.mm10$qname) - dedup.mm10.qname.unique  # [1] 3033332
#
# #  Numbers of unique seq (fragments)
# predup.CAST.seq.unique <- predup.CAST$seq %>% unique() %>% length()  # [1] 5516920
# length(predup.CAST$seq) - predup.CAST.seq.unique  # [1] 2650438
# 
# dedup.CAST.seq.unique <- dedup.CAST$seq %>% unique() %>% length()  # [1] 5340687
# length(dedup.CAST$seq) - dedup.CAST.seq.unique  # [1] 398247
# 
# predup.mm10.seq.unique <- predup.mm10$seq %>% unique() %>% length()  # [1] 5845476
# length(predup.mm10$seq) - predup.mm10.seq.unique  # [1] 2776858
# 
# dedup.mm10.seq.unique <- dedup.mm10$seq %>% unique() %>% length()  # [1] 5654808
# length(dedup.mm10$seq) - dedup.mm10.seq.unique  # [1] 411832


#  Check distributions of .bam flags ------------------------------------------
# predup.CAST$flag %>% table()
# .
# 83      99     147     163    1107    1123    1171    1187 
# 1416048 1416858 1416858 1416048  624616  626157  626157  624616

dedup.CAST$flag %>% table()
# .
# 83      99     147     163    1107    1123    1171    1187 
# 1416029 1416877 1416877 1416029   18318   18243   18243   18318

# predup.mm10_CAST_Nmasked$flag %>% table()
# .
# 83      99     147     163    1107    1123    1171    1187 
# 1491317 1492138 1492138 1491317  652071  653668  653668  652071

dedup.mm10_CAST_Nmasked$flag %>% table()
# .
# 83      99     147     163    1107    1123    1171    1187 
# 1491322 1492133 1492133 1491322   17624   17587   17587   17624

# predup.mm10$flag %>% table()
# .
# 83      99     147     163    1107    1123    1171    1187 
# 1498245 1499637 1499637 1498245  655776  657509  657509  655776

dedup.mm10$flag %>% table()
# .
# 83      99     147     163    1107    1123    1171    1187 
# 1498254 1499628 1499628 1498254   17741   17697   17697   17741


#  Remove "duplicate" status rows ---------------------------------------------
`%notin%` <- Negate(`%in%`)
duplicate_flags <- c(1107, 1123, 1171, 1187)
# filter.predup.CAST <- predup.CAST %>% dplyr::filter(flag %notin% duplicate_flags)
# filter.predup.CAST$flag %>% table()
# .
# 83      99     147     163 
# 1416048 1416858 1416858 1416048

filter.dedup.CAST <- dedup.CAST %>% dplyr::filter(flag %notin% duplicate_flags)
filter.dedup.CAST$flag %>% table()
# .
# 83      99     147     163 
# 1416029 1416877 1416877 1416029

# length(filter.predup.CAST$qname) - length(filter.dedup.CAST$qname)  # [1] 0


# filter.predup.mm10_CAST_Nmasked <- predup.mm10_CAST_Nmasked %>% dplyr::filter(flag %notin% duplicate_flags)
# filter.predup.mm10_CAST_Nmasked$flag %>% table()
# .
# 83      99     147     163 
# 1491317 1492138 1492138 1491317

filter.dedup.mm10_CAST_Nmasked <- dedup.mm10_CAST_Nmasked %>% dplyr::filter(flag %notin% duplicate_flags)
filter.dedup.mm10_CAST_Nmasked$flag %>% table()
# .
# 83      99     147     163 
# 1491322 1492133 1492133 1491322

# length(filter.predup.mm10_CAST_Nmasked$qname) - length(filter.dedup.mm10_CAST_Nmasked$qname)  # [1] 0


# filter.predup.mm10 <- predup.mm10 %>% dplyr::filter(flag %notin% duplicate_flags)
# filter.predup.mm10$flag %>% table()
# .
# 83      99     147     163 
# 1498245 1499637 1499637 1498245

filter.dedup.mm10 <- dedup.mm10 %>% dplyr::filter(flag %notin% duplicate_flags)
filter.dedup.mm10$flag %>% table()
# .
# 83      99     147     163 
# 1498254 1499628 1499628 1498254

# length(filter.predup.mm10$qname) - length(filter.dedup.mm10$qname)  # [1] 0


# #  Numbers of unique qname (barcodes)
# filter.predup.CAST.qname.unique <- filter.predup.CAST$qname %>% unique() %>% length()  # [1] 2832899
# length(filter.predup.CAST$qname) - filter.predup.CAST.qname.unique  # [1] 2832913
# 
# filter.dedup.CAST.qname.unique <- filter.dedup.CAST$qname %>% unique() %>% length()  # [1] 2832895
# length(filter.dedup.CAST$qname) - filter.dedup.CAST.qname.unique  # [1] 2832917
# 
# filter.predup.mm10.qname.unique <- filter.predup.mm10$qname %>% unique() %>% length()  # [1] 2997876
# length(filter.predup.mm10$qname) - filter.predup.mm10.qname.unique  # [1] 2997888
# 
# filter.dedup.mm10.qname.unique <- filter.dedup.mm10$qname %>% unique() %>% length()  # [1] 2997873
# length(filter.dedup.mm10$qname) - filter.dedup.mm10.qname.unique  # [1] 2997891
# 
# #  Numbers of unique seq (fragments)
# filter.predup.CAST.seq.unique <- filter.predup.CAST$seq %>% unique() %>% length()  # [1] 5327789
# length(filter.predup.CAST$seq) - filter.predup.CAST.seq.unique  # [1] 338023
# 
# filter.dedup.CAST.seq.unique <- filter.dedup.CAST$seq %>% unique() %>% length()  # [1] 5333351
# length(filter.dedup.CAST$seq) - filter.dedup.CAST.seq.unique  # [1] 332461
# 
# filter.predup.mm10.seq.unique <- filter.predup.mm10$seq %>% unique() %>% length()  # [1] 5642015
# length(filter.predup.mm10$seq) - filter.predup.mm10.seq.unique  # [1] 353749
# 
# filter.dedup.mm10.seq.unique <- filter.dedup.mm10$seq %>% unique() %>% length()  # [1] 5648188
# length(filter.dedup.mm10$seq) - filter.dedup.mm10.seq.unique  # [1] 347576


# -------------------------------------
# dedup.CAST.seq ----------------------
# -----------------
# dedup.CAST.seq.n_occur <- data.frame(table(dedup.CAST$seq))  #TODO Faster way to do this
# dedup.CAST.seq.n_occur %>% head()
# dedup.CAST.seq.n_occur %>% nrow()  # [1] 5340687
# 
# nrow(dedup.CAST) - nrow(dedup.CAST.seq.n_occur)  # [1] 398247
# dedup.CAST$seq %>% unique() %>% length()  # [1] 5340687
# 
# # -----------------
# dedup.CAST.seq.n_occur_gt1 <- dedup.CAST.seq.n_occur[dedup.CAST.seq.n_occur$Freq > 1, ]
# dedup.CAST.seq.n_occur_gt1 %>% head()
# dedup.CAST.seq.n_occur_gt1 %>% nrow()  # [1] 323968
# dedup.CAST.seq.n_occur_gt1$Freq %>% sum()  # [1] 722215
# 
# sum(dedup.CAST.seq.n_occur_gt1$Freq) - nrow(dedup.CAST.seq.n_occur_gt1)  # [1] 398247
# 
# # -----------------
# dedup.CAST.seq.n_seq_gt1 <- dedup.CAST[dedup.CAST$seq %in% dedup.CAST.seq.n_occur_gt1$Var1, ]
# dedup.CAST.seq.n_seq_gt1 %>% head()
# dedup.CAST.seq.n_seq_gt1 %>% nrow()  # [1] 722215
# 
# dedup.CAST.seq.n_seq_gt1$seq %>% unique() %>% length()  # [1] 323968
# length(dedup.CAST.seq.n_seq_gt1$seq) - length(unique(dedup.CAST.seq.n_seq_gt1$seq))  # [1] 398247

# # -----------------
b_f.dedup.CAST <- dedup.CAST %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    )
# b_f.dedup.CAST$b_f %>% unique() %>% length()  # [1] 5738932
# nrow(dedup.CAST)  # [1] 5738934
# 
# nrow(dedup.CAST) - length(unique(b_f.dedup.CAST$b_f))  # [1] 2
# 
# dedup.CAST.b_f.n_occur <- b_f.dedup.CAST$b_f %>%
#     table() %>%
#     data.frame()
# dedup.CAST.b_f.n_gt1 <- dedup.CAST.b_f.n_occur[dedup.CAST.b_f.n_occur$Freq > 1, ]
# dedup.CAST.b_f.n_seq_gt1 <- b_f.dedup.CAST[b_f.dedup.CAST$b_f %in% dedup.CAST.b_f.n_gt1$., ]


# -------------------------------------
# dedup.CAST.qname --------------------
# -----------------
# dedup.CAST.qname.n_occur <- data.frame(table(dedup.CAST$qname))  #TODO Faster way to do this
# dedup.CAST.qname.n_occur %>% head()
# dedup.CAST.qname.n_occur %>% nrow()  # [1] 2869453
# 
# nrow(dedup.CAST) - nrow(dedup.CAST.qname.n_occur)  # [1] 2869481
# dedup.CAST$qname %>% unique() %>% length()  # [1] 2869453
# 
# # -----------------
# dedup.CAST.qname.n_occur_gt1 <- dedup.CAST.qname.n_occur[dedup.CAST.qname.n_occur$Freq > 1, ]
# dedup.CAST.qname.n_occur_gt1 %>% head()
# dedup.CAST.qname.n_occur_gt1 %>% nrow()  # [1] 2869453
# dedup.CAST.qname.n_occur_gt1$Freq %>% sum()  # [1] 5738934
# 
# sum(dedup.CAST.qname.n_occur_gt1$Freq) - nrow(dedup.CAST.qname.n_occur_gt1)  # [1] 2869481
# 
# # -----------------
# dedup.CAST.qname.n_qname_gt1 <- dedup.CAST[dedup.CAST$qname %in% dedup.CAST.qname.n_occur_gt1$Var1, ]
# dedup.CAST.qname.n_qname_gt1 %>% head()
# dedup.CAST.qname.n_qname_gt1 %>% nrow()  # [1] 5738934
# 
# dedup.CAST.qname.n_qname_gt1$qname %>% unique() %>% length()  # [1] 2869453
# length(dedup.CAST.qname.n_qname_gt1$qname) - length(unique(dedup.CAST.qname.n_qname_gt1$qname))  # [1] 2869481

# -----------------
# b_f.dedup.CAST <- dedup.CAST %>%
#     tidyr::unite(
#         b_f,
#         c("qname", "seq"),
#         sep = "_",
#         remove = FALSE
#     )
# b_f.dedup.CAST$b_f %>% unique() %>% length()  # [1] 5738932
# nrow(dedup.CAST)  # [1] 5738934
# 
# nrow(dedup.CAST) - length(unique(b_f.dedup.CAST$b_f))  # [1] 2
# 
# dedup.CAST.b_f.n_occur <- b_f.dedup.CAST$b_f %>%
#     table() %>%
#     data.frame()
# dedup.CAST.b_f.n_gt1 <- dedup.CAST.b_f.n_occur[dedup.CAST.b_f.n_occur$Freq > 1, ]
# dedup.CAST.b_f.n_qname_gt1 <- b_f.dedup.CAST[b_f.dedup.CAST$b_f %in% dedup.CAST.b_f.n_gt1$., ]

# -------------------------------------
# filter.dedup.CAST.seq ---------------------
# -----------------
filter.dedup.CAST.seq.n_occur <- data.frame(table(filter.dedup.CAST$seq))
filter.dedup.CAST.seq.n_occur %>% head()
filter.dedup.CAST.seq.n_occur %>% nrow()  # [1] 5333351

nrow(filter.dedup.CAST) - nrow(filter.dedup.CAST.seq.n_occur)  # [1] 332461
filter.dedup.CAST$seq %>% unique() %>% length()  # [1] 5333351

# -----------------
filter.dedup.CAST.seq.n_occur_gt1 <- filter.dedup.CAST.seq.n_occur[filter.dedup.CAST.seq.n_occur$Freq > 1, ]
filter.dedup.CAST.seq.n_occur_gt1 %>% head()
filter.dedup.CAST.seq.n_occur_gt1 %>% nrow()  # [1] 271487
filter.dedup.CAST.seq.n_occur_gt1$Freq %>% sum()  # [1] 603948

sum(filter.dedup.CAST.seq.n_occur_gt1$Freq) - nrow(filter.dedup.CAST.seq.n_occur_gt1)  # [1] 332461

# -----------------
filter.dedup.CAST.seq.n_seq_gt1 <- filter.dedup.CAST[filter.dedup.CAST$seq %in% filter.dedup.CAST.seq.n_occur_gt1$Var1, ]
filter.dedup.CAST.seq.n_seq_gt1 %>% head()
filter.dedup.CAST.seq.n_seq_gt1 %>% nrow()  # [1] 603948

filter.dedup.CAST.seq.n_seq_gt1$seq %>% unique() %>% length()  # [1] 271487
length(filter.dedup.CAST.seq.n_seq_gt1$seq) - length(unique(filter.dedup.CAST.seq.n_seq_gt1$seq))  # [1] 332461

# -----------------
b_f.filter.dedup.CAST <- filter.dedup.CAST %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    )
b_f.filter.dedup.CAST$b_f %>% unique() %>% length()  # [1] 5665810
nrow(filter.dedup.CAST)  # [1] 5665812

nrow(filter.dedup.CAST) - length(unique(b_f.filter.dedup.CAST$b_f))  # [1] 2

filter.dedup.CAST.b_f.n_occur <- b_f.filter.dedup.CAST$b_f %>%
    table() %>%
    data.frame()
filter.dedup.CAST.b_f.n_gt1 <- filter.dedup.CAST.b_f.n_occur[filter.dedup.CAST.b_f.n_occur$Freq > 1, ]
filter.dedup.CAST.b_f.n_seq_gt1 <- b_f.filter.dedup.CAST[b_f.filter.dedup.CAST$b_f %in% filter.dedup.CAST.b_f.n_gt1$., ]

# -------------------------------------
# filter.predup.CAST.seq ---------------------
# -----------------
# filter.predup.CAST.seq.n_occur <- data.frame(table(filter.predup.CAST$seq))
# filter.predup.CAST.seq.n_occur %>% head()
# filter.predup.CAST.seq.n_occur %>% nrow()  # [1] 5327789
# 
# nrow(filter.predup.CAST) - nrow(filter.predup.CAST.seq.n_occur)  # [1] 338023
# filter.predup.CAST$seq %>% unique() %>% length()  # [1] 5327789
# 
# # -----------------
# filter.predup.CAST.seq.n_occur_gt1 <- filter.predup.CAST.seq.n_occur[filter.predup.CAST.seq.n_occur$Freq > 1, ]
# filter.predup.CAST.seq.n_occur_gt1 %>% head()
# filter.predup.CAST.seq.n_occur_gt1 %>% nrow()  # [1] 275969
# filter.predup.CAST.seq.n_occur_gt1$Freq %>% sum()  # [1] 613992
# 
# sum(filter.predup.CAST.seq.n_occur_gt1$Freq) - nrow(filter.predup.CAST.seq.n_occur_gt1)  # [1] 338023
# 
# # -----------------
# filter.predup.CAST.seq.n_seq_gt1 <- filter.predup.CAST[filter.predup.CAST$seq %in% filter.predup.CAST.seq.n_occur_gt1$Var1, ]
# filter.predup.CAST.seq.n_seq_gt1 %>% head()
# filter.predup.CAST.seq.n_seq_gt1 %>% nrow()  # [1] 613992
# 
# filter.predup.CAST.seq.n_seq_gt1$seq %>% unique() %>% length()  # [1] 275969
# length(filter.predup.CAST.seq.n_seq_gt1$seq) - length(unique(filter.predup.CAST.seq.n_seq_gt1$seq))  # [1] 338023

# -----------------
# b_f.filter.predup.CAST <- filter.predup.CAST %>%
#     tidyr::unite(
#         b_f,
#         c("qname", "seq"),
#         sep = "_",
#         remove = FALSE
#     )
# b_f.filter.predup.CAST$b_f %>% unique() %>% length()  # [1] 5665810
# nrow(filter.predup.CAST)  # [1] 5665812
# 
# nrow(filter.predup.CAST) - length(unique(b_f.filter.predup.CAST$b_f))  # [1] 2
# 
# filter.predup.CAST.b_f.n_occur <- b_f.filter.predup.CAST$b_f %>%
#     table() %>%
#     data.frame()
# filter.predup.CAST.b_f.n_gt1 <- filter.predup.CAST.b_f.n_occur[filter.predup.CAST.b_f.n_occur$Freq > 1, ]
# filter.predup.CAST.b_f.n_seq_gt1 <- b_f.filter.predup.CAST[b_f.filter.predup.CAST$b_f %in% filter.predup.CAST.b_f.n_gt1$., ]

# -------------------------------------
# dedup.mm10_CAST_Nmasked.seq ---------
# -----------------
# dedup.mm10_CAST_Nmasked.seq.n_occur <- data.frame(table(dedup.mm10_CAST_Nmasked$seq))
# dedup.mm10_CAST_Nmasked.seq.n_occur %>% head()
# dedup.mm10_CAST_Nmasked.seq.n_occur %>% nrow()  # [1] 5627892
# 
# nrow(dedup.mm10_CAST_Nmasked) - nrow(dedup.mm10_CAST_Nmasked.seq.n_occur)  # [1] 409440
# dedup.mm10_CAST_Nmasked$seq %>% unique() %>% length()  # [1] 5627892
# 
# # -----------------
# dedup.mm10_CAST_Nmasked.seq.n_occur_gt1 <- dedup.mm10_CAST_Nmasked.seq.n_occur[dedup.mm10_CAST_Nmasked.seq.n_occur$Freq > 1, ]
# dedup.mm10_CAST_Nmasked.seq.n_occur_gt1 %>% head()
# dedup.mm10_CAST_Nmasked.seq.n_occur_gt1 %>% nrow()  # [1] 338977
# dedup.mm10_CAST_Nmasked.seq.n_occur_gt1$Freq %>% sum()  # [1] 748417
# 
# sum(dedup.mm10_CAST_Nmasked.seq.n_occur_gt1$Freq) - nrow(dedup.mm10_CAST_Nmasked.seq.n_occur_gt1)  # [1] 409440
# 
# # -----------------
# dedup.mm10_CAST_Nmasked.seq.n_seq_gt1 <- dedup.mm10_CAST_Nmasked[dedup.mm10_CAST_Nmasked$seq %in% dedup.mm10_CAST_Nmasked.seq.n_occur_gt1$Var1, ]
# dedup.mm10_CAST_Nmasked.seq.n_seq_gt1 %>% head()
# dedup.mm10_CAST_Nmasked.seq.n_seq_gt1 %>% nrow()  # [1] 748417
# 
# dedup.mm10_CAST_Nmasked.seq.n_seq_gt1$seq %>% unique() %>% length()  # [1] 338977
# length(dedup.mm10_CAST_Nmasked.seq.n_seq_gt1$seq) - length(unique(dedup.mm10_CAST_Nmasked.seq.n_seq_gt1$seq))  # [1] 409440

# -----------------
b_f.dedup.mm10_CAST_Nmasked <- dedup.mm10_CAST_Nmasked %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    )
# b_f.dedup.mm10_CAST_Nmasked$b_f %>% unique() %>% length()  # [1] 6037331
# nrow(dedup.mm10_CAST_Nmasked)  # [1] 6037332
# 
# nrow(dedup.mm10_CAST_Nmasked) - length(unique(b_f.dedup.mm10_CAST_Nmasked$b_f))  # [1] 1

# -------------------------------------
# filter.predup.mm10_CAST_Nmasked.seq -
# -----------------
# filter.predup.mm10_CAST_Nmasked.seq.n_occur <- data.frame(table(filter.predup.mm10_CAST_Nmasked$seq))
# filter.predup.mm10_CAST_Nmasked.seq.n_occur %>% head()
# filter.predup.mm10_CAST_Nmasked.seq.n_occur %>% nrow()  # [1] 5615296
# 
# nrow(filter.predup.mm10_CAST_Nmasked) - nrow(filter.predup.mm10_CAST_Nmasked.seq.n_occur)  # [1] 351614
# filter.predup.mm10_CAST_Nmasked$seq %>% unique() %>% length()  # [1] 5615296
# 
# # -----------------
# filter.predup.mm10_CAST_Nmasked.seq.n_occur_gt1 <- filter.predup.mm10_CAST_Nmasked.seq.n_occur[filter.predup.mm10_CAST_Nmasked.seq.n_occur$Freq > 1, ]
# filter.predup.mm10_CAST_Nmasked.seq.n_occur_gt1 %>% head()
# filter.predup.mm10_CAST_Nmasked.seq.n_occur_gt1 %>% nrow()  # [1] 289697
# filter.predup.mm10_CAST_Nmasked.seq.n_occur_gt1$Freq %>% sum()  # [1] 641311
# 
# sum(filter.predup.mm10_CAST_Nmasked.seq.n_occur_gt1$Freq) - nrow(filter.predup.mm10_CAST_Nmasked.seq.n_occur_gt1)  # [1] 351614
# 
# # -----------------
# filter.predup.mm10_CAST_Nmasked.seq.n_seq_gt1 <- filter.predup.mm10_CAST_Nmasked[filter.predup.mm10_CAST_Nmasked$seq %in% filter.predup.mm10_CAST_Nmasked.seq.n_occur_gt1$Var1, ]
# filter.predup.mm10_CAST_Nmasked.seq.n_seq_gt1 %>% head()
# filter.predup.mm10_CAST_Nmasked.seq.n_seq_gt1 %>% nrow()  # [1] 641311
# 
# filter.predup.mm10_CAST_Nmasked.seq.n_seq_gt1$seq %>% unique() %>% length()  # [1] 289697
# length(filter.predup.mm10_CAST_Nmasked.seq.n_seq_gt1$seq) - length(unique(filter.predup.mm10_CAST_Nmasked.seq.n_seq_gt1$seq))  # [1] 351614

# -----------------
# b_f.filter.predup.mm10_CAST_Nmasked <- filter.predup.mm10_CAST_Nmasked %>%
#     tidyr::unite(
#         b_f,
#         c("qname", "seq"),
#         sep = "_",
#         remove = FALSE
#     )
# b_f.filter.predup.mm10_CAST_Nmasked$b_f %>% unique() %>% length()  # [1] 5966909
# nrow(filter.predup.mm10_CAST_Nmasked)  # [1] 5966910
# 
# nrow(filter.predup.mm10_CAST_Nmasked) - length(unique(b_f.filter.predup.mm10_CAST_Nmasked$b_f))  # [1] 1

# -------------------------------------
# dedup.mm10.seq ----------------------
# -----------------
# dedup.mm10.seq.n_occur <- data.frame(table(dedup.mm10$seq))
# dedup.mm10.seq.n_occur %>% head()
# dedup.mm10.seq.n_occur %>% nrow()  # [1] 5654808
# 
# nrow(dedup.mm10) - nrow(dedup.mm10.seq.n_occur)  # [1] 411832
# dedup.mm10$seq %>% unique() %>% length()  # [1] 5654808
# 
# # -----------------
# dedup.mm10.seq.n_occur_gt1 <- dedup.mm10.seq.n_occur[dedup.mm10.seq.n_occur$Freq > 1, ]
# dedup.mm10.seq.n_occur_gt1 %>% head()
# dedup.mm10.seq.n_occur_gt1 %>% nrow()  # [1] 340377
# dedup.mm10.seq.n_occur_gt1$Freq %>% sum()  # [1] 752209
# 
# sum(dedup.mm10.seq.n_occur_gt1$Freq) - nrow(dedup.mm10.seq.n_occur_gt1)  # [1] 411832
# 
# # -----------------
# dedup.mm10.seq.n_seq_gt1 <- dedup.mm10[dedup.mm10$seq %in% dedup.mm10.seq.n_occur_gt1$Var1, ]
# dedup.mm10.seq.n_seq_gt1 %>% head()
# dedup.mm10.seq.n_seq_gt1 %>% nrow()  # [1] 752209
# 
# dedup.mm10.seq.n_seq_gt1$seq %>% unique() %>% length()  # [1] 340377
# length(dedup.mm10.seq.n_seq_gt1$seq) - length(unique(dedup.mm10.seq.n_seq_gt1$seq))  # [1] 411832

# -----------------
b_f.dedup.mm10 <- dedup.mm10 %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    )
# b_f.dedup.mm10$b_f %>% unique() %>% length()  # [1] 6066640
# nrow(dedup.mm10)  # [1] 6066640
# 
# nrow(dedup.mm10) - length(unique(b_f.dedup.mm10$b_f))  # [1] 0

# -------------------------------------
# filter.predup.mm10.seq ---------------------
# -----------------
# filter.predup.mm10.seq.n_occur <- data.frame(table(filter.predup.mm10$seq))
# filter.predup.mm10.seq.n_occur %>% head()
# filter.predup.mm10.seq.n_occur %>% nrow()  # [1] 5642015
# 
# nrow(filter.predup.mm10) - nrow(filter.predup.mm10.seq.n_occur)  # [1] 353749
# filter.predup.mm10$seq %>% unique() %>% length()  # [1] 5642015
# 
# # -----------------
# filter.predup.mm10.seq.n_occur_gt1 <- filter.predup.mm10.seq.n_occur[filter.predup.mm10.seq.n_occur$Freq > 1, ]
# filter.predup.mm10.seq.n_occur_gt1 %>% head()
# filter.predup.mm10.seq.n_occur_gt1 %>% nrow()  # [1] 290837
# filter.predup.mm10.seq.n_occur_gt1$Freq %>% sum()  # [1] 644586
# 
# sum(filter.predup.mm10.seq.n_occur_gt1$Freq) - nrow(filter.predup.mm10.seq.n_occur_gt1)  # [1] 353749
# 
# # -----------------
# filter.predup.mm10.seq.n_seq_gt1 <- filter.predup.mm10[filter.predup.mm10$seq %in% filter.predup.mm10.seq.n_occur_gt1$Var1, ]
# filter.predup.mm10.seq.n_seq_gt1 %>% head()
# filter.predup.mm10.seq.n_seq_gt1 %>% nrow()  # [1] 644586
# 
# filter.predup.mm10.seq.n_seq_gt1$seq %>% unique() %>% length()  # [1] 290837
# length(filter.predup.mm10.seq.n_seq_gt1$seq) - length(unique(filter.predup.mm10.seq.n_seq_gt1$seq))  # [1] 353749

# -----------------
# b_f.filter.predup.mm10 <- filter.predup.mm10 %>%
#     tidyr::unite(
#         b_f,
#         c("qname", "seq"),
#         sep = "_",
#         remove = FALSE
#     )
# b_f.filter.predup.mm10$b_f %>% unique() %>% length()  # [1] 5995764
# nrow(filter.predup.mm10)  # [1] 5995764
# 
# nrow(filter.predup.mm10) - length(unique(b_f.filter.predup.mm10$b_f))  # [1] 0


# -----------------------------------------------------------------------------
# setDifference <- function(x, y) {
#     c(setdiff(x, y), setdiff(y, x)) %>% return()
# }

d.C.b_f <- b_f.dedup.CAST$b_f
d.C.b_f %>% length()  # [1] 5738934
d.C.b_f %>% unique() %>% length()  # [1] 5738932

d.Mn.b_f <- b_f.dedup.mm10_CAST_Nmasked$b_f
d.Mn.b_f %>% length()  # [1] 6037332
d.Mn.b_f %>% unique() %>% length()  # [1] 6037331

d.M.b_f <- b_f.dedup.mm10$b_f
d.M.b_f %>% length()  # [1] 6066640
d.M.b_f %>% unique() %>% length()  # [1] 6066640

sd.C.Mn <- setdiff(d.C.b_f, d.Mn.b_f)
sd.C.M <- setdiff(d.C.b_f, d.M.b_f)
sd.Mn.M <- setdiff(d.Mn.b_f, d.M.b_f)
sd.Mn.C <- setdiff(d.Mn.b_f, d.C.b_f)
sd.M.C <- setdiff(d.M.b_f, d.C.b_f)
sd.M.Mn <- setdiff(d.M.b_f, d.Mn.b_f)

sd.C.M %>% head() %>% as_tibble()
sd.C.M %>% length()  # [1] 622828
sd.C.M %>% unique() %>% length()  # [1] 622828

sd.C.Mn %>% head() %>% as_tibble()
sd.C.Mn %>% length()  # [1] 619253
sd.C.Mn %>% unique() %>% length()  # [1] 619253

sd.M.C %>% head() %>% as_tibble()
sd.M.C %>% length()  # [1] 950536
sd.M.C %>% unique() %>% length()  # [1] 950536

sd.M.Mn %>% head() %>% as_tibble()
sd.M.Mn %>% length()  # [1] 314362
sd.M.Mn %>% unique() %>% length()  # [1] 314362

sd.Mn.M %>% head() %>% as_tibble()
sd.Mn.M %>% length()  # [1] 285053
sd.Mn.M %>% unique() %>% length()  # [1] 285053

sd.Mn.C %>% head() %>% as_tibble()
sd.Mn.C %>% length()  # [1] 917652
sd.Mn.C %>% unique() %>% length()  # [1] 917652


check_in.CAST <- b_f.dedup.CAST[d.C.b_f %in% d.M.b_f, ]  # 5116104
check_notin.CAST <- b_f.dedup.CAST[d.C.b_f %notin% d.M.b_f, ]  # 622830

check_in.mm10 <- b_f.dedup.mm10[d.M.b_f %in% d.C.b_f, ]  # 5116104
check_notin.mm10 <- b_f.dedup.mm10[d.M.b_f %notin% d.C.b_f, ]  # 950536

length(check_notin.CAST$b_f) + length(check_notin.mm10$b_f) + length(check_in.CAST$b_f)  # [1] 6689470


# -----------------------------------------------------------------------------
#  https://stat545.com/join-cheatsheet.html#full_joinsuperheroes-publishers

zC <- b_f.dedup.CAST %>% select(b_f, tag.AS)
zM <- b_f.dedup.mm10 %>% select(b_f, tag.AS)
zMn <- b_f.dedup.mm10_CAST_Nmasked %>% select(b_f, tag.AS)

zC <- zC %>%
    dplyr::rename(
        CAST.AS = tag.AS,
        # CAST.XS = tag.XS
    )
zM <- zM %>%
    dplyr::rename(
        mm10.AS = tag.AS,
        # mm10.XS = tag.XS
    )
zMn <- zMn %>%
    dplyr::rename(
        mm10_CAST_Nmasked.AS = tag.AS,
        # mm10.XS = tag.XS
    )
zC.zM <- full_join(zC, zM, by = "b_f")
zC.zMn <- full_join(zC, zMn, by = "b_f")
zM.zMn <- full_join(zM, zMn, by = "b_f")
zC.zM.zMn <- full_join(zC, zM, by = "b_f") %>% 
    full_join(., zMn, by = "b_f")
zMn.zC.zM <- full_join(zMn, zC, by = "b_f") %>% 
    full_join(., zM, by = "b_f")

# -----------------
test1 <- zMn.zC.zM %>% 
    mutate(difference = CAST.AS - mm10.AS)

int <- 0 %>% as.integer()  #TODO Make the integer an argument

# -----------------
test2 <- test1 %>% 
    dplyr::mutate(
        assignment.initial = case_when(
            difference >= (-1 * int) & difference <= int ~ "Neutral",
            difference > int ~ "CAST-EiJ",
            difference < (-1 * int) ~ "mm10"
        )
    )

test2$assignment <- ifelse(
    is.na(test2$assignment.initial),
    ifelse(
        !is.na(test2$CAST.AS),
        "CAST-EiJ",
        ifelse(
            !is.na(test2$mm10.AS),
            "mm10",
            test2$assignment.initial
        )
    ),
    test2$assignment.initial
) %>% forcats::as_factor()

# -----------------
test3 <- test2 %>% 
    select(b_f, mm10_CAST_Nmasked.AS, CAST.AS, mm10.AS, assignment)

test3$tmp.mm10_CAST_Nmasked <- ifelse(
    is.na(test3$mm10_CAST_Nmasked.AS), "0", "1"
)

test3$tmp.CAST <- ifelse(
    is.na(test3$CAST.AS), "0", "1"
)

test3$tmp.mm10 <- ifelse(
    is.na(test3$mm10.AS), "0", "1"
)

test3$trinary <- paste0(test3$tmp.mm10_CAST_Nmasked, test3$tmp.CAST, test3$tmp.mm10) %>% forcats::as_factor()

test3$assignment_trinary <- paste(test3$assignment, test3$trinary) %>% as_factor()

# -----------------
test4 <- test3 %>%
    select(
        b_f,
        mm10_CAST_Nmasked.AS,
        CAST.AS,
        mm10.AS,
        assignment,
        trinary,
        assignment_trinary
    )

#  Make figures
# ------------------
ggplot(test4, aes(x = trinary)) +
    geom_bar(alpha = 0.5) +
    geom_text(stat = "count", aes(label = scales::comma(..count..)), size = 3, vjust = -0.5) +
    ylab("") +
    xlab("assembly-alignment combination") +
    scale_y_continuous(labels = comma)

ggplot(test4, aes(x = assignment)) +
    geom_bar(alpha = 0.5) +
    geom_text(stat = "count", aes(label = scales::comma(..count..)), size = 3, vjust = -0.5) +
    ylab("") +
    xlab("alignment assembly assignment") +
    scale_y_continuous(labels = comma)

ggplot(test4, aes(x = assignment_trinary)) +
    geom_bar(alpha = 0.5) +
    geom_text(stat = "count", aes(label = scales::comma(..count..)), size = 2.5, vjust = -0.5) +
    ylab("") +
    xlab("assembly-alignment assignment") +
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -----------------
ggplot(test4, aes(x = assignment, color = trinary, fill = trinary)) +
    geom_bar(alpha = 0.5, width = 0.9) +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

#  Cool
ggplot(test4, aes(x = assignment, color = trinary, fill = trinary)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

ggplot(test4, aes(x = trinary, color = assignment, fill = assignment)) +
    geom_bar(alpha = 0.5, width = 0.9) +
    ylab("") +
    xlab("assembly-alignment combination") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

#  Cool
ggplot(test4, aes(x = trinary, color = assignment, fill = assignment)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ylab("") +
    xlab("assembly-alignment combination") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

# -----------------
# -------
ggplot(test4, aes(x = assignment, fill = as_factor(CAST.AS))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle("AS from alignment to CAST") +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

ggplot(test4, aes(x = assignment, fill = as_factor(CAST.AS))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle("AS from alignment to CAST") +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma) +
    facet_grid(trinary ~ .)

# -------
ggplot(test4, aes(x = assignment, fill = as_factor(mm10.AS))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle("AS from alignment to mm10") +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

ggplot(test4, aes(x = assignment, fill = as_factor(mm10.AS))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle("AS from alignment to mm10") +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma) +
    facet_grid(trinary ~ .)

# -------
ggplot(test4, aes(x = assignment, fill = as_factor(mm10_CAST_Nmasked.AS))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle("AS from alignment to mm10-CAST-N-masked") +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

ggplot(test4, aes(x = assignment, fill = as_factor(mm10_CAST_Nmasked.AS))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle("AS from alignment to mm10-CAST-N-masked") +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma) +
    facet_grid(trinary ~ .)

# -----------------
# ggplot(test4, aes(x = CAST.AS, color = assignment)) +
#     geom_histogram(binwidth = 1, fill = "white", alpha = 0.5)
# ggplot(test4, aes(x = mm10.AS, color = assignment)) +
#     geom_histogram(binwidth = 1, fill = "white", alpha = 0.5)
# ggplot(test4, aes(x = mm10_CAST_Nmasked.AS, color = assignment)) +
#     geom_histogram(binwidth = 1, fill = "white", alpha = 0.5)

ggplot(test4, aes(x = CAST.AS, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to CAST")
ggplot(test4, aes(x = mm10.AS, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to mm10")
ggplot(test4, aes(x = mm10_CAST_Nmasked.AS, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to mm10-CAST-N-masked")

ggplot(test4, aes(x = CAST.AS, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to CAST") +
    facet_grid(assignment ~ .)
ggplot(test4, aes(x = mm10.AS, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to mm10") +
    facet_grid(assignment ~ .)
ggplot(test4, aes(x = mm10_CAST_Nmasked.AS, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to mm10-CAST-N-masked") +
    facet_grid(assignment ~ .)
