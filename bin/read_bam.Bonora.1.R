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
"/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora" %>% setwd()
# getwd()
# list.files()
# list.dirs()


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
(dedup.files <- list.files(pattern = "\\.dedup.bam$"))
(mark.files <- list.files(pattern = "\\Duplicates.sort.bam$"))
(files <- c(dedup.files, mark.files))
n <- length(files)
variables <- paste0("bam", 1:n)
variables_rname <- paste0("bam", 1:n, "$rname")
command_pipe <- paste0("<- bam", 1:n, " %>% ")
chromosomes <- c(paste0("chr", 1:19), "chrX")

mapply(
    assign, variables, files, MoreArgs = list(envir = parent.frame())
)

#  Assign .bam information as list
dedup.129 <- bam1 %>% Rsamtools::scanBam(., param = map_params)
dedup.CAST <- bam2 %>% Rsamtools::scanBam(., param = map_params)
dedup.mm10_CAST_129_Nmasked <- bam3 %>% Rsamtools::scanBam(., param = map_params)
dedup.mm10 <- bam4 %>% Rsamtools::scanBam(., param = map_params)
mark.129 <- bam5 %>% Rsamtools::scanBam(., param = map_params)
mark.CAST <- bam6 %>% Rsamtools::scanBam(., param = map_params)
mark.predup.CAST <- bam7 %>% Rsamtools::scanBam(., param = map_params)
mark.mm10_CAST_129_Nmasked <- bam8 %>% Rsamtools::scanBam(., param = map_params)
mark.mm10 <- bam9 %>% Rsamtools::scanBam(., param = map_params)

#  Convert .bam information from list to dataframe to tibble
variables <- c(
    "dedup.129",
    "dedup.CAST",
    "dedup.mm10_CAST_129_Nmasked",
    "dedup.mm10",
    "mark.129",
    "mark.CAST",
    "mark.predup.CAST",
    "mark.mm10_CAST_129_Nmasked",
    "mark.mm10"
)
command_pipe <- paste0("<- ", variables, " %>% ")

command <- paste0(command_pipe, "as.data.frame() %>% as_tibble()")
operation <- makeOperation(variables, command)
eval(parse(text = operation))
rm("bam1", "bam2", "bam3", "bam4", "bam5", "bam6", "bam7", "bam8", "bam9")

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

# dedup.129 %>% head()
# mark.129 %>% head()
dedup.CAST %>% head()
mark.CAST %>% head()
mark.predup.CAST %>% head()
# dedup.mm10_CAST_129_Nmasked %>% head()
# mark.mm10_CAST_129_Nmasked %>% head()
# dedup.mm10 %>% head()
# mark.mm10 %>% head()


#  Numbers of unique qname (barcodes)
dedup.129.qname.unique <- dedup.129$qname %>% unique() %>% length()  # [1] 2972
length(dedup.129$qname) - dedup.129.qname.unique  # [1] 13167654

dedup.CAST.qname.unique <- dedup.CAST$qname %>% unique() %>% length()  # [1] 2961
length(dedup.CAST$qname) - dedup.CAST.qname.unique  # [1] 13031773

mark.CAST.qname.unique <- mark.CAST$qname %>% unique() %>% length()  # [1] 2961
length(mark.CAST$qname) - mark.CAST.qname.unique  # [1] 13031773

mark.predup.CAST.qname.unique <- mark.predup.CAST$qname %>% unique() %>% length()  # [1] 3091
length(mark.predup.CAST$qname) - mark.predup.CAST.qname.unique  # [1] 13074289

dedup.mm10_CAST_129_Nmasked.qname.unique <- dedup.mm10_CAST_129_Nmasked$qname %>% unique() %>% length()  # [1] 2945
length(dedup.mm10_CAST_129_Nmasked$qname) - dedup.mm10_CAST_129_Nmasked.qname.unique  # [1] 12947283

dedup.mm10.qname.unique <- dedup.mm10$qname %>% unique() %>% length()  # [1] 2959
length(dedup.mm10$qname) - dedup.mm10.qname.unique  # [1] 13022048

#  Numbers of unique seq (fragments)
dedup.129.seq.unique <- dedup.129$seq %>% unique() %>% length()  # [1] 6326297
length(dedup.129$seq) - dedup.129.seq.unique  # [1] 6844329

dedup.CAST.seq.unique <- dedup.CAST$seq %>% unique() %>% length()  # [1] 6314847
length(dedup.CAST$seq) - dedup.CAST.seq.unique  # [1] 6719887

mark.CAST.seq.unique <- mark.CAST$seq %>% unique() %>% length()  # [1] 6314847
length(mark.CAST$seq) - mark.CAST.seq.unique  # 6719887

mark.predup.CAST.seq.unique <- mark.predup.CAST$seq %>% unique() %>% length()  # [1] 6326982
length(mark.predup.CAST$seq) - mark.predup.CAST.seq.unique  # 6750398

dedup.mm10_CAST_129_Nmasked.seq.unique <- dedup.mm10_CAST_129_Nmasked$seq %>% unique() %>% length()  # [1] 6337920
length(dedup.mm10_CAST_129_Nmasked$seq) - dedup.mm10_CAST_129_Nmasked.seq.unique  # [1] 6612308

dedup.mm10.seq.unique <- dedup.mm10$seq %>% unique() %>% length()  # [1] 6389895
length(dedup.mm10$seq) - dedup.mm10.seq.unique  # [1] 6635112


#  Check distributions of .bam flags ------------------------------------------
dedup.129$flag %>% table()
# .
# 83      99     147     163 
# 3291077 3294232 3294289 3291028

mark.129$flag %>% table()
# .
# 83      99     147     163    1107    1123    1171    1187 
# 2509245 2509951 2512854 2510982  781832  784281  781435  780046
# 2509245 + 781832  # [1] 3291077
# 2509951 + 784281  # [1] 3294232
# 2512854 + 781435  # [1] 3294289
# 2510982 + 780046  # [1] 3291028

dedup.CAST$flag %>% table()
# .
# 83      99     147     163 
# 3255709 3261656 3261696 3255673

mark.CAST$flag %>% table()
# .
# 83      99     147     163    1107    1123    1171    1187 
# 2503907 2506250 2507930 2505171  751802  755406  753766  750502

mark.predup.CAST$flag %>% table()
# .
# 83      99     147     163    1107    1123    1171    1187 
# 2506401 2508778 2510434 2507695  759953  763558  761902  758659

dedup.mm10_CAST_129_Nmasked$flag %>% table()
# .
# 83      99     147     163 
# 3235476 3239640 3239669 3235443

mark.mm10_CAST_129_Nmasked$flag %>% table()
# .
# 83      99     147     163    1107    1123    1171    1187 
# 2516917 2518993 2518571 2516503  718559  720647  721098  718940

dedup.mm10$flag %>% table()
# .
# 83      99     147     163 
# 3254309 3258198 3258226 3254274

mark.mm10$flag %>% table()
# .
# 83      99     147     163    1107    1123    1171    1187 
# 2527805 2529973 2529390 2527261  726504  728225  728836  727013


#  Remove "duplicate" status rows ---------------------------------------------
`%notin%` <- Negate(`%in%`)
duplicate_flags <- c(1107, 1123, 1171, 1187)
filter.mark.129 <- mark.129 %>% dplyr::filter(flag %notin% duplicate_flags)
filter.mark.129$flag %>% table()
# .
# 83      99     147     163 
# 2509245 2509951 2512854 2510982

filter.mark.CAST <- mark.CAST %>% dplyr::filter(flag %notin% duplicate_flags)
filter.mark.CAST$flag %>% table()
# .
# 83      99     147     163 
# 2503907 2506250 2507930 2505171

filter.mark.mm10_CAST_129_Nmasked <- mark.mm10_CAST_129_Nmasked %>% dplyr::filter(flag %notin% duplicate_flags)
filter.mark.mm10_CAST_129_Nmasked$flag %>% table()
# .
# 83      99     147     163 
# 2516917 2518993 2518571 2516503

filter.mark.mm10 <- mark.mm10 %>% dplyr::filter(flag %notin% duplicate_flags)
filter.mark.mm10$flag %>% table()
# .
# 83      99     147     163 
# 2527805 2529973 2529390 2527261


#  Numbers of unique qname (barcodes)
filter.mark.129.qname.unique <- filter.mark.129$qname %>% unique() %>% length()  # [1] 2878
length(filter.mark.129$qname) - filter.mark.129.qname.unique  # [1] 10040154

filter.mark.CAST.qname.unique <- filter.mark.CAST$qname %>% unique() %>% length()  # [1] 2864
length(filter.mark.CAST$qname) - filter.mark.CAST.qname.unique  # [1] 10020394

filter.mark.mm10_CAST_129_Nmasked.qname.unique <- filter.mark.mm10_CAST_129_Nmasked$qname %>% unique() %>% length()  # [1] 2856
length(filter.mark.mm10_CAST_129_Nmasked$qname) - filter.mark.mm10_CAST_129_Nmasked.qname.unique  # [1] 10068128

filter.mark.mm10.qname.unique <- filter.mark.mm10$qname %>% unique() %>% length()  # [1] 2871
length(filter.mark.mm10$qname) - filter.mark.mm10.qname.unique  # [1] 10111558

#  Numbers of unique seq (fragments)
filter.mark.129.seq.unique <- filter.mark.129$seq %>% unique() %>% length()  # [1] 5929533
length(filter.mark.129$seq) - filter.mark.129.seq.unique  # [1] 4113499

filter.mark.CAST.seq.unique <- filter.mark.CAST$seq %>% unique() %>% length()  # [1] 5928689
length(filter.mark.CAST$seq) - filter.mark.CAST.seq.unique  # [1] 4094569

filter.mark.mm10_CAST_129_Nmasked.seq.unique <- filter.mark.mm10_CAST_129_Nmasked$seq %>% unique() %>% length()  # [1] 5975635
length(filter.mark.mm10_CAST_129_Nmasked$seq) - filter.mark.mm10_CAST_129_Nmasked.seq.unique  # [1] 4095349

filter.mark.mm10.seq.unique <- filter.mark.mm10$seq %>% unique() %>% length()  # [1] 6011517
length(filter.mark.mm10$seq) - filter.mark.mm10.seq.unique  # [1] 4102912


# -----------------------------------------------------------------------------
# dedup.129.seq (tidyr::unite as in read_bam.Disteche.R -----------------------
# -----------------------------------------------------------------------------
dedup.129.seq.n_occur <- data.frame(table(dedup.129$seq))  #TODO Faster way to do this
dedup.129.seq.n_occur %>% head()
dedup.129.seq.n_occur %>% nrow()  # [1] 6326297

nrow(dedup.129) - nrow(dedup.129.seq.n_occur)  # [1] 6844329
dedup.129$seq %>% unique() %>% length()  # [1] 6326297

# -----------------
dedup.129.seq.n_occur_gt1 <- dedup.129.seq.n_occur[dedup.129.seq.n_occur$Freq > 1, ]
dedup.129.seq.n_occur_gt1 %>% head()
dedup.129.seq.n_occur_gt1 %>% nrow()  # [1] 3041010
dedup.129.seq.n_occur_gt1$Freq %>% sum()  # [1] 9885339

sum(dedup.129.seq.n_occur_gt1$Freq) - nrow(dedup.129.seq.n_occur_gt1)  # [1] 6844329

# -----------------
dedup.129.seq.n_seq_gt1 <- dedup.129[dedup.129$seq %in% dedup.129.seq.n_occur_gt1$Var1, ]
dedup.129.seq.n_seq_gt1 %>% head()
dedup.129.seq.n_seq_gt1 %>% nrow()  # [1] 9885339

dedup.129.seq.n_seq_gt1$seq %>% unique() %>% length()  # [1] 3041010
length(dedup.129.seq.n_seq_gt1$seq) - length(unique(dedup.129.seq.n_seq_gt1$seq))  # [1] 6844329

# -----------------
dedup.129.barcode_fragment <- dedup.129 %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    )
dedup.129.barcode_fragment$b_f %>% unique() %>% length()  # [1] 7675828
nrow(dedup.129)  # [1] 13170626

nrow(dedup.129) - length(unique(dedup.129.barcode_fragment$b_f))  # [1] 5494798

dedup.129.b_f.n_occur <- dedup.129.barcode_fragment$b_f %>%
    table() %>%
    data.frame()
dedup.129.b_f.n_gt1 <- dedup.129.b_f.n_occur[dedup.129.b_f.n_occur$Freq > 1, ]
dedup.129.b_f.n_seq_gt1 <- dedup.129.barcode_fragment[dedup.129.barcode_fragment$b_f %in% dedup.129.b_f.n_gt1$., ]

nrow(dedup.129.b_f.n_gt1)  # [1] 2863367
nrow(dedup.129.b_f.n_seq_gt1)  # [1] 8358165

# -----------------------------------------------------------------------------
# dedup.129.seq (different tidyr::unite) --------------------------------------
# -----------------------------------------------------------------------------
dedup.129.barcode_fragment_etc <- dedup.129 %>%
    tidyr::unite(
        b_f,
        c("qname", "rname", "pos", "seq"),
        sep = "_",
        remove = FALSE
    )
dedup.129.barcode_fragment_etc$b_f %>% unique() %>% length()  # [1] 7680113
nrow(dedup.129)  # [1] 13170626

nrow(dedup.129) - length(unique(dedup.129.barcode_fragment_etc$b_f))  # [1] 5490513

dedup.129.b_f_etc.n_occur <- dedup.129.barcode_fragment_etc$b_f %>%
    table() %>%
    data.frame()
dedup.129.b_f_etc.n_gt1 <- dedup.129.b_f_etc.n_occur[dedup.129.b_f_etc.n_occur$Freq > 1, ]
dedup.129.b_f_etc.n_seq_gt1 <- dedup.129.barcode_fragment_etc[dedup.129.barcode_fragment_etc$b_f %in% dedup.129.b_f_etc.n_gt1$., ]

nrow(dedup.129.b_f_etc.n_gt1)  # [1] 2862096
nrow(dedup.129.b_f_etc.n_seq_gt1)  # [1] 8352609

# -----------------------------------------------------------------------------
# Some legit deduping ---------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------
# scrub.129 <- dedup.129 %>% dplyr::distinct(qname, rname, pos, seq, .keep_all = TRUE)
# 
# scrub.129.barcode_fragment_etc <- scrub.129 %>%
#     tidyr::unite(
#         b_f,
#         c("qname", "seq"),
#         sep = "_",
#         remove = FALSE
#     )
# scrub.129.barcode_fragment_etc$b_f %>% unique() %>% length()  # [1] 7675828
# nrow(scrub.129)  # [1] 7680113
# 
# nrow(scrub.129) - length(unique(scrub.129.barcode_fragment_etc$b_f))  # [1] 4285
# 
# scrub.129.b_f_etc.n_occur <- scrub.129.barcode_fragment_etc$b_f %>%
#     table() %>%
#     data.frame()
# scrub.129.b_f_etc.n_gt1 <- scrub.129.b_f_etc.n_occur[scrub.129.b_f_etc.n_occur$Freq > 1, ]
# scrub.129.b_f_etc.n_seq_gt1 <- scrub.129.barcode_fragment_etc[scrub.129.barcode_fragment_etc$b_f %in% scrub.129.b_f_etc.n_gt1$., ] %>% dplyr::arrange(b_f, rname, pos)
# View(scrub.129.b_f_etc.n_seq_gt1)
# 
# nrow(scrub.129.b_f_etc.n_gt1)  # [1] 4119
# nrow(scrub.129.b_f_etc.n_seq_gt1)  # [1] 8404

# -----------------
scrub.129A <- dedup.129 %>% dplyr::distinct(qname, seq, .keep_all = TRUE)

scrub.129A.barcode_fragment_etc <- scrub.129A %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    )
scrub.129A.barcode_fragment_etc$b_f %>% unique() %>% length()  # [1] 7675828
nrow(scrub.129A)  # [1] 7675828

nrow(scrub.129A) - length(unique(scrub.129A.barcode_fragment_etc$b_f))  # [1] 0

scrub.129A.b_f_etc.n_occur <- scrub.129A.barcode_fragment_etc$b_f %>%
    table() %>%
    data.frame()
scrub.129A.b_f_etc.n_gt1 <- scrub.129A.b_f_etc.n_occur[scrub.129A.b_f_etc.n_occur$Freq > 1, ]
scrub.129A.b_f_etc.n_seq_gt1 <- scrub.129A.barcode_fragment_etc[scrub.129A.barcode_fragment_etc$b_f %in% scrub.129A.b_f_etc.n_gt1$., ] %>% dplyr::arrange(b_f, rname, pos)
View(scrub.129A.b_f_etc.n_seq_gt1)

nrow(scrub.129A.b_f_etc.n_gt1)  # [1] 0
nrow(scrub.129A.b_f_etc.n_seq_gt1)  # [1] 0

# -----------------
scrub.129B <- dedup.129 %>% dplyr::distinct(qname, rname, pos, seq, .keep_all = TRUE)

scrub.129B.barcode_fragment_etc <- scrub.129B %>%
    tidyr::unite(
        b_f,
        c("qname", "rname", "pos", "seq"),
        sep = "_",
        remove = FALSE
    )
scrub.129B.barcode_fragment_etc$b_f %>% unique() %>% length()  # [1] 7680113
nrow(scrub.129B)  # [1] 7680113

nrow(scrub.129B) - length(unique(scrub.129B.barcode_fragment_etc$b_f))  # [1] 0

scrub.129B.b_f_etc.n_occur <- scrub.129B.barcode_fragment_etc$b_f %>%
    table() %>%
    data.frame()
scrub.129B.b_f_etc.n_gt1 <- scrub.129B.b_f_etc.n_occur[scrub.129B.b_f_etc.n_occur$Freq > 1, ]
scrub.129B.b_f_etc.n_seq_gt1 <- scrub.129B.barcode_fragment_etc[scrub.129B.barcode_fragment_etc$b_f %in% scrub.129B.b_f_etc.n_gt1$., ] %>% dplyr::arrange(b_f, rname, pos)
View(scrub.129B.b_f_etc.n_seq_gt1)

nrow(scrub.129B.b_f_etc.n_gt1)  # [1] 0
nrow(scrub.129B.b_f_etc.n_seq_gt1)  # [1] 0

# -----------------
scrub.129C <- dedup.129 %>% dplyr::distinct(qname, rname, pos, mpos, seq, .keep_all = TRUE)

scrub.129C.barcode_fragment_etc <- scrub.129C %>%
    tidyr::unite(
        b_f,
        c("qname", "rname", "pos", "mpos", "seq"),
        sep = "_",
        remove = FALSE
    )
scrub.129C.barcode_fragment_etc$b_f %>% unique() %>% length()  # [1] 7732477
nrow(scrub.129C)  # [1] 7732477

nrow(scrub.129C) - length(unique(scrub.129C.barcode_fragment_etc$b_f))  # [1] 0

scrub.129C.b_f_etc.n_occur <- scrub.129C.barcode_fragment_etc$b_f %>%
    table() %>%
    data.frame()
scrub.129C.b_f_etc.n_gt1 <- scrub.129C.b_f_etc.n_occur[scrub.129C.b_f_etc.n_occur$Freq > 1, ]
scrub.129C.b_f_etc.n_seq_gt1 <- scrub.129C.barcode_fragment_etc[scrub.129C.barcode_fragment_etc$b_f %in% scrub.129C.b_f_etc.n_gt1$., ] %>% dplyr::arrange(b_f, rname, pos)
View(scrub.129C.b_f_etc.n_seq_gt1)

nrow(scrub.129C.b_f_etc.n_gt1)  # [1] 0
nrow(scrub.129C.b_f_etc.n_seq_gt1)  # [1] 0


# -----------------------------------------------------------------------------
# filter.mark.129.seq (tidyr::unite as in read_bam.Disteche.R -----------------
# -----------------------------------------------------------------------------
filter.mark.129.seq.n_occur <- data.frame(table(filter.mark.129$seq))  #TODO Faster way to do this
filter.mark.129.seq.n_occur %>% head()
filter.mark.129.seq.n_occur %>% nrow()  # [1] 5929533

nrow(filter.mark.129) - nrow(filter.mark.129.seq.n_occur)  # [1] 4113499
filter.mark.129$seq %>% unique() %>% length()  # [1] 5929533

# -----------------
filter.mark.129.seq.n_occur_gt1 <- filter.mark.129.seq.n_occur[filter.mark.129.seq.n_occur$Freq > 1, ]
filter.mark.129.seq.n_occur_gt1 %>% head()
filter.mark.129.seq.n_occur_gt1 %>% nrow()  # [1] 2699347
filter.mark.129.seq.n_occur_gt1$Freq %>% sum()  # [1] 6812846

sum(filter.mark.129.seq.n_occur_gt1$Freq) - nrow(filter.mark.129.seq.n_occur_gt1)  # [1] 4113499

# -----------------
filter.mark.129.seq.n_seq_gt1 <- filter.mark.129[filter.mark.129$seq %in% filter.mark.129.seq.n_occur_gt1$Var1, ]
filter.mark.129.seq.n_seq_gt1 %>% head()
filter.mark.129.seq.n_seq_gt1 %>% nrow()  # [1] 6812846

filter.mark.129.seq.n_seq_gt1$seq %>% unique() %>% length()  # [1] 2699347
length(filter.mark.129.seq.n_seq_gt1$seq) - length(unique(filter.mark.129.seq.n_seq_gt1$seq))  # [1] 4113499

# -----------------
filter.mark.129.barcode_fragment <- filter.mark.129 %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    )
filter.mark.129.barcode_fragment$b_f %>% unique() %>% length()  # [1] 6585734
nrow(filter.mark.129)  # [1] 10043032

nrow(filter.mark.129) - length(unique(filter.mark.129.barcode_fragment$b_f))  # [1] 3457298

filter.mark.129.b_f.n_occur <- filter.mark.129.barcode_fragment$b_f %>%
    table() %>%
    data.frame()
filter.mark.129.b_f.n_gt1 <- filter.mark.129.b_f.n_occur[filter.mark.129.b_f.n_occur$Freq > 1, ]
filter.mark.129.b_f.n_seq_gt1 <- filter.mark.129.barcode_fragment[filter.mark.129.barcode_fragment$b_f %in% filter.mark.129.b_f.n_gt1$., ]

nrow(filter.mark.129.b_f.n_gt1)  # [1] 2650970
nrow(filter.mark.129.b_f.n_seq_gt1)  # [1] 6108268

# -----------------------------------------------------------------------------
# filter.mark.129.seq (different tidyr::unite) --------------------------------
# -----------------------------------------------------------------------------
filter.mark.129.barcode_fragment_etc <- filter.mark.129 %>%
    tidyr::unite(
        b_f,
        c("qname", "rname", "pos", "seq"),
        sep = "_",
        remove = FALSE
    )
filter.mark.129.barcode_fragment_etc$b_f %>% unique() %>% length()  # [1] 6588335
nrow(filter.mark.129)  # [1] 10043032

nrow(filter.mark.129) - length(unique(filter.mark.129.barcode_fragment_etc$b_f))  # [1] 3454697

filter.mark.129.b_f_etc.n_occur <- filter.mark.129.barcode_fragment_etc$b_f %>%
    table() %>%
    data.frame()
filter.mark.129.b_f_etc.n_gt1 <- filter.mark.129.b_f_etc.n_occur[filter.mark.129.b_f_etc.n_occur$Freq > 1, ]
filter.mark.129.b_f_etc.n_seq_gt1 <- filter.mark.129.barcode_fragment_etc[filter.mark.129.barcode_fragment_etc$b_f %in% filter.mark.129.b_f_etc.n_gt1$., ]

nrow(filter.mark.129.b_f_etc.n_gt1)  # [1] 2649376
nrow(filter.mark.129.b_f_etc.n_seq_gt1)  # [1] 6104073


# -----------------------------------------------------------------------------
# Some legit deduping ---------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------
scrub.mark.129A <- mark.129 %>% dplyr::distinct(qname, seq, .keep_all = TRUE)

scrub.mark.129A.barcode_fragment_etc <- scrub.mark.129A %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    )
scrub.mark.129A.barcode_fragment_etc$b_f %>% unique() %>% length()  # [1] 7675828
nrow(mark.129)
nrow(scrub.mark.129A)  # [1] 7675828

View(mark.129)

nrow(scrub.mark.129A) - length(unique(scrub.mark.129A.barcode_fragment_etc$b_f))  # [1] 0

scrub.mark.129A.b_f_etc.n_occur <- scrub.mark.129A.barcode_fragment_etc$b_f %>%
    table() %>%
    data.frame()
scrub.mark.129A.b_f_etc.n_gt1 <- scrub.mark.129A.b_f_etc.n_occur[scrub.mark.129A.b_f_etc.n_occur$Freq > 1, ]
scrub.mark.129A.b_f_etc.n_seq_gt1 <- scrub.mark.129A.barcode_fragment_etc[scrub.mark.129A.barcode_fragment_etc$b_f %in% scrub.mark.129A.b_f_etc.n_gt1$., ] %>% dplyr::arrange(b_f, rname, pos)
View(scrub.mark.129A.b_f_etc.n_seq_gt1)

nrow(scrub.mark.129A.b_f_etc.n_gt1)  # [1] 0
nrow(scrub.mark.129A.b_f_etc.n_seq_gt1)  # [1] 0

# -----------------
scrub.mark.129B <- mark.129 %>% dplyr::distinct(qname, rname, pos, seq, .keep_all = TRUE)

scrub.mark.129B.barcode_fragment_etc <- scrub.mark.129B %>%
    tidyr::unite(
        b_f,
        c("qname", "rname", "pos", "seq"),
        sep = "_",
        remove = FALSE
    )
scrub.mark.129B.barcode_fragment_etc$b_f %>% unique() %>% length()  # [1] 7680113
nrow(scrub.mark.129B)  # [1] 7680113

nrow(scrub.mark.129B) - length(unique(scrub.mark.129B.barcode_fragment_etc$b_f))  # [1] 0

scrub.mark.129B.b_f_etc.n_occur <- scrub.mark.129B.barcode_fragment_etc$b_f %>%
    table() %>%
    data.frame()
scrub.mark.129B.b_f_etc.n_gt1 <- scrub.mark.129B.b_f_etc.n_occur[scrub.mark.129B.b_f_etc.n_occur$Freq > 1, ]
scrub.mark.129B.b_f_etc.n_seq_gt1 <- scrub.mark.129B.barcode_fragment_etc[scrub.mark.129B.barcode_fragment_etc$b_f %in% scrub.mark.129B.b_f_etc.n_gt1$., ] %>% dplyr::arrange(b_f, rname, pos)
View(scrub.mark.129B.b_f_etc.n_seq_gt1)

nrow(scrub.mark.129B.b_f_etc.n_gt1)  # [1] 0
nrow(scrub.mark.129B.b_f_etc.n_seq_gt1)  # [1] 0

# -----------------
scrub.mark.129C <- mark.129 %>% dplyr::distinct(qname, rname, pos, mpos, seq, .keep_all = TRUE)

scrub.mark.129C.barcode_fragment_etc <- scrub.mark.129C %>%
    tidyr::unite(
        b_f,
        c("qname", "rname", "pos", "mpos", "seq"),
        sep = "_",
        remove = FALSE
    )
scrub.mark.129C.barcode_fragment_etc$b_f %>% unique() %>% length()  # [1] 7732477
nrow(scrub.mark.129C)  # [1] 7732477

nrow(scrub.mark.129C) - length(unique(scrub.mark.129C.barcode_fragment_etc$b_f))  # [1] 0

scrub.mark.129C.b_f_etc.n_occur <- scrub.mark.129C.barcode_fragment_etc$b_f %>%
    table() %>%
    data.frame()
scrub.mark.129C.b_f_etc.n_gt1 <- scrub.mark.129C.b_f_etc.n_occur[scrub.mark.129C.b_f_etc.n_occur$Freq > 1, ]
scrub.mark.129C.b_f_etc.n_seq_gt1 <- scrub.mark.129C.barcode_fragment_etc[scrub.mark.129C.barcode_fragment_etc$b_f %in% scrub.mark.129C.b_f_etc.n_gt1$., ] %>% dplyr::arrange(b_f, rname, pos)
View(scrub.mark.129C.b_f_etc.n_seq_gt1)

nrow(scrub.mark.129C.b_f_etc.n_gt1)  # [1] 0
nrow(scrub.mark.129C.b_f_etc.n_seq_gt1)  # [1] 0
