#!/usr/bin/env Rscript

library(BSgenome)
library(Mmusculus129S1SvImJ)
library(MmusculusCASTEiJ)
library(Mmusculus129Inserted)
library(MmusculusCASTInserted)
library(MmusculusCAST129Inserted)
library(Mmusculus129Nmasked)
library(MmusculusCASTNmasked)
library(MmusculusCAST129Nmasked)
library(stringr)
library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)


#  Run parts 1 and 2 first; load .Rdata into environment ----------------------
#  Temporary: Set up work directory (location TB∆)
directory_work <- "/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora"
setwd(directory_work)
rm(directory_work)

load("test.PE-processing.2021-1221.2022-0106.part-2.Rdata")
script <- "test.PE-processing.2021-1221.2022-0108.part-3.R"


#  Set up functions for liftOver conversions, etc. ----------------------------

#  BSgenome = Mmusculus129S1SvImJ
#+ BSgenome = MmusculusCAST129Nmasked
#+ BSgenome = MmusculusCAST129Inserted
compareEvenBSgenome <- function(tibble, row, BSgenome) {
    chr <- tibble$rname.even[row] %>% as.character()
    start <- tibble$pos.even[row]
    end <- tibble$pos_end.even[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- tibble$seq.even[row]  #TODO Fix this...
    string_2 <<- BSgenome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this...
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    comparison <- paste0(
        # paste0(chr, ":", start, "-", end, "\n"),
        p$string_1, "\n",
        p$result, "\n",
        p$string_2, "\n",
        p$n_diff, " difference(s)."
    )
    return(comparison)
}


compareOddBSgenome <- function(tibble, row, BSgenome) {
    chr <- tibble$rname.odd[row] %>% as.character()
    start <- tibble$pos.odd[row]
    end <- tibble$pos_end.odd[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- tibble$seq.odd[row]  #TODO Fix this...
    string_2 <<- BSgenome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this...
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    comparison <- paste0(
        # paste0(chr, ":", start, "-", end, "\n"),
        p$string_1, "\n",
        p$result, "\n",
        p$string_2, "\n",
        p$n_diff, " difference(s)."
    )
    return(comparison)
}


# comp129S1 <- function(bam, row, genome = Mmusculus129S1SvImJ) {
#     chr <- bam$rname.129S1[row] %>% as.character()
#     start <- bam$pos.129S1[row]
#     end <- bam$pos_end.129S1[row]
#     
#     gr <- GRanges(
#         seqnames = Rle(chr),
#         ranges = IRanges(start = start, end = end)
#     )
#     string_1 <<- bam$seq.129S1[row]  #TODO Fix this...
#     string_2 <<- genome %>%
#         BSgenome::getSeq(., names = gr) %>%
#         as.character()  #TODO Fix this...
#     
#     p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
#     comparison <- paste0(
#         # paste0(chr, ":", start, "-", end, "\n"),
#         p$string_1, "\n",
#         p$result, "\n",
#         p$string_2, "\n",
#         p$n_diff, " difference(s)."
#     )
#     return(comparison)
# }
# 
# 
# compNmasked <- function(bam, row, genome = MmusculusCAST129Nmasked) {
#     chr <- bam$rname.GB[row] %>% as.character()
#     start <- bam$pos.GB[row]
#     end <- bam$pos_end.GB[row]
#     
#     gr <- GRanges(
#         seqnames = Rle(chr),
#         ranges = IRanges(start = start, end = end)
#     )
#     string_1 <<- bam$seq.GB[row]  #TODO Fix this...
#     string_2 <<- genome %>%
#         BSgenome::getSeq(., names = gr) %>%
#         as.character()  #TODO Fix this...
#     
#     p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
#     comparison <- paste0(
#         # paste0(chr, ":", start, "-", end, "\n"),
#         p$string_1, "\n",
#         p$result, "\n",
#         p$string_2, "\n",
#         p$n_diff, " difference(s)."
#     )
#     return(comparison)
# }
# 
# 
# compNinserted <- function(bam, row, genome = MmusculusCAST129Inserted) {
#     chr <- bam$rname.GB[row] %>% as.character()
#     start <- bam$pos.GB[row]
#     end <- bam$pos_end.GB[row]
#     
#     gr <- GRanges(
#         seqnames = Rle(chr),
#         ranges = IRanges(start = start, end = end)
#     )
#     string_1 <<- bam$seq.GB[row]  #TODO Fix this...
#     string_2 <<- genome %>%
#         BSgenome::getSeq(., names = gr) %>%
#         as.character()  #TODO Fix this...
#     
#     p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
#     comparison <- paste0(
#         # paste0(chr, ":", start, "-", end, "\n"),
#         p$string_1, "\n",
#         p$result, "\n",
#         p$string_2, "\n",
#         p$n_diff, " difference(s)."
#     )
#     return(comparison)
# }


#  Set up liftOver chains -----------------------------------------------------
liftover_directory <- "../../2021-1105-1107/liftOver"

#  129S1
chain_129S1_to_mm10 <- paste0(
    liftover_directory, "/", "129S1-SvImJ-to-mm10.over.chain.munged"
) %>%
    rtracklayer::import.chain()

chain_mm10_to_129S1 <- paste0(
    liftover_directory, "/", "mm10-to-129S1-SvImJ.over.chain.munged"
) %>%
    rtracklayer::import.chain()

#  CAST
chain_CAST_to_mm10 <- paste0(
    liftover_directory, "/", "CAST-EiJ-to-mm10.over.chain.munged"
) %>%
    rtracklayer::import.chain()

chain_mm10_to_CAST <- paste0(
    liftover_directory, "/", "mm10-to-CAST-EiJ.over.chain.munged"
) %>%
    rtracklayer::import.chain()


#  Set up python code ---------------------------------------------------------
reticulate::use_condaenv(
    "/Users/kalavattam/miniconda3/envs/r41_env",
    required = TRUE
)

python_code <- "
def compare(string_1, string_2, no_match_c=' ', match_c='|'):
    if len(string_2) < len(string_1):
        string_1, string_2 = string_2, string_1
    result = ''
    n_diff = 0
    for c1, c2 in zip(string_1, string_2):
        if c1 == c2:
            result += match_c
        else:
            result += no_match_c
            n_diff += 1
    delta = len(string_2) - len(string_1)
    result += delta * no_match_c
    n_diff += delta
    return (result, n_diff)


string_1 = r.string_1
string_2 = r.string_2
result, n_diff = compare(string_1, string_2, no_match_c='*')
"

#  Test -----------------------------------------------------------------------
tbl <- uniline.a.ambiguous.dedup.129S1

int <- 50
tbl.sample <- tbl[sample(1:nrow(tbl), int, replace = FALSE), ]

tbl.sample

c.sample.even <- compareEvenBSgenome(
    tibble = tbl.sample,
    row = 1:nrow(tbl.sample),
    BSgenome = Mmusculus129S1SvImJ
)
c.sample.even %>% cat(paste0("\n\n", .))


# c.comp129S1 <- sapply(row, comp129S1, bam = bam)
# c.compNmasked <- sapply(row, compNmasked, bam = bam)
# c.compNinserted <- sapply(row, compNinserted, bam = bam)
