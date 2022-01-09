#!/usr/bin/env Rscript

#  Load libraries
library(BSgenome)
library(dplyr)
library(forcats)
library(GenomicRanges)
library(magrittr)
library(Mmusculus129S1SvImJ)
library(MmusculusCASTEiJ)
library(Mmusculus129Inserted)
library(MmusculusCASTInserted)
library(MmusculusCAST129Inserted)
library(Mmusculus129Nmasked)
library(MmusculusCASTNmasked)
library(MmusculusCAST129Nmasked)
library(parallel)
library(purrr)
library(reticulate)
library(Rsamtools)

options(pillar.sigfig = 8, scipen = 10000)


# -----------------------------------------------------------------------------
"/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora" %>% setwd()


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
files <- list.files(pattern = "\\.dedup.MarkDuplicates.sort.chr1.bam$")
n <- files %>% length()
variables <- paste0("bam", 1:n)
variables_rname <- paste0("bam", 1:n, "$rname")
command_pipe <- paste0("<- bam", 1:n, " %>% ")
chromosomes <- c(paste0("chr", 1:19), "chrX")

#  Read in .bam filenames
mapply(
    assign, variables, files, MoreArgs = list(envir = parent.frame())
)

#  Assign .bam information as list
dedup.129 <- bam1 %>% Rsamtools::scanBam(., param = map_params)
dedup.CAST <- bam2 %>% Rsamtools::scanBam(., param = map_params)
dedup.Nmasked <- bam3 %>% Rsamtools::scanBam(., param = map_params)
dedup.mm10 <- bam4 %>% Rsamtools::scanBam(., param = map_params)

#  Convert .bam information from list to dataframe to tibble
variables <- c(
    "dedup.129",
    "dedup.CAST",
    "dedup.Nmasked",
    "dedup.mm10"
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


#  Check the heads of the tibbles ---------------------------------------------
dedup.129 %>% head()
dedup.CAST %>% head()
dedup.Nmasked %>% head()
dedup.mm10 %>% head()



# -----------------------------------------------------------------------------
liftover_directory <- "../../2021-1105-1107/liftOver"

chain_129_to_mm10 <- paste0(liftover_directory, "/", "129S1-SvImJ-to-mm10.over.chain.munged") %>%
    rtracklayer::import.chain()
chain_mm10_to_129 <- paste0(liftover_directory, "/", "mm10-to-129S1-SvImJ.over.chain.munged") %>%
    rtracklayer::import.chain()

chain_CAST_to_mm10 <- paste0(liftover_directory, "/", "CAST-EiJ-to-mm10.over.chain.munged") %>%
    rtracklayer::import.chain()
chain_mm10_to_CAST <- paste0(liftover_directory, "/", "mm10-to-CAST-EiJ.over.chain.munged") %>%
    rtracklayer::import.chain()


# -----------------------------------------------------------------------------
reticulate::use_condaenv(
    "/Users/kalavattam/miniconda3/envs/r41_env",
    required = TRUE
)

python_code <- 
"
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


compareSequences <- function(bam, row, genome) {
    chr <- bam$rname[row] %>% as.character()
    start <- bam$pos[row]
    end <- bam$pos[row] + 49
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    gr_liftOver_129_to_mm10 <- liftOver(
        gr,
        chain_129_to_mm10
    ) %>% unlist()
    gr_liftOver_mm10_to_129 <- liftOver(
        gr,
        chain_mm10_to_129
    ) %>% unlist()
    
    string_1 <<- bam$seq[row]
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    comparison <- paste0(
        p$string_1, "\n",
        p$result, "\n",
        p$string_2, "\n",
        p$n_diff, " difference(s)."
    )
    return(comparison)
}


compareSequences129 <- function(bam, row) {
    chr <- bam$rname[row] %>% as.character()
    start <- bam$pos[row]
    end <- bam$pos[row] + 49
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    gr_liftOver_129_to_mm10 <- liftOver(
        gr,
        chain_129_to_mm10
    ) %>% unlist()
    gr_liftOver_mm10_to_129 <- liftOver(
        gr,
        chain_mm10_to_129
    ) %>% unlist()
    
    string_1 <<- bam$seq[row]
    string_2 <<- Mmusculus129S1SvImJ %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    comparison <- paste0(
        p$string_1, "\n",
        p$result, "\n",
        p$string_2, "\n",
        p$n_diff, " difference(s)."
    )
    return(comparison)
}


compareSequencesCAST <- function(bam, row) {
    chr <- bam$rname[row] %>% as.character()
    start <- bam$pos[row]
    end <- bam$pos[row] + 49
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    gr_liftOver_CAST_to_mm10 <- liftOver(
        gr,
        chain_CAST_to_mm10
    ) %>% unlist()
    gr_liftOver_mm10_to_CAST <- liftOver(
        gr,
        chain_mm10_to_CAST
    ) %>% unlist()
    
    string_1 <<- bam$seq[row]
    string_2 <<- MmusculusCASTEiJ %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    comparison <- paste0(
        p$string_1, "\n",
        p$result, "\n",
        p$string_2, "\n",
        p$n_diff, " difference(s)."
    )
    return(comparison)
}


# https://stackoverflow.com/questions/31050556/parallel-version-of-sapply
mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
    FUN <- match.fun(FUN)
    answer <- parallel::mclapply(X = X, FUN = FUN, ...)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) 
        names(answer) <- X
    if (!isFALSE(simplify) && length(answer)) 
        simplify2array(answer, higher = (simplify == "array"))
    else answer
}


printList <- function(list) {
    for (item in 1:length(list)) {
        cat(head(list[[item]]), "\n\n")
    }
}

# -----------------
row <- c(100:110)
bam <- dedup.CAST

out129 <- list()
out129 <- sapply(row, compareSequences129, bam = bam)
# out129 <- mcsapply(row, compareSequences129, bam = bam, mc.cores = 4)  # Freezes...

outCAST <- list()
outCAST <- sapply(row, compareSequencesCAST, bam = bam)
# outCAST <- mcsapply(row, compareSequencesCAST, bam = bam, mc.cores = 4)

out129[9] %>% cat()
outCAST[9] %>% cat()

# lapply(out129, cat)
# lapply(outCAST, cat)

out129 %>% printList()
outCAST %>% printList()

# -----------------
row <- c(100:110)
bam <- dedup.129

out129 <- list()
out129 <- sapply(row, compareSequences129, bam = bam)
# out129 <- mcsapply(row, compareSequences129, bam = bam, mc.cores = 4)  # Freezes...

outCAST <- list()
outCAST <- sapply(row, compareSequencesCAST, bam = bam)
# outCAST <- mcsapply(row, compareSequencesCAST, bam = bam, mc.cores = 4)

out129[9] %>% cat()
outCAST[9] %>% cat()

# lapply(out129, cat)
# lapply(outCAST, cat)

out129 %>% printList()
outCAST %>% printList()

# -----------------
row <- c(100:110)
bam <- dedup.129

dedup.129.out129 <- list()
dedup.129.out129 <- sapply(row, compareSequences, bam = bam, genome = Mmusculus129S1SvImJ)

dedup.129. <- list()
dedup.129.outCAST <- sapply(row, compareSequences, bam = bam, genome = MmusculusCASTEiJ)

dedup.129.out129Inserted <- list()
dedup.129.out129Inserted <- sapply(row, compareSequences, bam = bam, genome = Mmusculus129Inserted)

dedup.129.out129 %>% printList()
dedup.129.outCAST %>% printList()
dedup.129.out129Inserted %>% printList()


row <- c(100:110)
bam <- dedup.CAST

dedup.CAST.out129 <- list()
dedup.CAST.out129 <- sapply(row, compareSequences, bam = bam, genome = Mmusculus129S1SvImJ)

dedup.CAST.outCAST <- list()
dedup.CAST.outCAST <- sapply(row, compareSequences, bam = bam, genome = MmusculusCASTEiJ)

dedup.CASTout129Inserted <- list()
dedup.CASTout129Inserted <- sapply(row, compareSequences, bam = bam, genome = Mmusculus129Inserted)

dedup.CASTout129 %>% printList()
dedup.CASToutCAST %>% printList()
dedup.CASTout129Inserted %>% printList()


dedup.129.out129[9] %>% cat()
dedup.CAST.out129[9] %>% cat()

dedup.129.outCAST[9] %>% cat()
dedup.CAST.outCAST[9] %>% cat()
