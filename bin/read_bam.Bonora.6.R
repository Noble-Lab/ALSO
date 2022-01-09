#!/usr/bin/env Rscript

#  Load libraries
library(BSgenome)
library(dplyr)
library(forcats)
library(ggplot2)
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
library(Rsamtools)
library(rtracklayer)
library(scales)
library(treemap)
library(viridis)

options(pillar.sigfig = 8, scipen = 10000)

setwd("/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora")


tmp.final.129S1.rds <- readRDS("tmp.final.129S1.rds")
tmp.final.CAST.rds <- readRDS("tmp.final.CAST.rds")

tmp.final.129S1.rds$intersection_GB <- factor(
    tmp.final.129S1.rds$intersection_GB,
    levels = levels(tmp.final.129S1.rds$intersection_GB)[c(9, 5, 3, 2, 12:10, 7, 6, 4, 13, 8, 1)]
)

tmp.final.CAST.rds$intersection_GB <- factor(
    tmp.final.CAST.rds$intersection_GB,
    levels = levels(tmp.final.CAST.rds$intersection_GB)[c(8, 5, 3, 2, 11, 10, 9, 7, 6, 4, 12, 1)]
)

tmp.final.129S1.rds$liftOver_reason <- factor(
    tmp.final.129S1.rds$liftOver_reason,
    levels = levels(tmp.final.129S1.rds$liftOver_reason)[c(4, 1, 2, 3)]
)

tmp.final.CAST.rds$liftOver_reason <- factor(
    tmp.final.CAST.rds$liftOver_reason,
    levels = levels(tmp.final.CAST.rds$liftOver_reason)[c(4, 2, 1, 3)]
)


# -----------------------------------------------------------------------------
# tmp.final.129S1.rds %>% ggplot(., aes(x = intersection_GB)) +
#     geom_bar(alpha = 0.5) +
#     ggtitle("KA's 129S1 assignments with respect to\nGB's assignments") +
#     geom_text(stat = "count", aes(label = scales::comma(..count.., accuracy = 1)), size = 3, vjust = -0.5) +
#     ylab("") +
#     xlab("intersections with GB's assignments") +
#     scale_y_continuous(labels = comma)
# 
# tmp.final.CAST.rds %>% ggplot(., aes(x = intersection_GB)) +
#     geom_bar(alpha = 0.5) +
#     ggtitle("KA's CAST assignments with respect to\nGB's assignments") +
#     geom_text(stat = "count", aes(label = scales::comma(..count.., accuracy = 1)), size = 3, vjust = -0.5) +
#     ylab("") +
#     xlab("intersections with GB's assignments") +
#     scale_y_continuous(labels = comma)

# ----------------------
# tmp.final.129S1.rds$intersection_GB %>% table()
# tmp.final.CAST.rds$intersection_GB %>% table()

# df.129S1 <- tmp.final.129S1.rds$intersection_GB %>% as.data.frame() %>% dplyr::count(., tmp.final.129S1.rds$intersection_GB)
# group <- df.129S1$`tmp.final.129S1.rds$intersection_GB`
# value <- df.129S1$n
# data <- data.frame(group, value)
# treemap::treemap(
#     data,
#     index = "group",
#     vSize = "value",
#     type = "index",
#     title = "KA's 129S1 assignments with respect to\nGB's assignments",
#     palette = viridis(data$group %>% n_distinct()),
# )
# 
# df.CAST <- tmp.final.CAST.rds$intersection_GB %>% as.data.frame() %>% dplyr::count(., tmp.final.CAST.rds$intersection_GB)
# group <- df.CAST$`tmp.final.CAST.rds$intersection_GB`
# value <- df.CAST$n
# data <- data.frame(group, value)
# treemap::treemap(
#     data,
#     index = "group",
#     vSize = "value",
#     type = "index",
#     title = "KA's CAST assignments with respect to\nGB's assignments",
#     palette = viridis(data$group %>% n_distinct()),
# )

# ----------------------
tmp.final.129S1.rds %>% ggplot(., aes(x = assignment, color = intersection_GB, fill = intersection_GB)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ggtitle("KA's 129S1 assignments with respect to\nGB's assignments") +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

tmp.final.CAST.rds %>% ggplot(., aes(x = assignment, color = intersection_GB, fill = intersection_GB)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ggtitle("KA's CAST assignments with respect to\nGB's assignments") +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)



# -----------------------------------------------------------------------------
# tmp.final.129S1.rds %>% ggplot(., aes(x = liftOver_reason)) +
#     geom_bar(alpha = 0.5) +
#     ggtitle("KA's 129S1 assignments with respect to\nGB's assignments") +
#     geom_text(stat = "count", aes(label = scales::comma(..count.., accuracy = 1)), size = 3, vjust = -0.5) +
#     ylab("") +
#     xlab("liftOver status: CAST to mm10") +
#     theme(axis.text.x=element_text(angle = 30,hjust = 1)) +
#     scale_y_continuous(labels = comma)
# 
# tmp.final.CAST.rds %>% ggplot(., aes(x = liftOver_reason)) +
#     geom_bar(alpha = 0.5) +
#     ggtitle("KA's CAST assignments with respect to\nGB's assignments") +
#     geom_text(stat = "count", aes(label = scales::comma(..count.., accuracy = 1)), size = 3, vjust = -0.5) +
#     ylab("") +
#     xlab("liftOver status: CAST to mm10") +
#     theme(axis.text.x=element_text(angle = 30,hjust = 1)) +
#     scale_y_continuous(labels = comma)

# ----------------------
tmp.final.129S1.rds %>% ggplot(., aes(x = assignment, color = liftOver_reason, fill = liftOver_reason)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ggtitle("KA's 129S1 assignments") +
    ylab("") +
    xlab("liftOver status: 129S1 to mm10") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

tmp.final.CAST.rds %>% ggplot(., aes(x = assignment, color = liftOver_reason, fill = liftOver_reason)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ggtitle("KA's CAST assignments") +
    ylab("") +
    xlab("liftOver status: CAST to mm10") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)


# -----------------------------------------------------------------------------
# tmp.final.129S1.rds %>% ggplot(., aes(x = equal_pos_GB_liftOver)) +
#     geom_bar(alpha = 0.5) +
#     ggtitle("KA's 129S1 assignments with respect to\nGB's assignments") +
#     geom_text(stat = "count", aes(label = scales::comma(..count.., accuracy = 1)), size = 3, vjust = -0.5) +
#     ylab("") +
#     xlab("same genomic position following liftOver") +
#     scale_y_continuous(labels = comma)
# 
# tmp.final.CAST.rds %>% ggplot(., aes(x = equal_pos_GB_liftOver)) +
#     geom_bar(alpha = 0.5) +
#     ggtitle("KA's CAST assignments with respect to\nGB's assignments") +
#     geom_text(stat = "count", aes(label = scales::comma(..count.., accuracy = 1)), size = 3, vjust = -0.5) +
#     ylab("") +
#     xlab("same genomic position following liftOver") +
#     scale_y_continuous(labels = comma)

# ----------------------
tmp.final.129S1.rds %>% ggplot(., aes(x = assignment, color = equal_pos_GB_liftOver, fill = equal_pos_GB_liftOver)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ggtitle("KA's 129S1 assignments") +
    ylab("") +
    xlab("same genomic position as GB's assignments (following liftOver)") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

tmp.final.CAST.rds %>% ggplot(., aes(x = assignment, color = equal_pos_GB_liftOver, fill = equal_pos_GB_liftOver)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ggtitle("KA's CAST assignments") +
    ylab("") +
    xlab("same genomic position as GB's assignments (following liftOver)") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)


# -----------------------------------------------------------------------------
#  Set up liftOver chains
liftover_directory <- "../../2021-1105-1107/liftOver"

# --
#  129S1
chain_129S1_to_mm10 <- paste0(liftover_directory, "/", "129S1-SvImJ-to-mm10.over.chain.munged") %>%
    rtracklayer::import.chain()
chain_mm10_to_129S1 <- paste0(liftover_directory, "/", "mm10-to-129S1-SvImJ.over.chain.munged") %>%
    rtracklayer::import.chain()

# --
#  CAST
chain_CAST_to_mm10 <- paste0(liftover_directory, "/", "CAST-EiJ-to-mm10.over.chain.munged") %>%
    rtracklayer::import.chain()
chain_mm10_to_CAST <- paste0(liftover_directory, "/", "mm10-to-CAST-EiJ.over.chain.munged") %>%
    rtracklayer::import.chain()


# ----------------------
#  Set up python code
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

# ----------------------
#  Set up functions
compareSequencesKA <- function(bam, row, genome) {
    chr <- bam$rname[row] %>% as.character()
    start <- bam$pos[row]
    end <- bam$pos[row] + 49
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq[row]  #TODO Fix this...
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this...
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    comparison <- paste0(
        p$string_1, "\n",
        p$result, "\n",
        p$string_2, "\n",
        p$n_diff, " difference(s)."
    )
    return(comparison)
}


compareSequencesGB <- function(bam, row, genome) {
    chr <- bam$GB_rname[row] %>% as.character()
    start <- bam$GB_pos[row]
    end <- bam$GB_pos[row] + 49
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq[row]  #TODO Fix this...
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this...
    
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


# -----------------------------------------------------------------------------
# --
#  129S1
bam <- tmp.final.129S1.rds
equal_pos_GB_liftOver <- bam %>%
    tibble::rowid_to_column() %>%
    filter(., (equal_pos_GB_liftOver == FALSE)) %>%
    dplyr::rename(id = rowid)
row <- equal_pos_GB_liftOver$id %>% as.double()

comparison.GB.CAST_129S1_Nmasked <- sapply(row, compareSequencesGB, bam = bam, genome = MmusculusCAST129Nmasked)
comparison.GB.CAST_129S1_Inserted <- sapply(row, compareSequencesGB, bam = bam, genome = MmusculusCAST129Inserted)
comparison.KA <- sapply(row, compareSequencesKA, bam = bam, genome = Mmusculus129S1SvImJ)

WN.129S1 <- tibble::tibble(
    comparison.GB.CAST_129S1_Nmasked,
    comparison.GB.CAST_129S1_Inserted,
    comparison.KA
)

saveRDS(WN.129S1, file = "WN.129S1.rds", compress = FALSE)

# --
#  CAST
bam <- tmp.final.CAST.rds
equal_pos_GB_liftOver <- bam %>%
    tibble::rowid_to_column() %>%
    filter(., (equal_pos_GB_liftOver == FALSE)) %>%
    dplyr::rename(id = rowid)
row <- equal_pos_GB_liftOver$id %>% as.double()

comparison.GB.CAST_129S1_Nmasked <- sapply(row, compareSequencesGB, bam = bam, genome = MmusculusCAST129Nmasked)
comparison.GB.CAST_129S1_Inserted <- sapply(row, compareSequencesGB, bam = bam, genome = MmusculusCAST129Inserted)
comparison.KA <- sapply(row, compareSequencesKA, bam = bam, genome = MmusculusCASTEiJ)


WN.CAST <- tibble::tibble(
    comparison.GB.CAST_129S1_Nmasked,
    comparison.GB.CAST_129S1_Inserted,
    comparison.KA
)

saveRDS(WN.CAST, file = "WN.CAST.rds", compress = FALSE)


# -----------------------------------------------------------------------------
# --
#  129S1
set.seed(24)
sample <- sample(1:nrow(WN.129S1), 20, replace = FALSE)

paste0(
    "\n\n# ----------------------------------------------------\n",
    "mm10 assembly, N-masked at sites of 129S1-SvImJ and CAST-EiJ SNPs;",
    "\ncoordinates selected by GB's approach\n",
    "# ----------------------\n",
    WN.129S1$comparison.GB.CAST_129S1_Nmasked[sample],
    "\n\nmm10 assembly, 129S1-SvImJ and CAST-EiJ SNPs inserted;",
    "\ncoordinates selected by GB's approach\n",
    "# ----------------------\n",
    WN.129S1$comparison.GB.CAST_129S1_Inserted[sample],
    "\n\n129S1-SvImJ assembly;",
    "\ncoordinates selected by KA's approach\n",
    "# ----------------------\n",
    WN.129S1$comparison.KA[sample],
    "\n\n"
) %>% 
    cat()

# --
#  CAST
set.seed(24)
sample <- sample(1:nrow(WN.CAST), 20, replace = FALSE)

paste0(
    "\n\n# ----------------------------------------------------\n",
    "mm10 assembly, N-masked at sites of 129S1-SvImJ and CAST-EiJ SNPs;",
    "\ncoordinates selected by GB's approach\n",
    "# ----------------------\n",
    WN.CAST$comparison.GB.CAST_129S1_Nmasked[sample],
    "\n\nmm10 assembly, 129S1-SvImJ and CAST-EiJ SNPs inserted;",
    "\ncoordinates selected by GB's approach\n",
    "# ----------------------\n",
    WN.CAST$comparison.GB.CAST_129S1_Inserted[sample],
    "\n\nCAST-EiJ assembly;",
    "\ncoordinates selected by KA's approach\n",
    "# ----------------------\n",
    WN.CAST$comparison.KA[sample],
    "\n\n"
) %>% 
    cat()
