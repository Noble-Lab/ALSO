
library(BSgenome)
library(Mmusculus129S1SvImJ)
library(MmusculusCASTEiJ)
library(Mmusculus129Inserted)
library(MmusculusCASTInserted)
library(MmusculusCAST129Inserted)
library(Mmusculus129Nmasked)
library(MmusculusCASTNmasked)
library(MmusculusCAST129Nmasked)
library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)


#  Set up work directory (locations TB∆) --------------------------------------
directory_user <- "/Users/kalavattam"
directory_base <- "Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do"
directory_work <- "get_unique_fragments/Bonora"

setwd(
    paste(directory_user, directory_base, directory_work, sep = "/")
)
rm(directory_user, directory_base, directory_work)


#  Set up functions -----------------------------------------------------------
`%notin%` <- Negate(`%in%`)


evaluateOperation <- function(operation = operation) {
    return(eval(parse(text = operation), envir = .GlobalEnv))
}


loadRDS <- readRDS


loadTibbleFromRDS <- function(variable, file) {
    command <- paste0(variable, " <- loadRDS(file = \"", file, "\")") %>% as.list()
    return(eval(parse(text = command), envir = globalenv()))
}


makeOperation <- function(variable = variable, command = command) {
    operation <- paste(variable, command) #%>% paste(., collapse = "; ")
    return(operation)
}


mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
    # https://stackoverflow.com/questions/31050556/parallel-version-of-sapply
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


#  Set up functions for liftOver conversions, etc. ----------------------------
comp129S1 <- function(bam, row, genome = Mmusculus129S1SvImJ) {
    chr <- bam$rname.129S1[row] %>% as.character()
    start <- bam$pos.129S1[row]
    end <- bam$pos_end.129S1[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq.129S1[row]  #TODO Fix this
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this
    
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


comp129S1_lift <- function(bam, row, genome) {
    chr <- bam$rname.129S1[row] %>% as.character()
    start <- bam$pos.129S1[row]
    end <- bam$pos_end.129S1[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    gr_liftOver_129S1_to_mm10 <- liftOver(
        gr,
        chain_129S1_to_mm10
    ) %>%
        unlist()
    string_1 <<- bam$seq.129S1[row]  #TODO Fix this
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr_liftOver_129S1_to_mm10) %>%
        as.character()  #TODO Fix this
    
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


compNmasked <- function(bam, row, genome = MmusculusCAST129Nmasked) {
    chr <- bam$rname.GB[row] %>% as.character()
    start <- bam$pos.GB[row]
    end <- bam$pos_end.GB[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq.GB[row]  #TODO Fix this
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this
    
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


compNinserted <- function(bam, row, genome = MmusculusCAST129Inserted) {
    chr <- bam$rname.GB[row] %>% as.character()
    start <- bam$pos.GB[row]
    end <- bam$pos_end.GB[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq.GB[row]  #TODO Fix this
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this
    
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


diff129S1 <- function(bam, row, genome = Mmusculus129S1SvImJ) {
    chr <- bam$rname.129S1[row] %>% as.character()
    start <- bam$pos.129S1[row]
    end <- bam$pos_end.129S1[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq.129S1[row]  #TODO Fix this
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    difference <- paste0(
        # paste0(chr, ":", start, "-", end, "\n"),
        # p$string_1, "\n",
        # p$result, "\n",
        # p$string_2, "\n",
        # p$n_diff, " difference(s)."
        p$n_diff
    ) %>%
        as.numeric()
    return(difference)
}


diff129S1_lift <- function(bam, row, genome) {
    chr <- bam$rname.129S1[row] %>% as.character()
    start <- bam$pos.129S1[row]
    end <- bam$pos_end.129S1[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    gr_liftOver_129S1_to_mm10 <- liftOver(
        gr,
        chain_129S1_to_mm10
    ) %>%
        unlist()
    string_1 <<- bam$seq.129S1[row]  #TODO Fix this
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr_liftOver_129S1_to_mm10) %>%
        as.character()  #TODO Fix this
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    difference <- paste0(
        # paste0(chr, ":", start, "-", end, "\n"),
        # p$string_1, "\n",
        # p$result, "\n",
        # p$string_2, "\n",
        # p$n_diff, " difference(s)."
        p$n_diff
    ) %>%
        as.numeric()
    return(difference)
}


diffNmasked <- function(bam, row, genome = MmusculusCAST129Nmasked) {
    chr <- bam$rname.GB[row] %>% as.character()
    start <- bam$pos.GB[row]
    end <- bam$pos_end.GB[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq.GB[row]  #TODO Fix this
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    difference <- paste0(
        # paste0(chr, ":", start, "-", end, "\n"),
        # p$string_1, "\n",
        # p$result, "\n",
        # p$string_2, "\n",
        # p$n_diff, " difference(s)."
        p$n_diff
    ) %>%
        as.numeric()
    return(difference)
}


diffNinserted <- function(bam, row, genome = MmusculusCAST129Inserted) {
    chr <- bam$rname.GB[row] %>% as.character()
    start <- bam$pos.GB[row]
    end <- bam$pos_end.GB[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq.GB[row]  #TODO Fix this
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    difference <- paste0(
        # paste0(chr, ":", start, "-", end, "\n"),
        # p$string_1, "\n",
        # p$result, "\n",
        # p$string_2, "\n",
        # p$n_diff, " difference(s)."
        p$n_diff
    ) %>%
        as.numeric()
    return(difference)
}


#  If you were testing functions, then remove the variables
rm(chr, start, end, gr, gr_liftOver_129S1_to_mm10)


# Set up liftOver chains ------------------------------------------------------
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

rm(liftover_directory)


#  Set up python code
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


#  Load KA_GB.joint tibble from .rds file -------------------------------------
chromosome <- "chrX"
j <- loadRDS(file = paste0("KA_GB.joint.", chromosome, ".rds"))

j <- j %>% dplyr::rename(
    c(AS.GB = tag.AS.GB, GB_assignment = GB_assignment.2)
    ) %>%  # new_name = old_name
    dplyr::select(sort(dplyr::current_vars())) %>%
    dplyr::select(
        -in_alt_CAST,
        -in_ref_129S1,
        -in_ambig,
        -in_contra,
        -GB_assignment.1
    )
#TODO 1/2 Add at least the rename() to the preceding script,
#TODO 2/2 *.5.test-scraps.*.R; add perhaps the select() too

#TODO 1/2 Problem: mm10-N-masked values are most NAs, except for AS.mm10: Go
#TODO 2/2 back (to the preceding script, I think) to check on why/how this is


#  (Name of section TBD) ------------------------------------------------------
j.assign <- j %>%
    dplyr::relocate(
        c(GB_assignment, KA_assignment), .before = AS.129S1
    ) %>%
    dplyr::group_by(GB_KA_assignment) %>%
    dplyr::group_split()
j.assign

#  Indices beginning with KA.Ambiguous: [[11:15]]
#+ - 11: GB.Ambiguous
#+ - 12: GB.129S1
#+ - 13: GB.CAST
#+ - 14: GB.Not_present
#+ - 15: GB.Contra

#  Indices beginning with GB.Ambiguous: [[3, 8, 11, 18]]
#+ - 3: KA.129S1
#+ - 8: KA.CAST
#+ - 11: KA.Ambiguous
#+ - 18: KA.NA

# #  KA.Ambiguous × GB.129S1
# bam <- j.assign[[12]] %>%
#     dplyr::select(sort(dplyr::current_vars()))  # KA.Ambiguous by GB.129S1

#  GB.Ambiguous × KA.129S1
bam <- j.assign[[3]] %>%
    dplyr::select(sort(dplyr::current_vars()))  # GB.Ambiguous by KA.129S1

bam <- bam %>%
    dplyr::arrange(rname.liftOver_129S1_mm10, pos.liftOver_129S1_mm10) %>%
    dplyr::rename(c(AS.Nmasked = AS.mm10)) %>%  # new_name = old_name
    dplyr::select(-tidyselect::contains(c("CAST", ".mm10"))) %>%
    tibble::rowid_to_column() %>%
    dplyr::rename(ID = rowid)  # new_name = old_name

#  Back-ups for testing...
# bam.bak <- bam
bam <- bam.bak

#  If below line is commented out, then process all observations in tibble
bam <- bam[sample(1:nrow(bam), 500, replace = FALSE), ]  
bam <- bam %>%
    dplyr::arrange(rname.liftOver_129S1_mm10, pos.liftOver_129S1_mm10) %>%
    tibble::rowid_to_column() %>%
    dplyr::rename(ID.post_sample = rowid)

row <- bam$ID.post_sample %>% as.double()


#  Test comparisons -----------------------------------------------------------
c.comp129S1 <- sapply(
    row, comp129S1, bam = bam
)

c.compNmasked <- sapply(
    row, compNmasked, bam = bam
)
c.compNinserted <- sapply(
    row, compNinserted, bam = bam
)

# c.comp129S1_Nmasked <- sapply(
#     row, comp129S1_lift, bam = bam, genome = MmusculusCAST129Nmasked
# )
# c.comp129S1_Ninserted <- sapply(
#     row, comp129S1_lift, bam = bam, genome = MmusculusCAST129Inserted
# )


# -----------------
#  c.comp129S1 vs. c.compNmasked vs. c.compNinserted
#  Direct comparisonss
c.comp129S1 %>% paste0("\n\n", .) %>% cat()
c.compNmasked %>% paste0("\n\n", .) %>% cat()
c.compNinserted %>% paste0("\n\n", .) %>% cat()

#  Number of differences (vs. *Nmasked)
c.comp129S1[c.comp129S1 != c.compNmasked] %>% paste0("\n\n", .) %>% cat()
c.compNmasked[c.comp129S1 != c.compNmasked] %>% paste0("\n\n", .) %>% cat()

#  Number of differences (vs. *N-inserted)
c.comp129S1[c.comp129S1 != c.compNinserted] %>% paste0("\n\n", .) %>% cat()
c.compNinserted[c.comp129S1 != c.compNinserted] %>% paste0("\n\n", .) %>% cat()

#  Number of differences (w/r/t *Nmasked)
c.comp129S1[c.comp129S1 != c.compNmasked] %>% paste0("\n\n", .) %>% cat()
c.compNmasked[c.comp129S1 != c.compNmasked] %>% paste0("\n\n", .) %>% cat()
c.compNinserted[c.comp129S1 != c.compNmasked] %>% paste0("\n\n", .) %>% cat()


# # -----------------
# #  comp129S1 vs. comp129S1_lift
# #  Direct comparisonss
# c.comp129S1 %>% paste0("\n\n", .) %>% cat()
# c.comp129S1_Nmasked %>% paste0("\n\n", .) %>% cat()
# c.comp129S1_Ninserted %>% paste0("\n\n", .) %>% cat()
# 
# #  Number of differences (vs. *_Nmasked)
# c.comp129S1[c.comp129S1 != c.comp129S1_Nmasked] %>%
#     paste0("\n\n", .) %>%
#     cat()
# c.comp129S1_Nmasked[c.comp129S1 != c.comp129S1_Nmasked] %>%
#     paste0("\n\n", .) %>%
#     cat()
# 
# #  Number of differences (vs. *Ninserted)
# c.comp129S1[c.comp129S1 != c.compNinserted] %>%
#     paste0("\n\n", .) %>%
#     cat()
# c.compNinserted[c.comp129S1 != c.compNinserted] %>%
#     paste0("\n\n", .) %>%
#     cat()
# 
# #  Number of differences (vs. *_Ninserted)
# c.comp129S1[c.comp129S1 != c.comp129S1_Ninserted] %>%
#     paste0("\n\n", .) %>%
#     cat()
# c.comp129S1_Ninserted[c.comp129S1 != c.comp129S1_Ninserted] %>%
#     paste0("\n\n", .) %>%
#     cat()
# 
# #  Number of differences (*Ninserted vs. *_Ninserted)
# c.compNinserted[c.compNinserted != c.comp129S1_Ninserted] %>%
#     paste0("\n\n", .) %>%
#     cat()  # GB's method
# c.comp129S1_Ninserted[c.compNinserted != c.comp129S1_Ninserted] %>%
#     paste0("\n\n", .) %>%
#     cat()  # KA's method


paste0(
    "\n\n# ----------------------------------------------------\n",
    "mm10 assembly, N-masked at sites of 129S1-SvImJ and CAST-EiJ SNPs;",
    "\ncoordinates selected by GB's approach\n",
    "# ----------------------\n",
    c.compNmasked[c.comp129S1 != c.compNmasked],
    "\n\nmm10 assembly, 129S1-SvImJ and CAST-EiJ SNPs inserted;",
    "\ncoordinates selected by GB's approach\n",
    "# ----------------------\n",
    c.compNinserted[c.comp129S1 != c.compNmasked],
    "\n\n129S1-SvImJ assembly;",
    "\ncoordinates selected by KA's approach\n",
    "# ----------------------\n",
    c.comp129S1[c.comp129S1 != c.compNmasked],
    "\n\n"
) %>% 
    cat()

paste0(
    "\n\n# ----------------------------------------------------\n",
    "mm10 assembly, N-masked at sites of 129S1-SvImJ and CAST-EiJ SNPs;",
    "\ncoordinates selected by GB's approach\n",
    "# ----------------------\n",
    c.compNmasked[20],
    "\n\nmm10 assembly, 129S1-SvImJ and CAST-EiJ SNPs inserted;",
    "\ncoordinates selected by GB's approach\n",
    "# ----------------------\n",
    c.compNinserted[20],
    "\n\n129S1-SvImJ assembly;",
    "\ncoordinates selected by KA's approach\n",
    "# ----------------------\n",
    c.comp129S1[20],
    "\n\n"
) %>% 
    cat()


#  Test differences -----------------------------------------------------------
d.diff129S1 <- sapply(
    row, diff129S1, bam = bam
)
d.diffNmasked <- sapply(
    row, diffNmasked, bam = bam
)
d.diffNinserted <- sapply(
    row, diffNinserted, bam = bam
)

# suffix <- "_test_2021-1201.rds"
suffix <- "_test_2022-0109.rds"
saveRDS(d.diff129S1, file = paste0("d.diff129S1", suffix))
saveRDS(d.diffNmasked, file = paste0("d.diffNmasked", suffix))
saveRDS(d.diffNinserted, file = paste0("d.diffNinserted", suffix))

hist(
    d.diff129S1,
    breaks = 27,
    xlab = "differences",
    right = FALSE
)
hist(
    d.diffNmasked,
    breaks = 27,
    xlab = "differences",
    right = FALSE
)
hist(
    d.diffNinserted,
    breaks = 27,
    xlab = "differences",
    right = FALSE
)

# d.diff129S1_Nmasked <- sapply(
#     row, diff129S1_lift, bam = bam, genome = MmusculusCAST129Nmasked
# )
# d.diff129S1_Ninserted <- sapply(
#     row, diff129S1_lift, bam = bam, genome = MmusculusCAST129Inserted
# )


#  Script is completed; clean up environment ----------------------------------
rm(list = ls())
