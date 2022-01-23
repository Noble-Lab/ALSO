#!/usr/bin/env Rscript

library(BSgenome)
library(caret)
library(ggplot2)
library(Mmusculus129S1SvImJ)
library(MmusculusCASTEiJ)
library(Mmusculus129Inserted)
library(MmusculusCASTInserted)
library(MmusculusCAST129Inserted)
library(Mmusculus129Nmasked)
library(MmusculusCASTNmasked)
library(MmusculusCAST129Nmasked)
library(magrittr)
library(pheatmap)
library(Rsamtools)
library(scales)
library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)


#  Set up functions -----------------------------------------------------------
`%notin%` <- Negate(`%in%`)


compGenome <- function(rname, pos, pos_end, seq, row, genome) {
    chr <- rname[row] %>% as.character()
    start <- pos[row]
    end <- pos_end[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- seq[row]
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()
    
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


compGenomeDifferences <- function(rname, pos, pos_end, seq, row, genome) {
    chr <- rname[row] %>% as.character()
    start <- pos[row]
    end <- pos_end[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- seq[row]
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()
    
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


convertPercent <- function(x) {
    if(is.numeric(x)) {
        ifelse(is.na(x), x, paste0(round(x * 100L, 2), "%"))
    } else x
}


evaluateOperation <- function(operation = operation) {
    return(eval(parse(text = operation), envir = .GlobalEnv))
}


loadRDS <- readRDS


loadTibbleFromRDS <- function(variable, file) {
    command <- paste0(variable, " <- loadRDS(file = \"", file, "\")") %>% as.list()
    return(eval(parse(text = command), envir = globalenv()))
}


makeOperation <- function(variable, command) {
    operation <- paste(variable, command) #%>% paste(., collapse = "; ")
    return(operation)
}


#  Set up work directory, etc. (locations TB∆) --------------------------------
directory_user <- "/Users/kalavattam"
directory_project <- "Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
directory_bam <- "results/kga0/2022-0121_segregated_reads.thresh_SNP_1.thresh_Q_30"
directory_data <- "data/kga0"
path.1 <- paste0(directory_user, "/", directory_project)
path.2 <- paste0(directory_bam)
path.3 <- paste0(directory_data)

rm(directory_bam, directory_data, directory_project, directory_user)

setwd(paste0(path.1, "/", path.3))


#  Script name ----------------------------------------------------------------
script <- "test.PE-processing.part-7.R"


#  Load in the KA and GB datasets ---------------------------------------------

#  Having run parts 5 and 6 first, load *.Rdata files for them
load("test.PE-processing.part-5.Rdata")  # KA
load("test.PE-processing.part-6.Rdata")  # GB

# chromosome <- "chr1"
chromosome <- "chrX"

#  Clean up unneeded variables, etc. from the environment
operation <- paste0(
    "rm(",
        variable, ", ",
        variable_b, ", ",
        variable_m, ", ",
        variable_uniline,
    ")"
)
evaluateOperation(operation)

rm(suffix, variable, variable_b, variable_m, variable_uniline)


#  Rename variables to make things clearer ------------------------------------
bi.AS.pmin <- b.full
uni.AS.pmin <- AS.pmin

bi.joint.GB <- joint.GB
uni.joint.GB <- uniline.joint.GB

rm(AS.pmin, b.full, joint.GB, uniline.joint.GB)


#  Create factor levels for NA entries ----------------------------------------
uni.AS.pmin$assignment <- uni.AS.pmin$assignment %>% as.character()
uni.AS.pmin$assignment[is.na(uni.AS.pmin$assignment)] <- "NA"
uni.AS.pmin$assignment <- uni.AS.pmin$assignment %>% forcats::as_factor()
uni.AS.pmin$assignment %>% levels()

uni.AS.pmin$assignment <- uni.AS.pmin$assignment %>%
    plyr::revalue(c("Neutral" = "Ambiguous"))


#  Deconcatenate the uni.AS.pmin tibble ----------------------------------------
df <- uni.AS.pmin %>%
    dplyr::select(
        -trinary, -trinary.r,
        -assignment_trinary, -assignment_trinary.r
    )

df <- df %>% dplyr::mutate(
    AS.mm10.odd = AS.mm10,
    AS.129S1.odd = AS.129S1,
    AS.CAST.odd = AS.CAST,
    assignment.odd = assignment,
    AS.mm10.even = AS.mm10,
    AS.129S1.even = AS.129S1,
    AS.CAST.even = AS.CAST,
    assignment.even = assignment
) %>%
    dplyr::select(-AS.mm10, -AS.129S1, -AS.CAST, -assignment) %>%
    dplyr::relocate(
        c(
            AS.mm10.odd, AS.129S1.odd, AS.CAST.odd, assignment.odd,
            AS.mm10.even, AS.129S1.even, AS.CAST.even, assignment.even
        ),
        .before = coordinate.odd
    )

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
        "dplyr::select(-tibble)"
)
operation <- makeOperation(paste0("bi.AS.pmin"), command)
evaluateOperation(operation)

uni.AS.pmin <- df

# #  Again, check to make sure that there are no more than two entries per
# #+ group ID
# n_occur <- data.frame(table(eval(parse(text = paste0("AS.pmin")))$groupid))
# print(n_occur[n_occur$Freq > 2, ])

#  Clean up
rm(df, n_occur, even, odd)


###############################################################################
names(uni.joint.GB)[names(uni.joint.GB) == "coordinate.GB.even"] <-
    "coordinate.even"
names(uni.joint.GB)[names(uni.joint.GB) == "coordinate.GB.odd"] <-
    "coordinate.odd"

joint <- dplyr::full_join(
    uni.AS.pmin, uni.joint.GB, by = "coordinate.odd"
)
(joint$coordinate.odd %notin% uni.AS.pmin$coordinate.odd) %>% table()
(joint$coordinate.odd %notin% uni.joint.GB$coordinate.odd) %>% table()
(joint$coordinate.even.x == joint$coordinate.even.y) %>% table()

names(joint)[names(joint) == "coordinate.even.x"] <- "coordinate.KA.even"
names(joint)[names(joint) == "coordinate.even.y"] <- "coordinate.GB.even"
names(joint)[names(joint) == "coordinate.odd"] <- "coordinate.KA.odd"

names(joint)[names(joint) == "assignment.odd"] <- "assignment.KA.odd"
names(joint)[names(joint) == "assignment.even"] <- "assignment.KA.even"

joint$coordinate.GB.odd <- joint$coordinate.KA.odd
joint <- joint %>%
    dplyr::relocate(coordinate.GB.odd, .after = qmpos.GB.odd) %>%
    dplyr::select(-discrepancy.GB.odd, -discrepancy.GB.even)

joint$assignment.KA.odd <- paste0("KA.", joint$assignment.KA.odd)
joint$assignment.GB.odd[is.na(joint$assignment.GB.odd)] <- "GB.NA"
joint$assignment.GB.even[is.na(joint$assignment.GB.even)] <- "NA"

joint$assignment.KA.even <- joint$assignment.KA.even %>% as.character()

order <- c("129", "CAST", "Ambiguous", "NA")
joint$assignment.KA.even <- joint$assignment.KA.even %>%
    stringr::str_remove("-EiJ") %>%
    stringr::str_remove("S1-SvImJ") %>%
    factor(., levels = order)
joint$assignment.KA.odd <- joint$assignment.KA.odd %>%
    stringr::str_remove("-EiJ") %>%
    stringr::str_remove("S1-SvImJ") %>%
    factor(., levels = paste0("KA.", order))

order <- c("129", "CAST", "Ambiguous", "Contra", "NA")
joint$assignment.GB.even <- joint$assignment.GB.even %>%
    stringr::str_remove("GB.") %>%
    stringr::str_remove("S1") %>%
    factor(., levels = order)
joint$assignment.GB.odd <- joint$assignment.GB.odd %>%
    stringr::str_remove("S1") %>%
    factor(., levels = paste0("GB.", order))

# #  Check
# joint$assignment.KA.odd %>% head()
# joint$assignment.KA.odd %>% tail()
# joint$assignment.KA.odd %>% levels()
# 
# joint$assignment.KA.even %>% head()
# joint$assignment.KA.even %>% tail()
# joint$assignment.KA.even %>% levels()
# 
# joint$assignment.GB.odd %>% head()
# joint$assignment.GB.odd %>% tail()
# joint$assignment.GB.odd %>% levels()
# 
# joint$assignment.GB.even %>% head()
# joint$assignment.GB.even %>% tail()
# joint$assignment.GB.even %>% levels()

order <- c(
    "KA.129 × GB.Ambiguous",
    "KA.129 × GB.129",
    "KA.129 × GB.CAST",
    "KA.129 × GB.Contra",
    "KA.129 × GB.NA",
    "KA.CAST × GB.129",
    "KA.CAST × GB.CAST",
    "KA.CAST × GB.Ambiguous",
    "KA.CAST × GB.NA",
    "KA.Ambiguous × GB.129",
    "KA.Ambiguous × GB.CAST",
    "KA.Ambiguous × GB.Ambiguous",
    "KA.Ambiguous × GB.Contra",
    "KA.Ambiguous × GB.NA",
    "KA.NA × GB.129",
    "KA.NA × GB.CAST",
    "KA.NA × GB.Ambiguous",
    "KA.NA × GB.Contra"
)
joint$KA_by_GB <- paste0(
    joint$assignment.KA.odd, " × ", joint$assignment.GB.odd
) %>%
    factor(levels = order)

order <- c(
    "GB.129 × KA.129",
    "GB.129 × KA.CAST",
    "GB.129 × KA.Ambiguous",
    "GB.129 × KA.NA",
    "GB.CAST × KA.129",
    "GB.CAST × KA.CAST",
    "GB.CAST × KA.Ambiguous",
    "GB.CAST × KA.NA",
    "GB.Ambiguous × KA.129",
    "GB.Ambiguous × KA.CAST",
    "GB.Ambiguous × KA.Ambiguous",
    "GB.Ambiguous × KA.NA",
    "GB.Contra × KA.129",
    "GB.Contra × KA.Ambiguous",
    "GB.Contra × KA.NA",
    "GB.NA × KA.129",
    "GB.NA × KA.CAST",
    "GB.NA × KA.Ambiguous"
)
joint$GB_by_KA <- paste0(
    joint$assignment.GB.odd, " × ", joint$assignment.KA.odd
) %>%
    factor(levels = order)

# #  Check
# joint$KA_by_GB %>% head()
# joint$GB_by_KA %>% head()
# joint$KA_by_GB %>% levels()
# joint$GB_by_KA %>% levels()

joint <- joint %>%
    dplyr::relocate(c(KA_by_GB, GB_by_KA), .before = AS.mm10.odd)


#  Generate confusion matrix --------------------------------------------------

#  'data' is for predicted data, i.e., the data from my pipeline; 'reference'
#+ is for true results, i.e., the data from GB's pipeline
example <- confusionMatrix(
    data = joint$assignment.KA.even,
    reference = joint$assignment.GB.even,
    dnn = c("KA", "GB")
)
example.table <- example$table %>%
    as.data.frame() %>%
    tibble::as_tibble()
example.table <- example.table %>%
    mutate(KA = factor(KA, rev(unique(KA))))

ggplot(example.table, aes(x = GB, y = KA, fill = Freq)) +
    geom_tile(color = "white") +
    theme_bw() +
    coord_equal() +
    scale_fill_distiller(palette = "Spectral", direction = -1) +
    guides(fill = FALSE) +
    labs(title = "Value distribution") +
    geom_text(aes(label = Freq), color = "black")

rm(example, example.table)


#  Set up liftOver chains -----------------------------------------------------
liftover_directory <- "liftOver"

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

joint.assign <- joint %>%
    dplyr::group_by(KA_by_GB) %>%
    dplyr::group_split()


#  ----------------------------------------------------------------------------

#  [[5]] KA.129 × GB.NA
bam <- joint.assign[[5]] %>% dplyr::select(sort(dplyr::current_vars()))
bam <- bam[
    sample(
        1:nrow(bam),
        if (nrow(bam) < 500) {
            nrow(bam)
        } else {
            500
        },
        replace = FALSE
    ), 
]  # If commented out, then process all observations in tibble
bam <- bam %>%
    tibble::rowid_to_column() %>%
    dplyr::rename(ID.post_sample = rowid)

row <- bam$ID.post_sample %>% as.double()


#  Comparisons ----------------------------------------------------------------
test.129S1.even <- sapply(
    row, compGenome,
    rname = bam$rname.129S1.even,
    pos = bam$pos.129S1.even,
    pos_end = bam$pos_end.129S1.even,
    seq = bam$seq.129S1.even,
    genome = Mmusculus129S1SvImJ
)
test.129S1.odd <- sapply(
    row, compGenome,
    rname = bam$rname.129S1.odd,
    pos = bam$pos.129S1.odd,
    pos_end = bam$pos_end.129S1.odd,
    seq = bam$seq.129S1.odd,
    genome = Mmusculus129S1SvImJ
)

test.129S1.even %>% paste0("\n\n", .) %>% cat()
test.129S1.odd %>% paste0("\n\n", .) %>% cat()


test.CAST_129_Nmasked.even <- sapply(
    row, compGenome,
    rname = bam$rname.GB.even,
    pos = bam$pos.GB.even,
    pos_end = bam$pos_end.GB.even,
    seq = bam$seq.GB.even,
    genome = MmusculusCAST129Nmasked
)
test.CAST_129_Nmasked.odd <- sapply(
    row, compGenome,
    rname = bam$rname.GB.odd,
    pos = bam$pos.GB.odd,
    pos_end = bam$pos_end.GB.odd,
    seq = bam$seq.GB.odd,
    genome = MmusculusCAST129Nmasked
)

test.CAST_129_Nmasked.even %>% paste0("\n\n", .) %>% cat()
test.CAST_129_Nmasked.odd %>% paste0("\n\n", .) %>% cat()


test.CAST_129_inserted.even <- sapply(
    row, compGenome,
    rname = bam$rname.GB.even,
    pos = bam$pos.GB.even,
    pos_end = bam$pos_end.GB.even,
    seq = bam$seq.GB.even,
    genome = MmusculusCAST129Inserted
)
test.CAST_129_inserted.odd <- sapply(
    row, compGenome,
    rname = bam$rname.GB.odd,
    pos = bam$pos.GB.odd,
    pos_end = bam$pos_end.GB.odd,
    seq = bam$seq.GB.odd,
    genome = MmusculusCAST129Inserted
)


test.diff.129.even <- sapply(
    row, compGenomeDifferences,
    rname = bam$rname.129S1.even,
    pos = bam$pos.129S1.even,
    pos_end = bam$pos_end.129S1.even,
    seq = bam$seq.129S1.even,
    genome = Mmusculus129S1SvImJ
)
test.diff.129.odd <- sapply(
    row, compGenomeDifferences,
    rname = bam$rname.129S1.odd,
    pos = bam$pos.129S1.odd,
    pos_end = bam$pos_end.129S1.odd,
    seq = bam$seq.129S1.odd,
    genome = Mmusculus129S1SvImJ
)

hist(
    test.diff.129.even,
    breaks = max(test.diff.129.even),
    main = "129S1/SvImJ, even",
    xlab = "mismatches per seq",
    right = FALSE
)
hist(
    test.diff.129.odd,
    breaks = max(test.diff.129.odd),
    main = "129S1/SvImJ, odd",
    xlab = "mismatches per seq",
    right = FALSE
)

#TODO 1/2 Get the 'seq' associated with the pmin(), i.e., the worse alignment;
#TODO 2/2 put them in a separate column; then use them for the above analyses

#TODO 1/3 Test seq.even and seq.odd for most mismatches; select the one with;
#TODO 2/3 the most mismatches; put it in its own column; then use them for the
#TODO 3/3 above analyses
