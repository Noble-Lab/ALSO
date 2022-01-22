#!/usr/bin/env Rscript

library(BSgenome)
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

setwd(paste0(path.1))


#  Script name ----------------------------------------------------------------
script <- "test.PE-processing.2021-1221.2022-0117.part-5.R"


#  Set up what to query from .bam files ---------------------------------------
map_info <- c(
    "qname", "flag", "rname", "pos", "mapq", "cigar",
    "mrnm", "mpos", "isize", "seq", "qual"
)
map_params <- Rsamtools::ScanBamParam(what = map_info)


#  Munge the KA-generated GB .bam files ---------------------------------------
setwd(paste0(path.1, "/", path.2))

# chromosome <- "chr1"
chromosome <- "chrX"

file <- list.files(pattern = "\\.bam")
variable_GB <- c(
    "GB.alt.CAST",
    "GB.ambig",
    "GB.contra",
    "GB.ref.129S1"
)

#  Check that the correct files will be assigned to the correct variables
mapply(
    assign, variable_GB, file, MoreArgs = list(envir = parent.frame())
)


#  Assign .bam information as list --------------------------------------------
command <- paste0(
    "<- ", variable_GB, " %>% ",
        "Rsamtools::scanBam(., param = map_params)"
)
operation <- makeOperation(variable_GB, command)
evaluateOperation(operation)

rm(map_info, map_params)


#  Convert .bam information from list to dataframe to tibble ------------------
command <- paste0(
    "<- ", variable_GB, " %>% ",
        "as.data.frame() %>% ",
        "as_tibble()"
)
operation <- makeOperation(variable_GB, command)
evaluateOperation(operation)


#  Reorder rname factor levels ------------------------------------------------
chromosomes <- c(paste0("chr", c(1:19)), "chrX", "chrY", "chrM")
command <- paste0(
    "<- forcats::fct_relevel(", variable_GB, "$rname, chromosomes)"
)
operation <- makeOperation(paste0(variable_GB, "$rname"), command)
evaluateOperation(operation)


#  Drop unused rname factor levels --------------------------------------------
command <- paste0(
    "<- ", variable_GB, "$rname %>% forcats::fct_drop()"
)
operation <- makeOperation(paste0(variable_GB, "$rname"), command)
evaluateOperation(operation)


#  Drop rows that are not chr1-19, chrX ---------------------------------------
chromosomes <- c(paste0("chr", 1:19), "chrX")
command <- paste0(
    "<- ", variable_GB, " %>% ", "filter(., rname %in% chromosomes)"
)
operation <- makeOperation(variable_GB, command)
evaluateOperation(operation)


#  Create and append pos_end, mpos_end columns --------------------------------
command <- paste0("<- ", variable_GB, "$pos + 49")
operation <- makeOperation(paste0(variable_GB, "$pos_end"), command)
evaluateOperation(operation)

command <- paste0("<- ", variable_GB, "$mpos + 49")
operation <- makeOperation(paste0(variable_GB, "$mpos_end"), command)
evaluateOperation(operation)

command <- paste0(
    "<- ", variable_GB, " %>% ",
        "dplyr::relocate(pos_end, .after = pos) %>% ",
        "dplyr::relocate(mpos_end, .after = mpos)"
)
operation <- makeOperation(variable_GB, command)
evaluateOperation(operation)


#  Create a "temporary tag" for use in distinguishing entries
command <- paste0(
    "<- ", variable_GB, " %>% ",
        "tidyr::unite(",
            "coordinate, ",
            "c(\"qname\", \"rname\", \"pos\"), ",
            # "c(\"qname\", \"rname\", \"pos\", \"pos_end\"), ",
            "sep = \"_\", ",
            "remove = FALSE",
        ") %>% ",
            "dplyr::relocate(coordinate, .before = qname)"
)
operation <- makeOperation(variable_GB, command)
evaluateOperation(operation)


#  Assignments and then smart join of GB tibbles ------------------------------

#  Explicitly label the assignments from Giancarlo's pipeline, then join the
#+ related tibbles together

#NOTE Remember, "alt" is "CAST", "ref" is "129"

GB.alt.CAST$assignment <- "GB.CAST"
GB.ref.129S1$assignment <- "GB.129S1"
GB.ambig$assignment <- "GB.Ambiguous"
GB.contra$assignment <- "GB.Contra"

joint.GB <- dplyr::bind_rows(GB.alt.CAST, GB.ref.129S1, GB.ambig, GB.contra)
colnames(joint.GB) <- paste0(colnames(joint.GB), ".GB")
colnames(joint.GB)[1] <- "coordinate"
colnames(joint.GB)[colnames(joint.GB) %>% length()] <- "assignment_GB"


##############################################################################
GB.alt.CAST <- GB.alt.CAST %>% dplyr::relocate(mpos, .after = pos)
GB.ambig <- GB.ambig %>% dplyr::relocate(mpos, .after = pos)
GB.contra <- GB.contra %>% dplyr::relocate(mpos, .after = pos)
GB.ref.129S1 <- GB.ref.129S1 %>% dplyr::relocate(mpos, .after = pos)


##############################################################################


#  Load in the "KA assignments" -----------------------------------------------
setwd(paste0(path.1, "/", path.3))

#  Having run part 5 first, load part-5 .Rdata into environment
load("test.PE-processing.2021-1221.2022-0117.part-5.Rdata")

# chromosome <- "chr1"
chromosome <- "chrX"


#  Create factor levels for NA entries ----------------------------------------
AS.pmin$assignment <- AS.pmin$assignment %>% as.character()
AS.pmin$assignment[is.na(AS.pmin$assignment)] <- "NA"
AS.pmin$assignment <- AS.pmin$assignment %>% forcats::as_factor()
AS.pmin$assignment %>% levels()

AS.pmin$assignment <- AS.pmin$assignment %>% plyr::revalue(
    c("Neutral" = "Ambiguous")
)

#  Deconcatenate the AS.pmin tibble -------------------------------------------
df <- AS.pmin %>%
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
    dplyr::select(-AS.mm10, -AS.129S1, -AS.CAST, -assignment)

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
operation <- makeOperation(paste0("AS.pmin"), command)
evaluateOperation(operation)

# #  Again, check to make sure that there are no more than two entries per
# #+ group ID
# n_occur <- data.frame(table(eval(parse(text = paste0("AS.pmin")))$groupid))
# print(n_occur[n_occur$Freq > 2, ])

#  Clean up
rm(df, n_occur, even, odd)

AS.pmin.split <- AS.pmin %>%
    dplyr::group_by(assignment) %>%
    dplyr::group_split()
# [[1]]  NA
# A tibble: 13,856 x 6
# [[2]] Neutral
# A tibble: 139,840 x 6
# [[3]] 129S1-SvImJ
# A tibble: 50,460 x 6
# [[4]] CAST-EiJ
# A tibble: 47,726 x 6

variable_KA <- c(
    "KA.NA",
    "KA.Ambiguous",
    "KA.129S1",
    "KA.CAST"
)
KA.NA <- AS.pmin.split[[1]]
KA.Ambiguous <- AS.pmin.split[[2]]
KA.129S1 <- AS.pmin.split[[3]]
KA.CAST <- AS.pmin.split[[4]]


#  Join the KA and GB data sets -----------------------------------------------

#  Strip the final field from GB.*$coordinate (to match KA.*$coordinate)
GB.alt.CAST$coordinate <- gsub('(.*)_\\w+', '\\1', GB.alt.CAST$coordinate)
GB.ambig$coordinate <- gsub('(.*)_\\w+', '\\1', GB.ambig$coordinate)
GB.contra$coordinate <- gsub('(.*)_\\w+', '\\1', GB.contra$coordinate)
GB.ref.129S1$coordinate <- gsub('(.*)_\\w+', '\\1', GB.ref.129S1$coordinate)
joint.GB$coordinate <- gsub('(.*)_\\w+', '\\1', joint.GB$coordinate)

variable_joint <- variable_KA %>%
    strsplit(., "\\.") %>%
    lapply(., `[[`, 2) %>%
    unlist() %>%
    paste0("KA_GB.", .)
command <- paste0(
    "<- full_join(",
        variable_KA, ", ",
        "joint.GB, by = \"coordinate\"",
    ") %>% ",
        "tidyr::drop_na(\"assignment\")"
)
operation <- makeOperation(variable_joint, command)
evaluateOperation(operation)

operation <- paste0(
    variable_joint, "$assignment_GB[",
        "is.na(", variable_joint, "$assignment_GB)",
    "] <- \"GB.not_present\""
)
evaluateOperation(operation)

command <- paste0(
    "<- ", variable_joint, " %>% ",
        "dplyr::rename(., GB_assignment.1 = assignment_GB)"
)
operation <- makeOperation(variable_joint, command)
evaluateOperation(operation)

#  Clean up
operation <- paste0("rm(", variable_KA, ")")
evaluateOperation(operation)


#  Identify intersections between "assignment" and "dedup.joint" --------------
for (i in 1:length(variable_joint)) {
    command <- paste0(
        "<- ifelse(",
        variable_joint[i], "$coordinate %in% ", variable_GB[4], "$coordinate, ",
        "'1', ",
        "'0'",
        ")"
    )
    operation <- makeOperation(
        paste0(variable_joint[i], "$in_ref_129S1"), command
    )
    print(operation)
    evaluateOperation(operation)
    
    command <- paste0(
        "<- ifelse(",
        variable_joint[i], "$coordinate %in% ", variable_GB[1], "$coordinate, ",
        "'1', ",
        "'0'",
        ")"
    )
    operation <- makeOperation(
        paste0(variable_joint[i], "$in_alt_CAST"), command
    )
    print(operation)
    evaluateOperation(operation)
    
    command <- paste0(
        "<- ifelse(",
        variable_joint[i], "$coordinate %in% ", variable_GB[2], "$coordinate, ",
        "'1', ",
        "'0'",
        ")"
    )
    operation <- makeOperation(
        paste0(variable_joint[i], "$in_ambig"), command
    )
    print(operation)
    evaluateOperation(operation)
    
    command <- paste0(
        "<- ifelse(",
        variable_joint[i], "$coordinate %in% ", variable_GB[3], "$coordinate, ",
        "'1', ",
        "'0'",
        ")"
    )
    operation <- makeOperation(
        paste0(variable_joint[i], "$in_contra"), command
    )
    print(operation)
    evaluateOperation(operation)
    
    command <- paste0(
        "<- paste0(",
        variable_joint[i], "$in_ref_129S1, ",
        variable_joint[i], "$in_alt_CAST, ",
        variable_joint[i], "$in_ambig, ",
        variable_joint[i], "$in_contra",
        ")"
    )
    operation <- makeOperation(
        paste0(variable_joint[i], "$GB_intersection"), command
    )
    print(operation)
    evaluateOperation(operation)
    
    command <- paste0(
        "<- ", variable_joint[i], "$GB_intersection", " %>% ",
        "as.factor()"
    )
    operation <- makeOperation(
        paste0(variable_joint[i], "$GB_intersection"), command
    )
    print(operation)
    evaluateOperation(operation)
}

#  Reorder the factor levels for column/variable GB_intersection
order_intersection <- c("1000", "0100", "0010", "0001", "0000")
command <- paste0(
    "<- ", "forcats::fct_relevel(",
        variable_joint, "$GB_intersection, order_intersection",
    ")"
)
operation <- makeOperation(
    paste0(variable_joint, "$GB_intersection"), command
)
evaluateOperation(operation)

#  Copy the vector of binary names to a new variable; then, give them easier-
#+ to-read names 
command <- paste0("<- ", variable_joint, "$GB_intersection")
operation <- makeOperation(
    paste0(variable_joint, "$GB_assignment.2"), command
)
evaluateOperation(operation)

command <- paste0(
    "<- ", variable_joint, "$GB_assignment.2 %>% plyr::revalue(c(",
        "'0000' = 'GB.Not_present', ",
        "'0001' = 'GB.Contra', ",
        "'0010' = 'GB.Ambiguous', ",
        "'0100' = 'GB.CAST', ",
        "'1000' = 'GB.129S1'",
    "))"
)
operation <- makeOperation(
    paste0(variable_joint, "$GB_assignment.2"), command
)
evaluateOperation(operation)

command <- paste0(
    "<- ", variable_joint, "$assignment %>% plyr::revalue(c(",
        "'129S1-SvImJ' = 'KA.129S1', ",
        "'CAST-EiJ' = 'KA.CAST', ",
        "'Neutral' = 'KA.Ambiguous', ",
        "'NA' = 'KA.NA'",
    "))"
)
operation <- makeOperation(paste0(variable_joint, "$KA_assignment"), command)
evaluateOperation(operation)  # The errors reported here are fine...


#  Set up KA × GB assignment factors, KA × GB assignment factors --------------
KA_GB.129S1$KA_GB_assignment <- paste0(
    KA_GB.129S1$KA_assignment, ' × ', KA_GB.129S1$GB_assignment.2
)
KA_GB.CAST$KA_GB_assignment <- paste0(
    KA_GB.CAST$KA_assignment, ' × ', KA_GB.CAST$GB_assignment.2
)
KA_GB.NA$KA_GB_assignment <- paste0(
    KA_GB.NA$KA_assignment, ' × ', KA_GB.NA$GB_assignment.2
)
KA_GB.Ambiguous$KA_GB_assignment <- paste0(
    KA_GB.Ambiguous$KA_assignment, ' × ', KA_GB.Ambiguous$GB_assignment.2
)

KA_GB.129S1$GB_KA_assignment <- paste0(
    KA_GB.129S1$GB_assignment.2, ' × ', KA_GB.129S1$KA_assignment
)
KA_GB.CAST$GB_KA_assignment <- paste0(
    KA_GB.CAST$GB_assignment.2, ' × ', KA_GB.CAST$KA_assignment
)
KA_GB.NA$GB_KA_assignment <- paste0(
    KA_GB.NA$GB_assignment.2, ' × ', KA_GB.NA$KA_assignment
)
KA_GB.Ambiguous$GB_KA_assignment <- paste0(
    KA_GB.Ambiguous$GB_assignment.2, ' × ', KA_GB.Ambiguous$KA_assignment
)

command <- paste0("<- ", variable_joint, "$KA_GB_assignment %>% as_factor()")
operation <- makeOperation(
    paste0(variable_joint, "$KA_GB_assignment"), command
)
evaluateOperation(operation)

command <- paste0("<- ", variable_joint, "$GB_KA_assignment %>% as_factor()")
operation <- makeOperation(
    paste0(variable_joint, "$GB_KA_assignment"), command
)
evaluateOperation(operation)


#  Generate table -------------------------------------------------------------
command <- paste0(
    "<- table(", variable_joint, "$GB_assignment.2) %>% ",
        "as.data.frame() %>% ",
        "as_tibble() %>%",
        "dplyr::rename(c(GB.assignment = Var1, ", variable_joint, " = Freq))"
)
operation <- makeOperation(paste0(variable_joint, ".table"), command)
evaluateOperation(operation)

KA_GB.all.table <- full_join(
    KA_GB.129S1.table, KA_GB.CAST.table, by = "GB.assignment"
) %>%
    full_join(., KA_GB.Ambiguous.table, by = "GB.assignment") %>%
    full_join(., KA_GB.NA.table, by = "GB.assignment") %>%
    purrr::map_df(rev)


#  Generate table of proportions ----------------------------------------------
command <- paste0(
    "<- ", variable_joint, "$GB_assignment.2 %>% ",
        "table() %>% ",
        "prop.table() %>% ",
        "as.data.frame() %>% ",
        "as_tibble() %>% ",
        "dplyr::rename(c(GB.assignment = ., ", variable_joint, " = Freq))"
)
operation <- makeOperation(paste0(variable_joint, ".table.prop"), command)
evaluateOperation(operation)

table.prop <- full_join(
    KA_GB.129S1.table.prop, KA_GB.CAST.table.prop, by = "GB.assignment"
) %>%
    full_join(., KA_GB.Ambiguous.table.prop, by = "GB.assignment") %>%
    full_join(., KA_GB.NA.table.prop, by = "GB.assignment") %>%
    purrr::map_df(rev)


#  Generate table of percentages ----------------------------------------------
command <- paste0("<- ", variable_joint, ".table.prop")
operation <- makeOperation(paste0(variable_joint, ".table.percent"), command)
evaluateOperation(operation)

command <- paste0(
    "<- ", variable_joint, ".table.percent$", variable_joint, " %>% ",
        "convertPercent()"
)
operation <- makeOperation(
    paste0(variable_joint, ".table.percent$", variable_joint), command
)
evaluateOperation(operation)

table.percent <- full_join(
    KA_GB.129S1.table.percent, KA_GB.CAST.table.percent, by = "GB.assignment"
) %>%
    full_join(., KA_GB.Ambiguous.table.percent, by = "GB.assignment") %>%
    full_join(., KA_GB.NA.table.percent, by = "GB.assignment") %>%
    purrr::map_df(rev)


#  Create a full table of counts, including row and column sums ---------------
tmp.1 <- KA_GB.all.table %>% summarise(across(2:5, ~sum(.)))

tmp.2 <- KA_GB.all.table %>%
    select(colnames(.)[2:length(colnames(.))]) %>% 
    rowSums() %>%
    as_tibble_col()

tmp.3 <- tmp.2 %>%
    summarise(across(1, ~sum(.)))

tmp.4 <- bind_rows(tmp.2, tmp.3)

table.full <- bind_rows(KA_GB.all.table, tmp.1)
table.full <- bind_cols(table.full, tmp.4)
table.full$GB.assignment <- table.full$GB.assignment %>% as.character()
table.full[6, 1] <- "GB.sum"
table.full$GB.assignment <- table.full$GB.assignment %>% as_factor()
table.full <- table.full %>% dplyr::rename(., sum.row = value)
colnames(table.full)[6] <- "KA.sum"

rm(tmp.1, tmp.2, tmp.3, tmp.4)


#  Create a full table of proportions -----------------------------------------
tmp.1 <- table.full$KA.sum[-6] %>% as_tibble_col()
tmp.2 <- tmp.1 / colSums(tmp.1)

table.prop.full <- cbind(table.prop, tmp.2) %>% as_tibble()
colnames(table.prop.full)[6] <- "KA.sum"

rm(tmp.1, tmp.2)


#  Create a full table of percentages -----------------------------------------
column <- c(
    "KA_GB.129S1",
    "KA_GB.CAST",
    "KA_GB.Ambiguous",
    "KA_GB.NA",
    "KA.sum"
)

table.percent.full <- list()
for (i in column) {
    command <- paste0("<- table.prop.full$", i," %>% convertPercent()")
    operation <- makeOperation(paste0("table.percent.full$", i), command)
    print(operation)
    evaluateOperation(operation)
}
table.percent.full <- table.percent.full %>%
    as_tibble() %>%
    purrr::map_df(rev)

table.percent.full <- cbind(table.percent[, 1], table.percent.full) %>%
    as_tibble()

table.percent.full$GB.assignment <- table.percent.full$GB.assignment %>%
    as.factor()


#  Print the full tables ------------------------------------------------------
table.full
table.full.df <- data.frame(
    table.full[, 2:5], row.names = table.full$GB.assignment
)
table.full.df <- table.full.df[-6, ]
colnames(table.full.df) <- c("KA.129S1", "KA.CAST", "KA.Ambiguous", "KA.NA")
table.full.df <- table.full.df[order(nrow(table.full.df):1), ]

breaks <- seq(0, 50000, 100)
pheatmap::pheatmap(
    table.full.df,
    main = paste0("KA × GB assignments, ", chromosome),
    display_numbers = TRUE,
    number_format = "%.0f",
    color = colorRampPalette(c('white','red'))(length(breaks)),
    breaks = breaks,
    border_color = "#FFFFFF",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_number = 11
)

table.prop.full
table.prop.full.df <- data.frame(
    table.prop.full[, 2:5], row.names = table.prop.full$GB.assignment
)
colnames(table.prop.full.df) <- variable_KA
table.prop.full.df <- table.prop.full.df[order(nrow(table.prop.full.df):1), ]
pheatmap::pheatmap(
    table.prop.full.df %>% as_tibble(),
    main = paste0("KA × GB assignments, ", chromosome),
    labels_row = table.prop.full$GB.assignment,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = colorRampPalette(c('white','red'))(100),
    border_color = "#FFFFFF",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_number = 11
)


#  Set up functions for liftOver conversions, etc. ----------------------------
comp129S1 <- function(bam, row, genome = Mmusculus129S1SvImJ) {
    chr <- bam$rname.129S1[row] %>% as.character()
    start <- bam$pos.129S1[row]
    end <- bam$pos_end.129S1[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq.129S1[row]  #TODO Fix this...
    string_2 <<- genome %>%
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
    string_1 <<- bam$seq.129S1[row]  #TODO Fix this...
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr_liftOver_129S1_to_mm10) %>%
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


compCAST <- function(bam, row, genome = MmusculusCASTEiJ) {
    chr <- bam$rname.CAST[row] %>% as.character()
    start <- bam$pos.CAST[row]
    end <- bam$pos_end.CAST[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq.CAST[row]  #TODO Fix this...
    string_2 <<- genome %>%
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


compCAST_lift <- function(bam, row, genome) {
    chr <- bam$rname.CAST[row] %>% as.character()
    start <- bam$pos.CAST[row]
    end <- bam$pos_end.CAST[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    gr_liftOver_CAST_to_mm10 <- liftOver(
        gr,
        chain_CAST_to_mm10
    ) %>%
        unlist()
    string_1 <<- bam$seq.CAST[row]  #TODO Fix this...
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr_liftOver_CAST_to_mm10) %>%
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


compNmasked <- function(bam, row, genome = MmusculusCAST129Nmasked) {
    chr <- bam$rname.GB[row] %>% as.character()
    start <- bam$pos.GB[row]
    end <- bam$pos_end.GB[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq.GB[row]  #TODO Fix this...
    string_2 <<- genome %>%
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


compNinserted <- function(bam, row, genome = MmusculusCAST129Inserted) {
    chr <- bam$rname.GB[row] %>% as.character()
    start <- bam$pos.GB[row]
    end <- bam$pos_end.GB[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq.GB[row]  #TODO Fix this...
    string_2 <<- genome %>%
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


diff129S1 <- function(bam, row, genome = Mmusculus129S1SvImJ) {
    chr <- bam$rname.129S1[row] %>% as.character()
    start <- bam$pos.129S1[row]
    end <- bam$pos_end.129S1[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq.129S1[row]  #TODO Fix this...
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this...
    
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
    string_1 <<- bam$seq.129S1[row]  #TODO Fix this...
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr_liftOver_129S1_to_mm10) %>%
        as.character()  #TODO Fix this...
    
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


diffCAST <- function(bam, row, genome = MmusculusCASTEiJ) {
    chr <- bam$rname.CAST[row] %>% as.character()
    start <- bam$pos.CAST[row]
    end <- bam$pos_end.CAST[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq.CAST[row]  #TODO Fix this...
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this...
    
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


diffCAST_lift <- function(bam, row, genome) {
    chr <- bam$rname.CAST[row] %>% as.character()
    start <- bam$pos.CAST[row]
    end <- bam$pos_end.CAST[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    gr_liftOver_CAST_to_mm10 <- liftOver(
        gr,
        chain_CAST_to_mm10
    ) %>%
        unlist()
    string_1 <<- bam$seq.CAST[row]  #TODO Fix this...
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr_liftOver_CAST_to_mm10) %>%
        as.character()  #TODO Fix this...
    
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
    string_1 <<- bam$seq.GB[row]  #TODO Fix this...
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this...
    
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
    string_1 <<- bam$seq.GB[row]  #TODO Fix this...
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this...
    
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


#  Join the KA_GB.*, then munge the join --------------------------------------
j <- dplyr::bind_rows(KA_GB.129S1, KA_GB.CAST, KA_GB.Ambiguous, KA_GB.NA)

#  Check on the number of GB.* entries in the join
j$GB_assignment.1 %>% table()
# .
#       GB.129S1   GB.Ambiguous        GB.CAST      GB.Contra GB.not_present 
#          32680         171106          30817             20          17259
#  It checks out...

j$GB_assignment.2 %>% table()
# .
#       GB.129S1        GB.CAST   GB.Ambiguous      GB.Contra GB.Not_present 
#          32680          30817         171106             20          17259
#  It checks out...

j$KA_assignment <- j$KA_assignment %>%
    plyr::revalue(
        c("Ambiguous" = "KA.Ambiguous")
    )
j$KA_assignment %>% table()
# .
#        KA.NA KA.Ambiguous     KA.129S1      KA.CAST 
#        13856       139840        50460        47726

j <- j %>% dplyr::rename(
        c(GB_assignment = GB_assignment.2)
    ) %>%  # new_name = old_name
    dplyr::select(sort(dplyr::current_vars())) %>%
    dplyr::select(
        -in_alt_CAST,
        -in_ref_129S1,
        -in_ambig,
        -in_contra,
        -GB_assignment.1
    )
#TODO current_vars() is deprecated

j.assign <- j %>%
    dplyr::relocate(
        c(GB_assignment, KA_assignment), .before = AS.129S1
    ) %>%
    dplyr::group_by(GB_KA_assignment) %>%
    dplyr::group_split()
j.assign

#  Indices beginning with KA.Ambiguous: [[11:15]]
#+ - 11: KA.Ambiguous × GB.Ambiguous
#+ - 12: KA.Ambiguous × GB.129S1
#+ - 13: KA.Ambiguous × GB.Not_present
#+ - 14: KA.Ambiguous × GB.Contra
#+ - 15: KA.Ambiguous × GB.CAST

#  Indices beginning with GB.Ambiguous: [[3, 8, 11, 18]]
#+ - 01: GB.Ambiguous × KA.129S1
#+ - 07: GB.Ambiguous × KA.CAST
#+ - 11: GB.Ambiguous × KA.Ambiguous
#+ - 16: GB.Ambiguous × KA.NA


#  ----------------------------------------------------------------------------
#  12. KA.Ambiguous × GB.129S1
bam <- j.assign[[12]] %>% dplyr::select(sort(dplyr::current_vars()))

#  03. GB.Ambiguous × KA.129S1
# bam <- j.assign[[3]] %>% dplyr::select(sort(dplyr::current_vars()))

#  02. GB.129S1 × KA.129S1
# bam <- j.assign[[2]] %>% dplyr::select(sort(dplyr::current_vars()))

bam <- bam %>%
    # dplyr::arrange(rname.liftOver_129S1_mm10, pos.liftOver_129S1_mm10) %>%
    dplyr::rename(c(AS.Nmasked = AS.mm10)) %>%  # new_name = old_name
    dplyr::select(-tidyselect::contains(c("CAST", ".mm10"))) %>%
    tibble::rowid_to_column() %>%
    dplyr::rename(ID = rowid)  # new_name = old_name

#  Back-ups for testing...
# bam.bak <- bam
bam <- bam.bak

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
    # dplyr::arrange(rname.liftOver_129S1_mm10, pos.liftOver_129S1_mm10) %>%
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

