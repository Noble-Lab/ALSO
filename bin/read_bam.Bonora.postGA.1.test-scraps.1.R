#!/usr/bin/env Rscript

library(Rsamtools)
library(tidyverse)

# -----------------------------------------------------------------------------
#  ...from read_bam.Bonora.4.R

#  Temporary: Set up work directory (location TB∆)
"/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora/segregatedReads.SNPTHRESH1.Q30" %>% setwd()

#  Files are from
#+ /net/noble/vol2/home/gbonora/proj/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouseDiff_Nmasked/segregatedReads.SNPTHRESH1.Q30
#+ 
#+ For more information, see the following script:
#+ /net/noble/vol2/home/gbonora/proj/2019_sciATAC_analysis/results/gbonora/20191105_sciATAC_mouseDiff_Nmasked/20191105_sciATAC_mouseDiff_Nmasked_allelicSegregation_workflow.sh


# -----------------------------------------------------------------------------
#  Set up functions
`%notin%` <- Negate(`%in%`)

loadRDS <- readRDS

makeOperation <- function(variable, command) {
    operation <- paste(variable, command) %>% paste(., collapse = "; ")
    return(operation)
}


# -----------------------------------------------------------------------------
#  What to query from .bam files
# map_info <- c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar", "mrnm", "mpos", "isize", "seq", "qual")
map_info <- c("qname", "flag", "rname", "strand", "pos", "mapq", "cigar", "mrnm", "mpos", "isize", "seq")
tag_info <- c("AS", "XS", "NM")
map_params <- Rsamtools::ScanBamParam(what = map_info, tag = tag_info)

#  Assign variables necessary for subsequent commands -------------------------
#+ 
#+ Remember, "alt" is "CAST", "ref" is "129"
file <- list.files(pattern = "\\.chr1.bam$")
variable <- file %>% strsplit(., "\\.") %>% lapply(., `[[`, 5) %>% unlist() %>% paste0("GB.", .)
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

#  Assign .bam information as list --------------------------------------------
command_pipe <- paste0("<- ", variable, " %>% ")
command <- paste0(command_pipe, "Rsamtools::scanBam(., param = map_params)")
operation <- makeOperation(variable, command)
eval(parse(text = operation))
rm(map_info, tag_info, map_params)

#  Convert .bam information from list to dataframe to tibble ------------------
#+ 
#+ And add a row ID column to each tibble
command <- paste0(command_pipe, "as.data.frame() %>% tibble::as_tibble() %>% tibble::rowid_to_column(\"ID\")")
operation <- makeOperation(variable, command)
eval(parse(text = operation))

#  Save tibbles as uncompressed .rds files ------------------------------------
command <- paste0("saveRDS(", variable, ", file = \"", variable, ".rds\", compress = FALSE)")
eval(parse(text = command))

#  Load tibbles from .rds files -----------------------------------------------
command <- paste0(variable, " <- loadRDS(", variable, ", file = \"", variable, ".rds\")")
eval(parse(text = command))

#  Reorder rname factor levels ------------------------------------------------
chromosome <- c(paste0("chr", 1:19), "chrX", "chrY", "chrM")
command <- paste0("<- forcats::fct_relevel(", variable, "$rname, chromosome)")
operation <- makeOperation(paste0(variable, "$rname"), command)
eval(parse(text = operation))

command <- paste0("<- forcats::fct_relevel(", variable, "$mrnm, chromosome)")
operation <- makeOperation(paste0(variable, "$mrnm"), command)
eval(parse(text = operation))

#  Drop unused rname factor levels --------------------------------------------
command <- paste0("<- ", variable, "$rname %>% forcats::fct_drop()")
operation <- makeOperation(paste0(variable, "$rname"), command)
eval(parse(text = operation))

command <- paste0("<- ", variable, "$mrnm %>% forcats::fct_drop()")
operation <- makeOperation(paste0(variable, "$mrnm"), command)
eval(parse(text = operation))

#  Drop rows that are not chr1-19, chrX, chrY, chrM ---------------------------
# command <- paste0(command_pipe, "filter(., rname %in% chromosome)")
# operation <- makeOperation(variable, command)
# eval(parse(text = operation))
# 
# command <- paste0(command_pipe, "filter(., mrnm %in% chromosome)")
# operation <- makeOperation(variable, command)
# eval(parse(text = operation))

#  Create and append pos_end, mpos_end vectors --------------------------------
#+ 
#+ Use 49 b/c expression is $pos + 50 exclusive, i.e., "[pos, 50)"
command <- paste0("<- ", variable, "$pos + 49")
operation <- makeOperation(paste0(variable, "$pos_end"), command)
eval(parse(text = operation))

command <- paste0("<- ", variable, "$mpos + 49")
operation <- makeOperation(paste0(variable, "$mpos_end"), command)
eval(parse(text = operation))

# -------
# --
#  Pre-check: Create tibbles with less rows for tests of the duplicate-removal
variable_abbrev <- c(paste0(variable, ".abbrev"))
command <- paste0("<- ", variable, "[1:10000, ]")
operation <- makeOperation(variable_abbrev, command)
eval(parse(text = operation))
# GB.alt.abbrev %>% nrow()  # [1] 10000
# GB.ambig.abbrev %>% nrow()  # [1] 10000
# GB.contra.abbrev %>% nrow()  # [1] 10000
# GB.ref.abbrev %>% nrow()  # [1] 10000

# --
#  Store untouched backups
bak.GB.alt <- GB.alt
bak.GB.ambig <- GB.ambig
bak.GB.contra <- GB.contra
bak.GB.ref <- GB.ref

#  Assign the abbreivated tibbles to the main variables
GB.alt <- GB.alt.abbrev
GB.ambig <- GB.ambig.abbrev
GB.contra <- GB.contra.abbrev
GB.ref <- GB.ref.abbrev

#  If necessary, clean up the environment
rm(
    # GB.alt.abbrev,
    # GB.ambig.abbrev,
    # GB.contra.abbrev,
    # GB.ref.abbrev,
    GB.alt.all_dups,
    GB.ambig.all_dups,
    GB.contra.all_dups,
    GB.ref.all_dups,
    GB.alt.each_dup,
    GB.ambig.each_dup,
    GB.contra.each_dup,
    GB.ref.each_dup
)

#  Reorder the pos_end and mpos_end columns; create a column of concatenated
#+ barcode_sequence strings (i.e., "reads"); and then create a column of
#+ concatenated barcode_fragment strings (i.e., "coordinates," where a single
#+ coordinate is rname, pos, and pos_end)
#+ 
#+ #NOTE
#+ - Terminology change: "b_s" is now "read"
#+ - Terminology change: "b_f" is now "coordinate"
command_pipe <- paste0("<- ", variable, " %>% ")
command <- paste0(
    command_pipe,
    "dplyr::relocate(pos_end, .after = pos) %>% dplyr::relocate(mpos_end, .after = mpos)",
    " %>% tidyr::unite(read, c(\"qname\", \"seq\"), sep = \"_\", remove = FALSE)",
    " %>% tidyr::unite(coordinate, c(\"qname\", \"rname\", \"pos\", \"pos_end\"), sep = \"_\", remove = FALSE)"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))

#  Order columns by AS, rname, and pos
# command <- paste0(  # Option A
#     command_pipe,
#     "dplyr::arrange(-tag.AS, mapq, coordinate, flag)"
# )
# operation <- makeOperation(variable, command)
# eval(parse(text = operation))

command <- paste0(  # Option B  # Go with option B
    "<- ", variable, "[order(", variable, "$coordinate, -", variable, "$tag.AS, -", variable, "$mapq), ]"
)
# command <- paste0(  # Option C
#     "<- ", variable, "[order(", variable, "$coordinate, -", variable, "$mapq, ", variable, "$tag.AS), ]"
# )
operation <- makeOperation(variable, command)
eval(parse(text = operation))

#  Filter out all rows without unique read
# command <- paste0(command_pipe, "dplyr::distinct(read, .keep_all = TRUE)")
# operation <- makeOperation(variable, command)
# eval(parse(text = operation))
# # GB.alt %>% nrow()  # [1] 1796387
# # GB.ambig %>% nrow()  # [1] 5333701
# # GB.contra %>% nrow()  # [1] 5882
# # GB.ref %>% nrow()  # [1] 1834526

#  Filter out all rows without unique coordinates -----------------------------
#+ (Numbers are for data sets comprised of only chr1)
# -------
#  Pre-check: nrow() 
# GB.alt %>% nrow()  # [1] 2599984
# GB.ambig %>% nrow()  # [1] 9013879
# GB.contra %>% nrow()  # [1] 7900
# GB.ref %>% nrow()  # [1] 2655891

# -------
# --
#  Pre-check: Collect duplicate information: Get one of each duplicated element
command <- paste0(command_pipe, "dplyr::filter(duplicated(.[[\"coordinate\"]]))")
operation <- makeOperation(paste0(variable,".each_dup"), command)
eval(parse(text = operation))
# GB.alt.each_dup %>% nrow()  # [1] 880904
# GB.ambig.each_dup %>% nrow()  # [1] 4090693
# GB.contra.each_dup %>% nrow()  # [1] 2319
# GB.ref.each_dup %>% nrow()  # [1] 901236

GB.alt.each_dup <- GB.alt.each_dup %>% select(ID, coordinate, rname, pos, pos_end, flag, tag.AS, mapq, seq)
GB.ambig.each_dup <- GB.ambig.each_dup %>% select(ID, coordinate, rname, pos, pos_end, flag, tag.AS, mapq, seq)
GB.contra.each_dup <- GB.contra.each_dup %>% select(ID, coordinate, rname, pos, pos_end, flag, tag.AS, mapq, seq)
GB.ref.each_dup <- GB.ref.each_dup %>% select(ID, coordinate, rname, pos, pos_end, flag, tag.AS, mapq, seq)

View(GB.ambig.each_dup)

# --
#  Pre-check: Collect duplicate information: Get all of the duplicated elements
command <- paste0(command_pipe, "dplyr::filter(coordinate %in% unique(.[[\"coordinate\"]][duplicated(.[[\"coordinate\"]])]))")
operation <- makeOperation(paste0(variable, ".all_dups"), command)
eval(parse(text = operation))
# GB.alt.all_dups %>% nrow()  # [1] 1344095
# GB.ambig.all_dups %>% nrow()  # [1] 6223797
# GB.contra.all_dups %>% nrow()  # [1] 3554
# GB.ref.all_dups %>% nrow()  # [1] 1374738

GB.alt.all_dups <- GB.alt.all_dups %>% select(ID, coordinate, rname, pos, pos_end, flag, tag.AS, mapq, seq)
GB.ambig.all_dups <- GB.ambig.all_dups %>% select(ID, coordinate, rname, pos, pos_end, flag, tag.AS, mapq, seq)
GB.contra.all_dups <- GB.contra.all_dups %>% select(ID, coordinate, rname, pos, pos_end, flag, tag.AS, mapq, seq)
GB.ref.all_dups <- GB.ref.all_dups %>% select(ID, coordinate, rname, pos, pos_end, flag, tag.AS, mapq, seq)

View(GB.ambig.all_dups)

# --
#  Pre-check: See what the dedupilcated "mints" look like...
command <- paste0("<- ", variable, "[!duplicated(", variable, "$coordinate), ]")
operation <- makeOperation(paste0(variable, ".dedup"), command)
eval(parse(text = operation))

GB.alt.dedup <- GB.alt.dedup %>% select(ID, coordinate, rname, pos, pos_end, flag, tag.AS, mapq, seq)
GB.ambig.dedup <- GB.ambig.dedup %>% select(ID, coordinate, rname, pos, pos_end, flag, tag.AS, mapq, seq)
GB.contra.dedup <- GB.contra.dedup %>% select(ID, coordinate, rname, pos, pos_end, flag, tag.AS, mapq, seq)
GB.ref.dedup <- GB.ref.dedup %>% select(ID, coordinate, rname, pos, pos_end, flag, tag.AS, mapq, seq)

View(GB.ambig.dedup)  # Examining, e.g., 215 and 9743 in the different GB.ambig tibbles... Yes, Option B is suitable
View(GB.ambig)

# -------
# --
#  Pre-check: Create tibbles with less rows for tests of the duplicate-removal
variable_abbrev <- c(
    paste0(variable, ".abbrev")
)
command <- paste0("<- ", variable, "[1:10000, ]")
operation <- makeOperation(variable_abbrev, command)
eval(parse(text = operation))
# GB.alt.abbrev %>% nrow()  # [1] 10000
# GB.ambig.abbrev %>% nrow()  # [1] 10000
# GB.contra.abbrev %>% nrow()  # [1] 10000
# GB.ref.abbrev %>% nrow()  # [1] 10000

# --
command_pipe <- paste0("<- ", variable_abbrev, " %>% ")
command <- paste0(command_pipe, "dplyr::filter(duplicated(.[[\"coordinate\"]]))")
operation <- makeOperation(paste0(variable,".each_dup"), command)
eval(parse(text = operation))
# GB.alt.each_dup %>% nrow()  # [1] 880904
# GB.ambig.each_dup %>% nrow()  # [1] 4090693
# GB.contra.each_dup %>% nrow()  # [1] 2319
# GB.ref.each_dup %>% nrow()  # [1] 901236

#  Test: Duplicate removal: Are we left w/deduplicated rows w/the best AS's?
command <- paste0(command_pipe, "dplyr::distinct(coordinate, .keep_all = TRUE)")
operation <- makeOperation(variable_abbrev, command)
eval(parse(text = operation))
# GB.alt.abbrev %>% nrow()  # [1] 9795
# GB.ambig.abbrev %>% nrow()  # [1] 5380
# GB.contra.abbrev %>% nrow()  # [1] 5582
# GB.ref.abbrev %>% nrow()  # [1] 9820




# GB.alt %>% nrow()  # [1] 1719080
# GB.ambig %>% nrow()  # [1] 4923186
# GB.contra %>% nrow()  # [1] 5581
# GB.ref %>% nrow()  # [1] 1754655

#  Remember, "alt" is "CAST", "ref" is "129"
GB.alt.CAST <- GB.alt
GB.ref.129S1 <- GB.ref

#  Clean up
#TODO Determine better spots fot this...
rm(GB.alt, GB.ref)
rm(chromosome, file, suffix, variable, variable_assignment)
rm(command, command_pipe, int, operation)