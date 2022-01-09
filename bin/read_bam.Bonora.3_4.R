#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
#  ...from read_bam.Bonora.3.R

#  Load libraries
library(dplyr)
library(forcats)
library(magrittr)
library(parallel)
library(purrr)
library(Rsamtools)

options(pillar.sigfig = 8, scipen = 10000)


# -----------------------------------------------------------------------------
#  Temporary: Set up work directory (location TB∆)
"/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora" %>% setwd()


# -----------------------------------------------------------------------------
#  Set up functions
`%notin%` <- Negate(`%in%`)

makeOperation <- function(variable, command) {
    operation <- paste(variable, command) %>% paste(., collapse = "; ")
    return(operation)
}


# -----------------------------------------------------------------------------
#  Set up what to query from .bam files
# map_info <- c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar", "mrnm", "mpos", "isize", "seq", "qual")
map_info <- c("qname", "flag", "rname", "strand", "pos", "mapq", "cigar", "mrnm", "mpos", "isize", "seq")
tag_info <- c("AS", "XS", "NM")
map_params <- Rsamtools::ScanBamParam(what = map_info, tag = tag_info)

#  Assign variables necessary for subsequent commands
file <- list.files(pattern = "\\chr1.bam$")
file <- file[c(1:3)]
variable <- file %>% strsplit(., "-") %>% lapply(., `[[`, 1) %>% unlist() %>% paste0("dedup.", .)
command_pipe <- paste0("<- ", variable, " %>% ")
chromosome <- c(paste0("chr", 1:19), "chrX")

mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

#  Assign .bam information as list
command <- paste0(command_pipe, "Rsamtools::scanBam(., param = map_params)")
operation <- makeOperation(variable, command)
eval(parse(text = operation))
# rm(map_info, tag_info, map_params)

#  Convert .bam information from list to dataframe to tibble
command <- paste0(command_pipe, "as.data.frame() %>% as_tibble()")
operation <- makeOperation(variable, command)
eval(parse(text = operation))

#  Reorder rname factor levels
command <- paste0("<- forcats::fct_relevel(", variable, "$rname, chromosome)")
operation <- makeOperation(paste0(variable, "$rname"), command)
eval(parse(text = operation))

#  Drop unused rname factor levels
command <- paste0("<- ", variable, "$rname %>% forcats::fct_drop()")
operation <- makeOperation(paste0(variable, "$rname"), command)
eval(parse(text = operation))

#  Drop rows that are not chr1-19, chrX
command <- paste0(command_pipe, "filter(., rname %in% chromosome)")
operation <- makeOperation(variable, command)
eval(parse(text = operation))

#  Create and append pos_end vector
command <- paste0("<- ", variable, "$pos + 49")
operation <- makeOperation(paste0(variable, "$pos_end"), command)
eval(parse(text = operation))

#  Create and append mpos_end vector
command <- paste0("<- ", variable, "$mpos + 49")
operation <- makeOperation(paste0(variable, "$mpos_end"), command)
eval(parse(text = operation))

#  Reorder the pos_end and mpos_end columns; create a column of concatenated
#+ barcode_sequence strings; create a column of concatenated barcode_fragment
#+ (where "fragment" is rname, pos, and pos_end) strings; and then arrange
#+ columns by AS, rname, and pos
command <- paste0(
    command_pipe,
    "dplyr::relocate(pos_end, .after = pos) %>% dplyr::relocate(mpos_end, .after = mpos)",
    " %>% tidyr::unite(b_s, c(\"qname\", \"seq\"), sep = \"_\", remove = FALSE)",
    " %>% tidyr::unite(b_f, c(\"qname\", \"rname\", \"pos\", \"pos_end\"), sep = \"_\", remove = FALSE)",
    " %>% dplyr::arrange(-tag.AS, rname, pos, flag)"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))

#  Filter out all rows without unique b_s
command <- paste0(command_pipe, "dplyr::distinct(b_s, .keep_all = TRUE)")
operation <- makeOperation(variable, command)
eval(parse(text = operation))
# dedup.129S1 %>% nrow()  # [1] 530589
# dedup.CAST %>% nrow()  # [1] 528371
# dedup.mm10 %>% nrow()  # [1] 526011

# #  Filter out all rows without unique b_f
# command <- paste0(command_pipe, "dplyr::distinct(b_f, .keep_all = TRUE)")
# operation <- makeOperation(variable, command)
# eval(parse(text = operation))
# # dedup.129S1 %>% nrow()  # [1] 494694
# # dedup.CAST %>% nrow()  # [1] 492644
# # dedup.mm10 %>% nrow()  # [1] 491416

#  Create "minimal tibbles" ("mint") of only ... and AS values, and rename the
#+ AS columns to identify their alignment of origin
variable_assignment <- paste0("mint.", variable)
suffix <- variable_assignment %>% strsplit(., "\\.") %>% lapply(., `[[`, 3) %>% unlist()
command <- paste0(command_pipe, "dplyr::select(b_s, tag.AS) %>% dplyr::rename(AS.", suffix, " = tag.AS)")  #TODO b_f or b_s
operation <- makeOperation(variable_assignment, command)
eval(parse(text = operation))

#  Fully join the "mints" by ...; e.g., see the following URL:
#  https://stat545.com/join-cheatsheet.html#full_joinsuperheroes-publishers
mint.full <- full_join(mint.dedup.mm10, mint.dedup.129S1, by = "b_s") %>% 
    full_join(., mint.dedup.CAST, by = "b_s")  #TODO b_f or b_s

#  Calculate the differences between 129 and CAST alignment scores, then
#+ initially assign categories to the difference based on values w/r/t/the
#+ variable int (see below)
int <- 0 %>% as.integer()  #TODO Make the integer an argument
mint.full <- mint.full %>% 
    mutate(difference = AS.129S1 - AS.CAST) %>% 
    dplyr::mutate(
        assignment.initial = case_when(
            difference >= (-1 * int) & difference <= int ~ "Neutral",
            difference > int ~ "129S1-SvImJ",
            difference < (-1 * int) ~ "CAST-EiJ"
        )
    )

#  Formal assignments of categories: See the block of code
mint.full$assignment <- ifelse(
    is.na(mint.full$assignment.initial),
    ifelse(
        !is.na(mint.full$AS.129S1),
        "129S1-SvImJ",
        ifelse(
            !is.na(mint.full$AS.CAST),
            "CAST-EiJ",
            mint.full$assignment.initial
        )
    ),
    mint.full$assignment.initial
) %>%
    forcats::as_factor()

#  Remove all but b_f, AS, and formal assignment columns
mint.full <- mint.full %>% 
    select(b_s, AS.mm10, AS.129S1, AS.CAST, assignment)  #TODO b_f or b_s

#  Assign 0/1 categories for what aligned to what
mint.full$tmp.mm10 <- ifelse(
    is.na(mint.full$AS.mm10), "0", "1"
)

mint.full$tmp.129S1 <- ifelse(
    is.na(mint.full$AS.129S1), "0", "1"
)

mint.full$tmp.CAST <- ifelse(
    is.na(mint.full$AS.CAST), "0", "1"
)

mint.full$trinary <- paste0(
    mint.full$tmp.mm10,
    mint.full$tmp.129S1,
    mint.full$tmp.CAST
) %>%
    forcats::as_factor()

mint.full$assignment_trinary <- paste(
    mint.full$assignment,
    mint.full$trinary
) %>% as_factor()

#  Remove all but b_f, AS, and formal assignment, and "trinary" columns
mint.full <- mint.full %>%
    select(
        b_s,  #TODO b_f or b_s
        AS.mm10,
        AS.129S1,
        AS.CAST,
        assignment,
        trinary,
        assignment_trinary
    )

#  Create individual tibbles for each possible "assignment" value
mint.full$assignment <- mint.full$assignment %>% as.character()
mint.full$assignment[is.na(mint.full$assignment)] <- "NA"
variable <- list(
    c("assignment.129S1", "129S1-SvImJ"),
    c("assignment.CAST", "CAST-EiJ"),
    c("assignment.Neutral", "Neutral"),
    c("assignment.NA", "NA")
)
command <- paste0(
    "<- mint.full %>% dplyr::filter_at(vars(assignment), any_vars(. %in% \"",
    lapply(variable, `[[`, 2), "\"))"
)
operation <- makeOperation(lapply(variable, `[[`, 1), command)
eval(parse(text = operation))

#  Make $assignment factors again
command <- paste0("<- ", lapply(variable, `[[`, 1), "$assignment %>% as.factor()")
operation <- makeOperation(paste0(lapply(variable, `[[`, 1), "$assignment"), command)
eval(parse(text = operation))

assignment.NA$assignment <- NA
assignment.NA$assignment <- NA_character_

# rm(variable, mint.dedup.129S1, mint.dedup.CAST, mint.dedup.mm10, mint.full)

#  Full joins for 129S1 and CAST  #TODO mm10, NA, and Neutral
assignment.129S1 <- full_join(dedup.129S1, assignment.129S1, by = "b_s") %>%  #TODO b_f or b_s
    tidyr::drop_na("AS.129S1") %>% 
    dplyr::select(-dplyr::one_of("AS.129S1")) %>%
    dplyr::rename(
        AS.129S1 = tag.AS,
        XS.129S1 = tag.XS,
        NM.129S1 = tag.NM
    )

assignment.CAST <- full_join(dedup.CAST, assignment.CAST, by = "b_s") %>%  #TODO b_f or b_s
    tidyr::drop_na("AS.CAST") %>% 
    dplyr::select(-dplyr::one_of("AS.CAST")) %>%
    dplyr::rename(
        AS.CAST = tag.AS,
        XS.CAST = tag.XS,
        NM.CAST = tag.NM
    )


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
#  What to query from .bam files
# map_info <- c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar", "mrnm", "mpos", "isize", "seq", "qual")
map_info <- c("qname", "flag", "rname", "strand", "pos", "mapq", "cigar", "mrnm", "mpos", "isize", "seq")
tag_info <- c("AS", "XS", "NM")
map_params <- Rsamtools::ScanBamParam(what = map_info, tag = tag_info)

#  Assign variables necessary for subsequent commands; remember, "alt" is
#+ "CAST", "ref" is "129"
file <- list.files(pattern = "\\.chr1.bam$")
variable <- file %>% strsplit(., "\\.") %>% lapply(., `[[`, 5) %>% unlist() %>% paste0("GB.", .)
command_pipe <- paste0("<- ", variable, " %>% ")
chromosome <- c(paste0("chr", 1:19), "chrX")

mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

#  Assign .bam information as list
command <- paste0(command_pipe, "Rsamtools::scanBam(., param = map_params)")
operation <- makeOperation(variable, command)
eval(parse(text = operation))
# rm(map_info, tag_info, map_params)

#  Convert .bam information from list to dataframe to tibble
command <- paste0(command_pipe, "as.data.frame() %>% as_tibble()")
operation <- makeOperation(variable, command)
eval(parse(text = operation))

#  Reorder rname factor levels
command <- paste0("<- forcats::fct_relevel(", variable, "$rname, chromosome)")
operation <- makeOperation(paste0(variable, "$rname"), command)
eval(parse(text = operation))

#  Drop unused rname factor levels
command <- paste0("<- ", variable, "$rname %>% forcats::fct_drop()")
operation <- makeOperation(paste0(variable, "$rname"), command)
eval(parse(text = operation))

#  Drop rows that are not chr1-19, chrX
command <- paste0(command_pipe, "filter(., rname %in% chromosome)")
operation <- makeOperation(variable, command)
eval(parse(text = operation))

#  Create and append pos_end vector
command <- paste0("<- ", variable, "$pos + 49")
operation <- makeOperation(paste0(variable, "$pos_end"), command)
eval(parse(text = operation))

#  Create and append mpos_end vector
command <- paste0("<- ", variable, "$mpos + 49")
operation <- makeOperation(paste0(variable, "$mpos_end"), command)
eval(parse(text = operation))

#  Reorder the pos_end and mpos_end columns; create a column of concatenated
#+ barcode_sequence strings; create a column of concatenated barcode_fragment
#+ (where "fragment" is rname, pos, and pos_end) strings; and then arrange
#+ columns by AS, rname, and pos
command <- paste0(
    command_pipe,
    "dplyr::relocate(pos_end, .after = pos) %>% dplyr::relocate(mpos_end, .after = mpos)",
    " %>% tidyr::unite(b_s, c(\"qname\", \"seq\"), sep = \"_\", remove = FALSE)",
    " %>% tidyr::unite(b_f, c(\"qname\", \"rname\", \"pos\", \"pos_end\"), sep = \"_\", remove = FALSE)",
    " %>% dplyr::arrange(-tag.AS, rname, pos, flag)"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))

#  Filter out all rows without unique b_s
command <- paste0(command_pipe, "dplyr::distinct(b_s, .keep_all = TRUE)")
operation <- makeOperation(variable, command)
eval(parse(text = operation))
# GB.alt %>% nrow()  # [1] 1796387
# GB.ambig %>% nrow()  # [1] 5333701
# GB.contra %>% nrow()  # [1] 5882
# GB.ref %>% nrow()  # [1] 1834526

# #  Filter out all rows without unique b_f
# command <- paste0(command_pipe, "dplyr::distinct(b_f, .keep_all = TRUE)")
# operation <- makeOperation(variable, command)
# eval(parse(text = operation))
# # GB.alt %>% nrow()  # [1] 1719080
# # GB.ambig %>% nrow()  # [1] 4923186
# # GB.contra %>% nrow()  # [1] 5581
# # GB.ref %>% nrow()  # [1] 1754655

#  Remember, "alt" is "CAST", "ref" is "129"
GB.alt.CAST <- GB.alt
GB.ref.129S1 <- GB.ref

#  Clean up
#TODO Determine better spots fot this...
# rm(GB.alt, GB.ref)
# rm(chromosome, file, suffix, variable, variable_assignment)
# rm(command, command_pipe, int, operation)
