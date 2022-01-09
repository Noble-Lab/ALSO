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
`%notin%` <- Negate(`%in%`)

makeOperation <- function(variable, command) {
    operation <- paste(variable, command) %>% paste(., collapse = "; ")
    return(operation)
}

#  What to query from .bam files
map_info <- c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar", "mrnm", "mpos", "isize", "seq", "qual")
tag_info <- c("AS", "XS", "NM")
map_params <- Rsamtools::ScanBamParam(what = map_info, tag = tag_info)

#  Assign variables necessary for subsequent commands
files <- list.files(pattern = "\\chr1.bam$")
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
dedup.Nmasked <- bam4 %>% Rsamtools::scanBam(., param = map_params)

#  Convert .bam information from list to dataframe to tibble
variables <- c(
    "dedup.129",
    "dedup.CAST",
    "dedup.Nmasked"
)

command_pipe <- paste0("<- ", variables, " %>% ")

command <- paste0(command_pipe, "as.data.frame() %>% as_tibble()")
operation <- makeOperation(variables, command)
eval(parse(text = operation))
rm("bam1", "bam2", "bam3", "bam4", "bam5")

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

dedup.129 %>% head()
dedup.CAST %>% head()
dedup.Nmasked %>% head()


# -----------------------------------------------------------------------------
#  Create a column of concatenated barcode_fragment strings; arrange data
#+ frames by rname, pos, and AS; then, filter out all rows without non-unique
#+ combinations of qname (barcode) and seq (read sequence) values
b_f.dedup.129 <- dedup.129 %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    ) %>%
    dplyr::arrange(rname, pos, tag.AS) %>% 
    dplyr::distinct(qname, seq, .keep_all = TRUE)

b_f.dedup.CAST <- dedup.CAST %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    ) %>%
    dplyr::arrange(rname, pos, tag.AS) %>% 
    dplyr::distinct(qname, seq, .keep_all = TRUE)

b_f.dedup.Nmasked <- dedup.Nmasked %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    ) %>%
    dplyr::arrange(rname, pos, tag.AS) %>% 
    dplyr::distinct(qname, seq, .keep_all = TRUE)

# setdiff.129.CAST <- setdiff(b_f.dedup.129$b_f, b_f.dedup.CAST$b_f)
# setdiff.129.Nmasked <- setdiff(b_f.dedup.129$b_f, b_f.dedup.Nmasked$b_f)
# setdiff.CAST.129 <- setdiff(b_f.dedup.CAST$b_f, b_f.dedup.129$b_f)
# setdiff.CAST.Nmasked <- setdiff(b_f.dedup.CAST$b_f, b_f.dedup.Nmasked$b_f)
# setdiff.Nmasked.CAST <- setdiff(b_f.dedup.Nmasked$b_f, b_f.dedup.CAST$b_f)
# setdiff.Nmasked.129 <- setdiff(b_f.dedup.Nmasked$b_f, b_f.dedup.129$b_f)


# -----------------------------------------------------------------------------
#  Create tables of only barcode_fragment and AS values
b_f.AS.dedup.129 <- b_f.dedup.129 %>% select(b_f, tag.AS)
b_f.AS.dedup.CAST <- b_f.dedup.CAST %>% select(b_f, tag.AS)
b_f.AS.dedup.Nmasked <- b_f.dedup.Nmasked %>% select(b_f, tag.AS)

#  Rename the AS columns to identify their alignment of origin
b_f.AS.dedup.129 <- b_f.AS.dedup.129 %>%
    dplyr::rename(
        AS.129 = tag.AS,
        # mm10.XS = tag.XS
    )
b_f.AS.dedup.CAST <- b_f.AS.dedup.CAST %>%
    dplyr::rename(
        AS.CAST = tag.AS,
        # CAST.XS = tag.XS
    )
b_f.AS.dedup.Nmasked <- b_f.AS.dedup.Nmasked %>%
    dplyr::rename(
        AS.Nmasked = tag.AS,
        # mm10.XS = tag.XS
    )

#  Fully join the tables of only barcode_fragment and AS values; e.g., see
#  https://stat545.com/join-cheatsheet.html#full_joinsuperheroes-publishers
b_f.AS.dedup.Nmasked.129.CAST <- full_join(b_f.AS.dedup.Nmasked, b_f.AS.dedup.129, by = "b_f") %>% 
    full_join(., b_f.AS.dedup.CAST, by = "b_f")

# -----------------
int <- 0 %>% as.integer()  #TODO Make the integer an argument

#  Calculate the differences between 129 and CAST alignment scores, then
#+ initially assign categories to the difference based on values w/r/t/the
#+ variable int (see below)
b_f.AS.dedup.Nmasked.129.CAST <- b_f.AS.dedup.Nmasked.129.CAST %>% 
    mutate(difference = AS.129 - AS.CAST) %>% 
    dplyr::mutate(
        assignment.initial = case_when(
            difference >= (-1 * int) & difference <= int ~ "Neutral",
            difference > int ~ "129S1-SvImJ",
            difference < (-1 * int) ~ "CAST-EiJ"
        )
    )

#  Formal assignments of categories: See the block of code
b_f.AS.dedup.Nmasked.129.CAST$assignment <- ifelse(
    is.na(b_f.AS.dedup.Nmasked.129.CAST$assignment.initial),
    ifelse(
        !is.na(b_f.AS.dedup.Nmasked.129.CAST$AS.129),
        "129S1-SvImJ",
        ifelse(
            !is.na(b_f.AS.dedup.Nmasked.129.CAST$AS.CAST),
            "CAST-EiJ",
            b_f.AS.dedup.Nmasked.129.CAST$assignment.initial
        )
    ),
    b_f.AS.dedup.Nmasked.129.CAST$assignment.initial
) %>%
    forcats::as_factor()

#  Remove all but b_f, AS, and formal assignment columns
b_f.AS.dedup.Nmasked.129.CAST <- b_f.AS.dedup.Nmasked.129.CAST %>% 
    select(b_f, AS.Nmasked, AS.129, AS.CAST, assignment)

#  Assign 0/1 categories for what aligned to what
b_f.AS.dedup.Nmasked.129.CAST$tmp.mm10_CAST_129_Nmasked <- ifelse(
    is.na(b_f.AS.dedup.Nmasked.129.CAST$AS.Nmasked), "0", "1"
)

b_f.AS.dedup.Nmasked.129.CAST$tmp.129 <- ifelse(
    is.na(b_f.AS.dedup.Nmasked.129.CAST$AS.129), "0", "1"
)

b_f.AS.dedup.Nmasked.129.CAST$tmp.CAST <- ifelse(
    is.na(b_f.AS.dedup.Nmasked.129.CAST$AS.CAST), "0", "1"
)

b_f.AS.dedup.Nmasked.129.CAST$trinary <- paste0(
    b_f.AS.dedup.Nmasked.129.CAST$tmp.mm10_CAST_129_Nmasked,
    b_f.AS.dedup.Nmasked.129.CAST$tmp.129,
    b_f.AS.dedup.Nmasked.129.CAST$tmp.CAST
) %>%
    forcats::as_factor()

b_f.AS.dedup.Nmasked.129.CAST$assignment_trinary <- paste(
    b_f.AS.dedup.Nmasked.129.CAST$assignment,
    b_f.AS.dedup.Nmasked.129.CAST$trinary
) %>% as_factor()

#  Remove all but b_f, AS, and formal assignment, and "trinary" columns
b_f.AS.dedup.Nmasked.129.CAST <- b_f.AS.dedup.Nmasked.129.CAST %>%
    select(
        b_f,
        AS.Nmasked,
        AS.129,
        AS.CAST,
        assignment,
        trinary,
        assignment_trinary
    )

#  Split data frame by "assignment" value
split.b_f.AS.dedup.Nmasked.129.CAST <- b_f.AS.dedup.Nmasked.129.CAST %>% 
    group_split(assignment)

assignment.129 <- split.b_f.AS.dedup.Nmasked.129.CAST[[2]]
assignment.CAST <- split.b_f.AS.dedup.Nmasked.129.CAST[[3]]
assignment.Neutral <- split.b_f.AS.dedup.Nmasked.129.CAST[[1]]
assignment.NA <- split.b_f.AS.dedup.Nmasked.129.CAST[[4]]

assignment.129 <- full_join(b_f.dedup.129, assignment.129, by = "b_f")
assignment.129 <- assignment.129 %>% 
    tidyr::drop_na("AS.129") %>% 
    dplyr::select(-dplyr::one_of("AS.129")) %>%
    dplyr::rename(
        AS.129 = tag.AS,
        XS.129 = tag.XS,
        NM.129 = tag.NM
    )

assignment.CAST <- full_join(b_f.dedup.CAST, assignment.CAST, by = "b_f")
assignment.CAST <- assignment.CAST %>% 
    tidyr::drop_na("AS.CAST") %>% 
    dplyr::select(-dplyr::one_of("AS.CAST")) %>%
    dplyr::rename(
        AS.CAST = tag.AS,
        XS.CAST = tag.XS,
        NM.CAST = tag.NM
    )

# assignment.Neutral.Nmasked <- full_join(b_f.dedup.Nmasked, assignment.Neutral, by = "b_f")
# assignment.Neutral.Nmasked <- assignment.Neutral.Nmasked %>% 
#     tidyr::drop_na("AS.Nmasked") %>% 
#     dplyr::select(-dplyr::one_of("AS.Nmasked")) %>%
#     dplyr::rename(
#         AS.Nmasked = tag.AS,
#         XS.Nmasked = tag.XS,
#         NM.Nmasked = tag.NM
#     ) %>%
#     dplyr::arrange(rname, pos)
# 
# assignment.Neutral.129 <- full_join(b_f.dedup.129, assignment.Neutral, by = "b_f")
# assignment.Neutral.129 <- assignment.Neutral.129 %>% 
#     tidyr::drop_na("AS.129") %>% 
#     dplyr::select(-dplyr::one_of("AS.129")) %>%
#     dplyr::rename(
#         AS.129 = tag.AS,
#         XS.129 = tag.XS,
#         NM.129 = tag.NM
#     ) %>%
#     dplyr::arrange(rname, pos)
# 
# assignment.Neutral.CAST <- full_join(b_f.dedup.CAST, assignment.Neutral, by = "b_f")
# assignment.Neutral.CAST <- assignment.Neutral.CAST %>% 
#     tidyr::drop_na("AS.CAST") %>% 
#     dplyr::select(-dplyr::one_of("AS.CAST")) %>%
#     dplyr::rename(
#         AS.CAST = tag.AS,
#         XS.CAST = tag.XS,
#         NM.CAST = tag.NM
#     ) %>%
#     dplyr::arrange(rname, pos)
# 
# assignment.NA.Nmasked <- full_join(b_f.dedup.Nmasked, assignment.NA, by = "b_f")
# assignment.NA.Nmasked <- assignment.NA.Nmasked %>% 
#     tidyr::drop_na("AS.Nmasked") %>% 
#     dplyr::select(-dplyr::one_of("AS.Nmasked")) %>%
#     dplyr::rename(
#         AS.Nmasked = tag.AS,
#         XS.Nmasked = tag.XS,
#         NM.Nmasked = tag.NM
#     ) %>%
#     dplyr::arrange(rname, pos)
# 
# assignment.NA.129 <- full_join(b_f.dedup.129, assignment.NA, by = "b_f")
# assignment.NA.129 <- assignment.NA.129 %>% 
#     tidyr::drop_na("AS.129") %>% 
#     dplyr::select(-dplyr::one_of("AS.129")) %>%
#     dplyr::rename(
#         AS.129 = tag.AS,
#         XS.129 = tag.XS,
#         NM.129 = tag.NM
#     ) %>%
#     dplyr::arrange(rname, pos)
# 
# assignment.NA.CAST <- full_join(b_f.dedup.CAST, assignment.NA, by = "b_f")
# assignment.NA.CAST <- assignment.NA.CAST %>% 
#     tidyr::drop_na("AS.CAST") %>% 
#     dplyr::select(-dplyr::one_of("AS.CAST")) %>%
#     dplyr::rename(
#         AS.CAST = tag.AS,
#         XS.CAST = tag.XS,
#         NM.CAST = tag.NM
#     ) %>%
#     dplyr::arrange(rname, pos)
