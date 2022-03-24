#!/usr/bin/env Rscript

#  deduplicate-with-respect-to-AS-MAPQ.R
#  KA


#  Source packages, sey options -----------------------------------------------
library(argparser)
library(Rsamtools)
library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)


#  Load in functions ----------------------------------------------------------
assignFilesToVariables <- function(file, string) {
    split <- file %>% strsplit(., "/")
    name <- split %>% 
        lapply(., `[[`, length(split[[1]])) #%>%
        # unlist() %>%
        # strsplit(., "-") %>%
        # lapply(., `[[`, 1) %>%
        # unlist() %>%
        # paste0(string, .) %>%
        # stringr::str_remove("S1")
    return(name)
}


evaluateOperation <- function(operation = operation) {
    return(eval(parse(text = operation), envir = .GlobalEnv))
}


makeOperation <- function(variable = variable, command = command) {
    operation <- paste0(variable, " <- ", command)
    return(operation)
}


#  Parse arguments ------------------------------------------------------------
script <- "deduplicate-with-respect-to-AS-MAPQ.R"

#  Create a parser
ap <- arg_parser(
    name = script,
    description = "",
    hide.opts = TRUE
)

#  Add command line arguments
ap <- add_argument(
    ap,
    short = "-i",
    arg = "--infile",
    type = "character",
    default = NULL,
    help = "bam infile, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-b",
    arg = "--bai",
    type = "character",
    default = NULL,
    help = "bam index, including path <chr>"
)
ap <- add_argument(
    ap,
    short = "-o",
    arg = "--outdir",
    type = "character",
    default = NULL,
    help = "directory for saving outfile(s), including path <chr>"
)

#  Parse the command line arguments
directory_base <- "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
directory_data <- "data"
directory_out <- paste0(directory_base, "/", directory_data)
infile <- paste0(
    directory_out, "/",
    "CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam"
)
bai <- paste0(
    directory_out, "/",
    "CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam.bai"
)
cl <- c(
    #  Arguments for analysis
    "--infile", infile,
    "--bai", bai,
    "--outdir", directory_out
)
arguments <- parse_args(ap, cl)  # RStudio-interactive work
# arguments <- parse_args(ap)  # Command-line calls

rm(directory_base, directory_data, directory_out, infile, bai)


#  Load in .bam information, including mate information -----------------------
file <- c(arguments$infile)
variable <- assignFilesToVariables(file)
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

file <- c(arguments$bai)
index <- assignFilesToVariables(file)
mapply(
    assign, index, file, MoreArgs = list(envir = parent.frame())
)

#  Using Rsamtools, load in standard .bam fields
command <- paste0(
    "Rsamtools::BamFile(",
        variable, ", index = ", index, ", asMates = TRUE",
    ")", " %>% ",
    "Rsamtools::scanBam()", " %>% ",
    "as.data.frame()", " %>% ",
    "tibble::as_tibble()"  #TODO Speed up via parallelization?
)
operation <- makeOperation("full", command)
evaluateOperation(operation)

#  Using Rsamtools, load in .bam AS and MD fields
map_params <- Rsamtools::ScanBamParam(tag = c("AS", "MD"))
command <- paste0(
    variable, " %>% ",
        "Rsamtools::BamFile(",
            variable, ", index = ", index, ", asMates = TRUE",
        ")", " %>% ",
        "Rsamtools::scanBam(param = map_params)", " %>% ",
        "as.data.frame()", " %>% ",
        "tibble::as_tibble()"  #TODO Speed up via parallelization?
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

rm(map_params)

#  Column-bind the standard, AS, and MD fields
command <- paste0(
    "dplyr::bind_cols(", variable, ".full, ", variable, ")"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  Remove unneeded variables
operation <- paste0("rm(", variable, ".full, ", index, ")")
evaluateOperation(operation)

rm(file, index)


#  Goal: In - m; out - df.2, df.3 ---------------------------------------------
#+ 
#+ If, for a group of duplicate rows, there are no differences in tag.AS or
#+ mapq, then pick one row from the group at random; if there are differences,
#+ then pick the row with the better tag.AS; if tag.AS are equal, then pick the
#+ row with the better mapq

#  Based on coordinate value, identify if a given row is a duplicate 
df.main$duplicated <- ave(
    df.main$coordinate, df.main$coordinate, FUN = length
) > 1L

#  Create separate df's for duplicated and non-duplicated observations
df.main.split <- df.main %>%
    dplyr::group_by(duplicated) %>%
    dplyr::group_split()
df.nondup <- df.main.split[[1]]
df <- df.main.split[[2]]
rm(df.main, df.main.split)  # Free up memory

#  Create categories based on the number of rows comprising a "duplicate group"
df <- df %>% 
    group_by(coordinate) %>% 
    dplyr::mutate(group_n_members = n()) %>%
    dplyr::mutate_at(vars(group_n_members), dplyr::funs(factor)) %>%
    dplyr::mutate(group_n_members = recode(
        group_n_members, "1" = "one", "2" = "two", "3" = "three"
    ))
#TODO dplyr::funs() is deprecated; replace it with another function

#  Create separate m's for "duplicate groups" comprised of two members and
#+ "duplicate groups" comprised of three members
df.split <- df %>%
    dplyr::group_by(group_n_members) %>%
    dplyr::group_split()
df.2 <- df.split[[1]]
df.3 <- df.split[[2]]
rm(df, df.split)


#  Goal: In - df.3; out - df.3.same, df.3.diff --------------------------------
#+ 
#+ If tag.AS is the same for each member of the group, then "same;" else "diff"
df.3 <- df.3 %>%
    dplyr::group_by(coordinate) %>%
    dplyr::mutate(group_AS_same = dplyr::if_else(
        dplyr::n_distinct(tag.AS) == 1, "same", "diff"
    )) %>%
    dplyr::mutate_at(vars(group_AS_same), dplyr::funs(factor))

df.split <- df.3 %>%
    dplyr::group_by(group_AS_same) %>%
    dplyr::group_split()
df.3.same <- df.split[[1]]
df.3.diff <- df.split[[2]]
rm(df.3, df.split)


#  Goal: In - df.3.diff; out - df.3.diff.n2, df.3.diff.n3 ---------------------
#+ 
#+ Separate them by how many differences they have; there should only be two
#+ options for differences: n = 2 or n = 3 differences
df.3.diff <- df.3.diff %>%
    group_by(coordinate) %>%
    dplyr::mutate(group_AS_n_distinct = dplyr::if_else(
        dplyr::n_distinct(tag.AS) == 2, "two", "three"
    )) %>%
    dplyr::mutate_at(vars(group_AS_same), dplyr::funs(factor))

df.split <- df.3.diff %>%
    dplyr::group_by(group_AS_n_distinct) %>%
    dplyr::group_split()
df.3.diff.n2 <- df.split[[2]]  #TODO What to do if no "n2?" For example, this is causing "chr2" to fail...
df.3.diff.n3 <- df.split[[1]]  #TODO Some kind of "check" system for each of these splits and assignments
rm(df.3.diff, df.split)


#  Goal: In - df.3.diff.n3; out - df.3.diff.n3.dedup --------------------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has the
#+ best tag.AS score (and thus also the best or same mapq)
df.3.diff.n3.dedup <- df.3.diff.n3[!duplicated(df.3.diff.n3$coordinate), ]
rm(df.3.diff.n3)


#  Goal: In - df.3.diff.n2; out - df.3.diff.n2.cut, *.cut.n1, *.cut.n2 --------
#+ 
#+ Per group, remove the row with the worst (lowest) tag.AS score; then,
#+ categorize group_AS following the cut of 1/3 of the tibble
df.3.diff.n2.cut <- df.3.diff.n2 %>%
    group_by(coordinate) %>%
    slice_max(tag.AS, n = 2, with_ties = FALSE) %>%
    dplyr::mutate(group_AS_n_distinct.post_cut = dplyr::if_else(
        dplyr::n_distinct(tag.AS) == 1, "same", "diff"
    )) %>%
    dplyr::mutate_at(vars(group_AS_n_distinct.post_cut), dplyr::funs(factor))

df.split <- df.3.diff.n2.cut %>%
    dplyr::group_by(group_AS_n_distinct.post_cut) %>%
    dplyr::group_split()
df.3.diff.n2.cut.same <- df.split[[1]]
df.3.diff.n2.cut.diff <- df.split[[2]]
rm(df.3.diff.n2, df.3.diff.n2.cut, df.split)


#  Goal: In - df.3.diff.n2.cut.same; out - *.n1.mq_same, *.n1.mq_diff ---------
#+ 
#+ Select for mapq before picking at random
#TODO Better description
df.3.diff.n2.cut.same <- df.3.diff.n2.cut.same %>%
    dplyr::group_by(coordinate) %>%
    dplyr::mutate(group_MAPQ_same = dplyr::if_else(
        dplyr::n_distinct(mapq) == 1, "same", "diff"
    )) %>%
    dplyr::mutate_at(vars(group_MAPQ_same), dplyr::funs(factor))

df.split <- df.3.diff.n2.cut.same %>%
    dplyr::group_by(group_MAPQ_same) %>%
    dplyr::group_split()
df.3.diff.n2.cut.same.mq_same <- df.split[[1]]
df.3.diff.n2.cut.same.mq_diff <- df.split[[2]]
rm(df.3.diff.n2.cut.same, df.split)


#  Goal: In - *.n1.mq_diff; out - *.n1.mq_diff.dedup --------------------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has both
#+ the best tag.AS score and best mapq score
df.3.diff.n2.cut.same.mq_diff.dedup <-
    df.3.diff.n2.cut.same.mq_diff[
        !duplicated(df.3.diff.n2.cut.same.mq_diff$coordinate), 
    ]
rm(df.3.diff.n2.cut.same.mq_diff)


#  Goal: In - *.n1.mq_same; out - *.n1.mq_same.dedup --------------------------
#+ !DEDUP
#+ 
#+ Per group, pick one row at random; per group, there are no differences in
#+ tag.AS and mapq scores
set.seed(24)
df.3.diff.n2.cut.same.mq_same.dedup <- df.3.diff.n2.cut.same.mq_same %>%
    group_by(coordinate) %>%
    slice_sample(n = 1)

rm(df.3.diff.n2.cut.same.mq_same)


#  Goal: In - df.3.diff.n2.cut.diff; out - df.3.diff.n2.cut.diff.dedup --------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has the
#+ best tag.AS score (and thus also the best or same mapq)
df.3.diff.n2.cut.diff.dedup <-
    df.3.diff.n2.cut.diff[
        !duplicated(df.3.diff.n2.cut.diff$coordinate), 
    ]
rm(df.3.diff.n2.cut.diff)


#  Goal: In - df.3.same; out - df.3.same.mq_same, df.3.same.mq_diff -----------
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has the
#+ best tag.AS score
df.3.same <- df.3.same %>%
    dplyr::group_by(coordinate) %>%
    dplyr::mutate(group_MAPQ_same = dplyr::if_else(
        dplyr::n_distinct(mapq) == 1, "same", "diff"
    )) %>%
    dplyr::mutate_at(vars(group_MAPQ_same), dplyr::funs(factor))
#TODO dplyr::funs() is deprecated; replace it with another function

df.split <- df.3.same %>%
    dplyr::group_by(group_MAPQ_same) %>%
    dplyr::group_split()
df.3.same.mq_same <- df.split[[1]]
df.3.same.mq_diff <- df.split[[2]]
rm(df.3.same, df.split)


#  Goal: In - df.3.same.mq_diff; out - df.3.same.mq_diff.dedup ----------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has both
#+ the best tag.AS score and best mapq score
df.3.same.mq_diff.dedup <-
    df.3.same.mq_diff[
        !duplicated(df.3.same.mq_diff$coordinate), 
    ]
rm(df.3.same.mq_diff)


#  Goal: In - df.3.same.mq_same; out - df.3.same.mq_same.dedup ----------------
#+ !DEDUP
#+ 
#+ Per group, pick one row at random; now, per group, there are no differences
#+ in tag.AS and mapq scores
set.seed(24)
df.3.same.mq_same.dedup <- df.3.same.mq_same %>%
    group_by(coordinate) %>%
    slice_sample(n = 1)

# rm(df.3.same.mq_same)


#  Goal: In - df.2; out - df.2.same, df.2.diff --------------------------------
#+ 
#+ If tag.AS is the same for each member of the group, then "same;" else "diff"
df.2 <- df.2 %>%
    dplyr::group_by(coordinate) %>%
    dplyr::mutate(group_AS_same = dplyr::if_else(
        dplyr::n_distinct(tag.AS) == 1, "same", "diff"
    )) %>%
    dplyr::mutate_at(vars(group_AS_same), dplyr::funs(factor))

df.split <- df.2 %>%
    dplyr::group_by(group_AS_same) %>%
    dplyr::group_split()
df.2.same <- df.split[[1]]
df.2.diff <- df.split[[2]]
rm(df.2, df.split)


#  Goal: In - df.2.diff; out - df.2.diff.dedup --------------------------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has the
#+ best tag.AS score
df.2.diff.dedup <- df.2.diff[!duplicated(df.2.diff$coordinate), ]
rm(df.2.diff)


#  Goal: In - df.2.same; out - df.2.same.mq_same, df.2.same.mq_diff -----------
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has the
#+ best tag.AS score
df.2.same <- df.2.same %>%
    dplyr::group_by(coordinate) %>%
    dplyr::mutate(group_MAPQ_same = dplyr::if_else(
        dplyr::n_distinct(mapq) == 1, "same", "diff"
    )) %>%
    dplyr::mutate_at(vars(group_MAPQ_same), dplyr::funs(factor))

df.split <- df.2.same %>%
    dplyr::group_by(group_MAPQ_same) %>%
    dplyr::group_split()
df.2.same.mq_same <- df.split[[1]]
df.2.same.mq_diff <- df.split[[2]]
rm(df.2.same, df.split)


#  Goal: In - df.2.same.mq_diff; out - df.2.same.mq_diff.dedup ----------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has both
#+ the best tag.AS score and best mapq score
df.2.same.mq_diff.dedup <-
    df.2.same.mq_diff[!duplicated(df.2.same.mq_diff$coordinate), ]
rm(df.2.same.mq_diff)


#  Goal: In - df.2.same.mq_same; out - df.2.same.mq_same.dedup ----------------
#+ !DEDUP
#+ 
#+ Per group, pick one row at random; per group, there are no differences in
#+ tag.AS and mapq scores
set.seed(24)
df.2.same.mq_same.dedup <- df.2.same.mq_same %>%
    group_by(coordinate) %>%
    slice_sample(n = 1)
rm(df.2.same.mq_same)
