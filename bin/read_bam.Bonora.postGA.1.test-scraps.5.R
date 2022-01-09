#!/usr/bin/env Rscript

packages <- c(
    "Rsamtools",
    "tidyverse"
)

for(i in packages){
    suppressPackageStartupMessages(library(i, character.only = TRUE))
}
rm(i, packages)

options(pillar.sigfig = 8, scipen = 10000)


# Set up work directory (location TB∆) ----------------------------------------
directory_user <- "/Users/kalavattam"
directory_base <- "Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do"
directory_work <- "get_unique_fragments/Bonora/segregatedReads.SNPTHRESH1.Q30"
setwd(
    paste0(directory_user, "/", directory_base, "/", directory_work)
)
rm(directory_base, directory_work)


#  Files are from...
#+ /net/noble/vol2/home/gbonora/proj/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouseDiff_Nmasked/segregatedReads.SNPTHRESH1.Q30
#+ 
#+ For more information, see the following script:
#+ /net/noble/vol2/home/gbonora/proj/2019_sciATAC_analysis/results/gbonora/20191105_sciATAC_mouseDiff_Nmasked/20191105_sciATAC_mouseDiff_Nmasked_allelicSegregation_workflow.sh


#  Set up functions -----------------------------------------------------------
`%notin%` <- Negate(`%in%`)


loadRDS <- readRDS


loadTibbleFromRDS <- function(variable, file) {
    command <- paste0(variable, " <- loadRDS(file = \"", file, "\")") %>% as.list()
    return(eval(parse(text = command), envir = globalenv()))
}


evaluateOperation <- function(operation = operation) {
    return(eval(parse(text = operation), envir = .GlobalEnv))
}


makeOperation <- function(variable, command) {
    operation <- paste(variable, command) #%>% paste(., collapse = "; ")
    return(operation)
}


# for (chromosome in c(paste0("chr", 2:4), "chrX")) {

#  Load tibbles (already sorted) from .rds files ------------------------------
#+ 
#+ .rds files come from convert_bam-GB_tibble.R

chromosome <- "chrX"
file <- list.files(pattern = paste0("\\.", chromosome, ".rds$"))
variable <- file %>% gsub("\\.rds$", "", .) %>% str_subset(., "^GB")
loadTibbleFromRDS(variable = variable, file = file)
#NOTE Remember, "alt" is "CAST", "ref" is "129"


#  Append information about the GB assignment ---------------------------------
command <- paste0("<- \"", gsub("\\.chr*(.+)$", "", variable), "\"")
operation <- makeOperation(paste0(variable, "$assignment_GB"), command)
evaluateOperation(operation)


#  Removal of Picard MarkDuplicates flags -------------------------------------
duplicate_flags <- c(1107, 1123, 1171, 1187)
command <- paste0(
    "<- ", variable,
        " %>% dplyr::filter(flag %notin% duplicate_flags)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Deduplication of the individual tibbles (already appropriately sorted) -----
variable_dedup <- paste0("dedup.", variable)
command <- paste0("<- ", variable, "[!duplicated(", variable, "$coordinate), ]")
operation <- makeOperation(variable_dedup, command)
evaluateOperation(operation)


#  Free up memory by removing the non-deduplicated tibbles --------------------
command <- paste0("rm(", variable, ")")
evaluateOperation(command)


#  Joint deduplication of the tibbles -----------------------------------------
#+ Create a vector of the four tibble variables
command <- paste0(
    "<- dplyr::bind_rows(",
        variable_dedup[1], ", ",
        variable_dedup[2], ", ",
        variable_dedup[3], ", ",
        variable_dedup[4],
    ")"
)
operation <- makeOperation("joint.GB", command) 
evaluateOperation(operation)

#+ Free up memory by unloading large tibbles
command <- paste0("rm(", variable_dedup,")")
evaluateOperation(command)

#+ Order the joint tibble by coordinate, AS, then mapq
joint.GB <- joint.GB[order(joint.GB$coordinate, -joint.GB$tag.AS -joint.GB$mapq), ]


#  Goal: In - m; out - j.2, j.3 -----------------------------------------------
#+ 
#+ If, for a group of duplicate rows, there are no differences in tag.AS or
#+ mapq, then pick one row from the group at random; if there are differences,
#+ then pick the row with the better tag.AS; if tag.AS are equal, then pick the
#+ row with the better mapq

#  Based on coordinate value, identify if a given row is a duplicate 
joint.GB$duplicated <- ave(
    joint.GB$coordinate, joint.GB$coordinate, FUN = length
) > 1L

#  Create separate j's for duplicated and non-duplicated observations
joint.GB.split <- joint.GB %>%
    dplyr::group_by(duplicated) %>%
    dplyr::group_split()
j.nondup <- joint.GB.split[[1]]
j <- joint.GB.split[[2]]
# rm(joint.GB, joint.GB.split)  # Free up memory

#  Create categories based on the number of rows comprising a "duplicate group"
j <- j %>% 
    group_by(coordinate) %>% 
    dplyr::mutate(group_n_members = n())  %>%
    dplyr::mutate_at(vars(group_n_members), dplyr::funs(factor)) %>%
    dplyr::mutate(group_n_members = recode(
        group_n_members, "1" = "one", "2" = "two", "3" = "three"
    ))
#TODO dplyr::funs() is deprecated; replace it with another function

#  Create separate m's for "duplicate groups" comprised of two members and
#+ "duplicate groups" comprised of three members
j.split <- j %>%
    dplyr::group_by(group_n_members) %>%
    dplyr::group_split()
j.2 <- j.split[[1]]
j.3 <- j.split[[2]]
# rm(j, j.split)


#  Goal: In - j.3; out - j.3.same, j.3.diff ----------------------------------- 
#+ 
#+ If tag.AS is the same for each member of the group, then "same;" else "diff"
j.3 <- j.3 %>%
    dplyr::group_by(coordinate) %>%
    dplyr::mutate(group_AS_same = dplyr::if_else(
        dplyr::n_distinct(tag.AS) == 1, "same", "diff"
    )
    ) %>%
    dplyr::mutate_at(vars(group_AS_same), dplyr::funs(factor))

j.split <- j.3 %>%
    dplyr::group_by(group_AS_same) %>%
    dplyr::group_split()
j.3.same <- j.split[[1]]
j.3.diff <- j.split[[2]]
# rm(j.3, j.split)


#  Goal: In - j.3.diff; out - j.3.diff.n2, j.3.diff.n3 ------------------------
#+ 
#+ Separate them by how many differences they have; there should only be two
#+ options for differences: n = 2 or n = 3 differences
j.3.diff <- j.3.diff %>%
    group_by(coordinate) %>%
    dplyr::mutate(group_AS_n_distinct = dplyr::if_else(
        dplyr::n_distinct(tag.AS) == 2, "two", "three"
    )) %>%
    dplyr::mutate_at(vars(group_AS_same), dplyr::funs(factor))

j.split <- j.3.diff %>%
    dplyr::group_by(group_AS_n_distinct) %>%
    dplyr::group_split()
j.3.diff.n2 <- j.split[[2]]  #TODO What to do if no "n2?" For example, this is causing "chr2" to fail...
j.3.diff.n3 <- j.split[[1]]  #TODO Some kind of "check" system for each of these splits and assignments
# rm(j.3.diff, j.split)


#  Goal: In - j.3.diff.n3; out - j.3.diff.n3.dedup ----------------------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has the
#+ best tag.AS score (and thus also the best or same mapq)
j.3.diff.n3.dedup <- j.3.diff.n3[!duplicated(j.3.diff.n3$coordinate), ]

# rm(j.3.diff.n3)


#  Goal: In - j.3.diff.n2; out - j.3.diff.n2.cut, *.cut.n1, *.cut.n2 ----------
#+ 
#+ Per group, remove the row with the worst (lowest) tag.AS score; then,
#+ categorize group_AS following the cut of 1/3 of the tibble
j.3.diff.n2.cut <- j.3.diff.n2 %>%
    group_by(coordinate) %>%
    slice_max(tag.AS, n = 2, with_ties = FALSE) %>%
    dplyr::mutate(group_AS_n_distinct.post_cut = dplyr::if_else(
        dplyr::n_distinct(tag.AS) == 1, "same", "diff"
    )) %>%
    dplyr::mutate_at(vars(group_AS_n_distinct.post_cut), dplyr::funs(factor))

j.split <- j.3.diff.n2.cut %>%
    dplyr::group_by(group_AS_n_distinct.post_cut) %>%
    dplyr::group_split()
j.3.diff.n2.cut.same <- j.split[[1]]
j.3.diff.n2.cut.diff <- j.split[[2]]
# rm(j.3.diff.n2, j.3.diff.n2.cut, j.split)


#  Goal: In - j.3.diff.n2.cut.same; out - *.n1.mq_same, *.n1.mq_diff ------------
#+ 
#+ Select for mapq before picking at random
#TODO Better description
j.3.diff.n2.cut.same <- j.3.diff.n2.cut.same %>%
    dplyr::group_by(coordinate) %>%
    dplyr::mutate(group_MAPQ_same = dplyr::if_else(
        dplyr::n_distinct(mapq) == 1, "same", "diff"
    )) %>%
    dplyr::mutate_at(vars(group_MAPQ_same), dplyr::funs(factor))

j.split <- j.3.diff.n2.cut.same %>%
    dplyr::group_by(group_MAPQ_same) %>%
    dplyr::group_split()
j.3.diff.n2.cut.same.mq_same <- j.split[[1]]
j.3.diff.n2.cut.same.mq_diff <- j.split[[2]]
# rm(j.3.diff.n2.cut.same, j.split)


#  Goal: In - *.n1.mq_diff; out - *.n1.mq_diff.dedup --------------------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has both
#+ the best tag.AS score and best mapq score
j.3.diff.n2.cut.same.mq_diff.dedup <-
    j.3.diff.n2.cut.same.mq_diff[
        !duplicated(j.3.diff.n2.cut.same.mq_diff$coordinate), 
    ]

# rm(j.3.diff.n2.cut.same.mq_diff)


#  Goal: In - *.n1.mq_same; out - *.n1.mq_same.dedup --------------------------
#+ !DEDUP
#+ 
#+ Per group, pick one row at random; per group, there are no differences in
#+ tag.AS and mapq scores
set.seed(24)
j.3.diff.n2.cut.same.mq_same.dedup <- j.3.diff.n2.cut.same.mq_same %>%
    group_by(coordinate) %>%
    slice_sample(n = 1)

# rm(j.3.diff.n2.cut.same.mq_same)


#  Goal: In - j.3.diff.n2.cut.diff; out - j.3.diff.n2.cut.diff.dedup --------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has the
#+ best tag.AS score (and thus also the best or same mapq)
j.3.diff.n2.cut.diff.dedup <-
    j.3.diff.n2.cut.diff[
        !duplicated(j.3.diff.n2.cut.diff$coordinate), 
    ]

# rm(j.3.diff.n2.cut.diff)


#  Goal: In - j.3.same; out - j.3.same.mq_same, j.3.same.mq_diff --------------
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has the
#+ best tag.AS score
j.3.same <- j.3.same %>%
    dplyr::group_by(coordinate) %>%
    dplyr::mutate(group_MAPQ_same = dplyr::if_else(
        dplyr::n_distinct(mapq) == 1, "same", "diff"
    )) %>%
    dplyr::mutate_at(vars(group_MAPQ_same), dplyr::funs(factor))
#TODO dplyr::funs() is deprecated; replace it with another function

j.split <- j.3.same %>%
    dplyr::group_by(group_MAPQ_same) %>%
    dplyr::group_split()
j.3.same.mq_same <- j.split[[1]]
j.3.same.mq_diff <- j.split[[2]]

# rm(j.3.same, j.split)


#  Goal: In - j.3.same.mq_diff; out - j.3.same.mq_diff.dedup ------------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has both
#+ the best tag.AS score and best mapq score
j.3.same.mq_diff.dedup <-
    j.3.same.mq_diff[
        !duplicated(j.3.same.mq_diff$coordinate), 
    ]

# rm(j.3.same.mq_diff)


#  Goal: In - j.3.same.mq_same; out - j.3.same.mq_same.dedup ------------------
#+ !DEDUP
#+ 
#+ Per group, pick one row at random; now, per group, there are no differences
#+ in tag.AS and mapq scores
set.seed(24)
j.3.same.mq_same.dedup <- j.3.same.mq_same %>%
    group_by(coordinate) %>%
    slice_sample(n = 1)

# rm(j.3.same.mq_same)


#  Goal: In - j.2; out - j.2.same, j.2.diff -----------------------------------
#+ 
#+ If tag.AS is the same for each member of the group, then "same;" else "diff"
j.2 <- j.2 %>%
    dplyr::group_by(coordinate) %>%
    dplyr::mutate(group_AS_same = dplyr::if_else(
        dplyr::n_distinct(tag.AS) == 1, "same", "diff"
    )) %>%
    dplyr::mutate_at(vars(group_AS_same), dplyr::funs(factor))

j.split <- j.2 %>%
    dplyr::group_by(group_AS_same) %>%
    dplyr::group_split()
j.2.same <- j.split[[1]]
j.2.diff <- j.split[[2]]

# rm(j.2, j.split)


#  Goal: In - j.2.diff; out - j.2.diff.dedup ----------------------------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has the
#+ best tag.AS score
j.2.diff.dedup <- j.2.diff[!duplicated(j.2.diff$coordinate), ]

# rm(j.2.diff)


#  Goal: In - j.2.same; out - j.2.same.mq_same, j.2.same.mq_diff --------------
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has the
#+ best tag.AS score
j.2.same <- j.2.same %>%
    dplyr::group_by(coordinate) %>%
    dplyr::mutate(group_MAPQ_same = dplyr::if_else(
        dplyr::n_distinct(mapq) == 1, "same", "diff"
    )) %>%
    dplyr::mutate_at(vars(group_MAPQ_same), dplyr::funs(factor))

j.split <- j.2.same %>%
    dplyr::group_by(group_MAPQ_same) %>%
    dplyr::group_split()
j.2.same.mq_same <- j.split[[1]]
j.2.same.mq_diff <- j.split[[2]]

# rm(j.2.same, j.split)


#  Goal: In - j.2.same.mq_diff; out - j.2.same.mq_diff.dedup ------------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has both
#+ the best tag.AS score and best mapq score
j.2.same.mq_diff.dedup <-
    j.2.same.mq_diff[!duplicated(j.2.same.mq_diff$coordinate), ]

# rm(j.2.same.mq_diff)


#  Goal: In - j.2.same.mq_same; out - j.2.same.mq_same.dedup ------------------
#+ !DEDUP
#+ 
#+ Per group, pick one row at random; per group, there are no differences in
#+ tag.AS and mapq scores
set.seed(24)
j.2.same.mq_same.dedup <- j.2.same.mq_same %>%
    group_by(coordinate) %>%
    slice_sample(n = 1)

# rm(j.2.same.mq_same)


# -----------------------------------------------------------------------------
dedup.joint.GB <- dplyr::bind_rows(
    j.nondup,
    j.2.diff.dedup,
    j.2.same.mq_diff.dedup,
    j.2.same.mq_same.dedup,
    j.3.diff.n2.cut.diff.dedup,
    j.3.diff.n2.cut.same.mq_diff.dedup,
    j.3.diff.n2.cut.same.mq_same.dedup,
    j.3.diff.n3.dedup,
    j.3.same.mq_diff.dedup,
    j.3.same.mq_same.dedup
) %>% 
    select(
        -duplicated,
        -group_n_members,
        -group_AS_same,
        -group_MAPQ_same,
        -group_AS_n_distinct,
        -group_AS_n_distinct.post_cut
    )

# rm(
#     j.2.diff.dedup,
#     j.2.same.mq_diff.dedup,
#     j.2.same.mq_same.dedup,
#     j.3.diff.n2.cut.diff.dedup,
#     j.3.diff.n2.cut.same.mq_diff.dedup,
#     j.3.diff.n2.cut.same.mq_same.dedup,
#     j.3.diff.n3.dedup,
#     j.3.same.mq_diff.dedup,
#     j.3.same.mq_same.dedup,
#     j.nondup
# )


#  Split and save the jointly deduplicated data sets from GB's analyses -------
split.GB <- dedup.joint.GB %>%
    group_by(assignment_GB) %>%
    group_split()

variable_dedup_joint <- paste0("dedup.joint.", variable)
command <- paste0("<- split.GB[[", c(1:4), "]]")
operation <- makeOperation(variable_dedup_joint, command)
evaluateOperation(operation)

# rm(dedup.joint.GB, split.GB)


#  Order data-set columns by rname, pos ---------------------------------------
command <- paste0(
    "<- ", variable_dedup_joint, "[",
        "order(",
            variable_dedup_joint, "$rname, ", variable_dedup_joint, "$pos), ",
    "]"
)
operation <- makeOperation(variable_dedup_joint, command)
evaluateOperation(operation)


#  Save tibbles as compressed .rds files --------------------------------------
command <- paste0("saveRDS(",
        variable_dedup_joint, ", file = \"", variable_dedup_joint, ".rds\"",
")")
evaluateOperation(command)
#TODO 1/2 Example of what these output files look like and/or some systematic
#TODO 2/2 way to tell they're from this script/step in the pipeline

#  Output is, e.g., "dedup.joint.GB.alt.chrX.rds" in the same working
#+ directory, "segregatedReads.SNPTHRESH1.Q30"

# }

rm(list = ls())
