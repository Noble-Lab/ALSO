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

# -----------------------------------------------------------------------------
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

loadTibbleFromRDS <- function(variable, file) {
    command <- paste0(variable, " <- loadRDS(file = \"", file, "\")") %>% as.list()
    return(eval(parse(text = command), envir = globalenv()))
}

makeOperation <- function(variable, command) {
    operation <- paste(variable, command) %>% paste(., collapse = "; ")
    return(operation)
}


#  Load tibbles (already sorted) from .rds files ------------------------------
chromosome <- "chr1"
file <- list.files(pattern = paste0("\\.", chromosome, ".rds$"))
variable <- file %>% gsub("\\.rds$", "", .)
loadTibbleFromRDS(variable = variable, file = file)
#NOTE Remember, "alt" is "CAST", "ref" is "129"


#  Append information about the GB assignment ---------------------------------
command <- paste0("<- \"", gsub("\\.chr*(.+)$", "", variable), "\"")
operation <- makeOperation(paste0(variable, "$assignment_GB"), command)
eval(parse(text = operation))


#  Deduplication of the individual tibbles ------------------------------------
variable_dedup <- paste0("dedup.", variable)
command <- paste0("<- ", variable, "[!duplicated(", variable, "$coordinate), ]")
operation <- makeOperation(variable_dedup, command)
eval(parse(text = operation))


#  Free up memory by removing the non-deduplicated tibbles --------------------
command <- paste0("rm(", variable, ")")
eval(parse(text = command))

rm(variable)


#  Joint deduplication of the tibbles -----------------------------------------
#+ Create a vector of the four tibble variables
variable_joint <- "joint.GB"
command <- paste0("<- dplyr::bind_rows(", variable_dedup[1], ", ", variable_dedup[2], ", ", variable_dedup[3], ", ", variable_dedup[4], ")")
operation <- makeOperation(variable_joint, command) 
eval(parse(text = operation))

#+ Order the joint tibble by coordinate, AS, then mapq
command <- paste0("<- ", variable_joint, "[order(", variable_joint, "$coordinate, -", variable_joint, "$tag.AS, -", variable_joint, "$mapq), ]")
operation <- makeOperation(variable_joint, command)
eval(parse(text = operation))
joint.GB %>% nrow()  # [1] 8402502

#+ Select all of duplicates within the tibble
command <- paste0("<- ", variable_joint, " %>% dplyr::filter(coordinate %in% unique(.[[\"coordinate\"]][duplicated(.[[\"coordinate\"]])]))")
operation <- makeOperation(paste0(variable_joint, ".all_dups"), command)
eval(parse(text = operation))
# joint.GB.all_dups %>%
#     dplyr::select(ID, pos, pos_end, tag.AS, mapq, assignment_GB) %>%
#     View()


#  Goal: In - m; out - m.2, m.3 -----------------------------------------------
#+ 
#+ If, for a group of duplicate rows, there are no differences in tag.AS or
#+ mapq, then pick one row from the group at random; if there are differences,
#+ then pick the row with the better tag.AS; if tag.AS are equal, then pick the
#+ row with the better mapq
#  Work with "minimal tibbles:" "m"
m <- joint.GB.all_dups %>%
    dplyr::select(ID, coordinate, pos, tag.AS, mapq, assignment_GB)

#  Based on coordinate value, identify if a given row is a duplicate 
m$duplicated <- ave(m$coordinate, m$coordinate, FUN = length) > 1L

#  Create categories based on the number of rows comprising a "duplicate group"
m <- m %>% 
    group_by(coordinate) %>% 
    dplyr::mutate(group_n_members = n())  %>%
    dplyr::mutate_at(vars(group_n_members), dplyr::funs(factor)) %>%
    dplyr::mutate(group_n_members = recode(
        group_n_members, "1" = "one", "2" = "two", "3" = "three")
    )

m$group_n_members %>% table()
m$group_n_members %>% table() %>% prop.table()

#  For better readability, remove the "duplicated" column
m <- m %>%
    # group_by(pos) %>%
    # mutate(group_name = pos, .before = group_n_members) %>%
    dplyr::select(-duplicated)

#  Create separate ms for "duplicate groups" comprised of two members and
#+ "duplicate groups" comprised of three members
m.split <- m %>%
    dplyr::group_by(group_n_members) %>%
    dplyr::group_split()
m.2 <- m.split[[1]]
m.3 <- m.split[[2]]
rm(m.split)


#  Goal: In - m.3; out - m.3.same, m.3.diff ----------------------------------- 
#+ 
#+ If tag.AS is the same for each member of the group, then "same;" else "diff"
m.3 <- m.3 %>%
    dplyr::group_by(coordinate) %>%
    dplyr::mutate(group_AS_same = dplyr::if_else(
        dplyr::n_distinct(tag.AS) == 1, "same", "diff"
        )
    ) %>%
    dplyr::mutate_at(vars(group_AS_same), dplyr::funs(factor))

m.3$group_AS_same %>% table()
m.3$group_AS_same %>% table() %>% prop.table()

m.split <- m.3 %>%
    dplyr::group_by(group_AS_same) %>%
    dplyr::group_split()
m.3.same <- m.split[[1]]
m.3.diff <- m.split[[2]]
rm(m.split)


#  Goal: In - m.3.diff; out - m.3.diff.n2, m.3.diff.n3 ------------------------
#+ 
#+ Separate them by how many differences they have; there should only be two
#+ options for differences: n = 2 or n = 3 differences
m.3.diff <- m.3.diff %>%
    group_by(coordinate) %>%
    dplyr::mutate(group_AS_n_distinct = dplyr::if_else(
        dplyr::n_distinct(tag.AS) == 2, "two", "three"
    )
) %>%
    dplyr::mutate_at(vars(group_AS_same), dplyr::funs(factor))

m.3.diff$group_AS_n_distinct %>% table()
m.3.diff$group_AS_n_distinct %>% table() %>% prop.table()

m.split <- m.3.diff %>%
    dplyr::group_by(group_AS_n_distinct) %>%
    dplyr::group_split()
m.3.diff.n2 <- m.split[[2]]
m.3.diff.n3 <- m.split[[1]]
rm(m.split)


#  Goal: In - m.3.diff.n3; out - m.3.diff.n3.dedup ----------------------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has the
#+ best tag.AS score (and thus also the best or same mapq)
m.3.diff.n3.dedup <- m.3.diff.n3[!duplicated(m.3.diff.n3$coordinate), ]

m.3.diff.n3$assignment_GB %>% table()
m.3.diff.n3.dedup$assignment_GB %>% table()

m.3.diff.n3$assignment_GB %>% table() %>% prop.table()
m.3.diff.n3.dedup$assignment_GB %>% table() %>% prop.table()


#  Goal: In - m.3.diff.n2; out - m.3.diff.n2.cut, *.cut.n1, *.cut.n2 ----------
#+ 
#+ Per group, remove the row with the worst (lowest) tag.AS score; then,
#+ categorize group_AS following the cut of 1/3 of the tibble
m.3.diff.n2.cut <- m.3.diff.n2 %>%
    group_by(coordinate) %>%
    slice_max(tag.AS, n = 2, with_ties = FALSE) %>%
    dplyr::mutate(group_AS_n_distinct.post_cut = dplyr::if_else(
        dplyr::n_distinct(tag.AS) == 1, "same", "diff"
    )
) %>%
    dplyr::mutate_at(vars(group_AS_n_distinct.post_cut), dplyr::funs(factor))

m.3.diff.n2.cut$group_AS_n_distinct.post_cut %>% table()
m.3.diff.n2.cut$group_AS_n_distinct.post_cut %>% table() %>% prop.table()

m.split <- m.3.diff.n2.cut %>%
    dplyr::group_by(group_AS_n_distinct.post_cut) %>%
    dplyr::group_split()
m.3.diff.n2.cut.same <- m.split[[1]]
m.3.diff.n2.cut.diff <- m.split[[2]]
rm(m.split)


#  Goal: In - m.3.diff.n2.cut.same; out - *.n1.mq_same, *.n1.mq_diff ------------
#+ 
#+ Select for mapq before picking at random
#TODO Better description
m.3.diff.n2.cut.same <- m.3.diff.n2.cut.same %>%
    dplyr::group_by(coordinate) %>%
    dplyr::mutate(group_MAPQ_same = dplyr::if_else(
        dplyr::n_distinct(mapq) == 1, "same", "diff"
    )
) %>%
    dplyr::mutate_at(vars(group_MAPQ_same), dplyr::funs(factor))

m.3.diff.n2.cut.same$group_MAPQ_same %>% table()
m.3.diff.n2.cut.same$group_MAPQ_same %>% table() %>% prop.table()

m.split <- m.3.diff.n2.cut.same %>%
    dplyr::group_by(group_MAPQ_same) %>%
    dplyr::group_split()
m.3.diff.n2.cut.same.mq_same <- m.split[[1]]
m.3.diff.n2.cut.same.mq_diff <- m.split[[2]]
rm(m.split)


#  Goal: In - *.n1.mq_diff; out - *.n1.mq_diff.dedup --------------------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has both
#+ the best tag.AS score and best mapq score
m.3.diff.n2.cut.same.mq_diff.dedup <-
    m.3.diff.n2.cut.same.mq_diff[!duplicated(m.3.diff.n2.cut.same.mq_diff$coordinate), ]

m.3.diff.n2.cut.same.mq_diff$assignment_GB %>% table()
m.3.diff.n2.cut.same.mq_diff.dedup$assignment_GB %>% table()

m.3.diff.n2.cut.same.mq_diff$assignment_GB %>% table() %>% prop.table()
m.3.diff.n2.cut.same.mq_diff.dedup$assignment_GB %>% table() %>% prop.table()


#  Goal: In - *.n1.mq_same; out - *.n1.mq_same.dedup --------------------------
#+ !DEDUP
#+ 
#+ Per group, pick one row at random; per group, there are no differences in
#+ tag.AS and mapq scores
set.seed(24)
m.3.diff.n2.cut.same.mq_same.dedup <- m.3.diff.n2.cut.same.mq_same %>%
    group_by(coordinate) %>%
    slice_sample(n = 1)

m.3.diff.n2.cut.same.mq_same$assignment_GB %>% table()
m.3.diff.n2.cut.same.mq_same.dedup$assignment_GB %>% table()

m.3.diff.n2.cut.same.mq_same$assignment_GB %>% table() %>% prop.table()
m.3.diff.n2.cut.same.mq_same.dedup$assignment_GB %>% table() %>% prop.table()


#  Goal: In - m.3.diff.n2.cut.diff; out - m.3.diff.n2.cut.diff.dedup --------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has the
#+ best tag.AS score (and thus also the best or same mapq)
m.3.diff.n2.cut.diff.dedup <-
    m.3.diff.n2.cut.diff[!duplicated(m.3.diff.n2.cut.diff$coordinate), ]

m.3.diff.n2.cut.diff$assignment_GB %>% table()
m.3.diff.n2.cut.diff$assignment_GB %>% table() %>% prop.table()

m.3.diff.n2.cut.diff.dedup$assignment_GB %>% table()
m.3.diff.n2.cut.diff.dedup$assignment_GB %>% table() %>% prop.table()


#  Goal: In - m.3.same; out - m.3.same.mq_same, m.3.same.mq_diff --------------
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has the
#+ best tag.AS score
m.3.same <- m.3.same %>%
    dplyr::group_by(coordinate) %>%
    dplyr::mutate(group_MAPQ_same = dplyr::if_else(
        dplyr::n_distinct(mapq) == 1, "same", "diff"
    )
) %>%
    dplyr::mutate_at(vars(group_MAPQ_same), dplyr::funs(factor))

m.3.same$group_MAPQ_same %>% table()
m.3.same$group_MAPQ_same %>% table() %>% prop.table()

m.split <- m.3.same %>%
    dplyr::group_by(group_MAPQ_same) %>%
    dplyr::group_split()
m.3.same.mq_same <- m.split[[1]]
m.3.same.mq_diff <- m.split[[2]]
rm(m.split)


#  Goal: In - m.3.same.mq_diff; out - m.3.same.mq_diff.dedup ------------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has both
#+ the best tag.AS score and best mapq score
m.3.same.mq_diff.dedup <-
    m.3.same.mq_diff[!duplicated(m.3.same.mq_diff$coordinate), ]

m.3.same.mq_diff$assignment_GB %>% table()
m.3.same.mq_diff$assignment_GB %>% table() %>% prop.table()

m.3.same.mq_diff.dedup$assignment_GB %>% table()
m.3.same.mq_diff.dedup$assignment_GB %>% table() %>% prop.table()


#  Goal: In - m.3.same.mq_same; out - m.3.same.mq_same.dedup ------------------
#+ !DEDUP
#+ 
#+ Per group, pick one row at random; per group, there are no differences in
#+ tag.AS and mapq scores
set.seed(24)
m.3.same.mq_same.dedup <- m.3.same.mq_same %>%
    group_by(coordinate) %>%
    slice_sample(n = 1)

m.3.same.mq_same$assignment_GB %>% table()
m.3.same.mq_same$assignment_GB %>% table() %>% prop.table()

m.3.same.mq_same.dedup$assignment_GB %>% table()
m.3.same.mq_same.dedup$assignment_GB %>% table() %>% prop.table()


#  Goal: In - m.2; out - m.2.same, m.2.diff -----------------------------------
#+ 
#+ If tag.AS is the same for each member of the group, then "same;" else "diff"
m.2 <- m.2 %>%
    dplyr::group_by(coordinate) %>%
    dplyr::mutate(group_AS_same = dplyr::if_else(
        dplyr::n_distinct(tag.AS) == 1, "same", "diff"
    )
) %>%
    dplyr::mutate_at(vars(group_AS_same), dplyr::funs(factor))

m.2$group_AS_same %>% table()
m.2$group_AS_same %>% table() %>% prop.table()

m.split <- m.2 %>%
    dplyr::group_by(group_AS_same) %>%
    dplyr::group_split()
m.2.same <- m.split[[1]]
m.2.diff <- m.split[[2]]
rm(m.split)


#  Goal: In - m.2.diff; out - m.2.diff.dedup ----------------------------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has the
#+ best tag.AS score
m.2.diff.dedup <- m.2.diff[!duplicated(m.2.diff$coordinate), ]

m.2.diff$assignment_GB %>% table()
m.2.diff.dedup$assignment_GB %>% table()

m.2.diff$assignment_GB %>% table() %>% prop.table()
m.2.diff.dedup$assignment_GB %>% table() %>% prop.table()


#  Goal: In - m.2.same; out - m.2.same.mq_same, m.2.same.mq_diff --------------
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has the
#+ best tag.AS score
m.2.same <- m.2.same %>%
    dplyr::group_by(coordinate) %>%
    dplyr::mutate(group_MAPQ_same = dplyr::if_else(
        dplyr::n_distinct(mapq) == 1, "same", "diff"
    )
) %>%
    dplyr::mutate_at(vars(group_MAPQ_same), dplyr::funs(factor))

m.2.same$group_MAPQ_same %>% table()
m.2.same$group_MAPQ_same %>% table() %>% prop.table()

m.split <- m.2.same %>%
    dplyr::group_by(group_MAPQ_same) %>%
    dplyr::group_split()
m.2.same.mq_same <- m.split[[1]]
m.2.same.mq_diff <- m.split[[2]]
rm(m.split)


#  Goal: In - m.2.same.mq_diff; out - m.2.same.mq_diff.dedup ------------------
#+ !DEDUP
#+ 
#+ Per group, take the first row, which, given the sorting up to now, has both
#+ the best tag.AS score and best mapq score
m.2.same.mq_diff.dedup <-
    m.2.same.mq_diff[!duplicated(m.2.same.mq_diff$coordinate), ]

m.2.same.mq_diff$assignment_GB %>% table()
m.2.same.mq_diff.dedup$assignment_GB %>% table()

m.2.same.mq_diff$assignment_GB %>% table() %>% prop.table()
m.2.same.mq_diff.dedup$assignment_GB %>% table() %>% prop.table()


#  Goal: In - m.2.same.mq_same; out - m.2.same.mq_same.dedup ------------------
#+ !DEDUP
#+ 
#+ Per group, pick one row at random; per group, there are no differences in
#+ tag.AS and mapq scores
set.seed(24)
m.2.same.mq_same.dedup <- m.2.same.mq_same %>%
    group_by(coordinate) %>%
    slice_sample(n = 1)

m.2.same.mq_same$assignment_GB %>% table()
m.2.same.mq_same$assignment_GB %>% table() %>% prop.table()

m.2.same.mq_same.dedup$assignment_GB %>% table()
m.2.same.mq_same.dedup$assignment_GB %>% table() %>% prop.table()


#  Check: Does everything add up (etc.)? --------------------------------------
nrow(m) == (nrow(m.2) + nrow(m.3))  # [1] TRUE

#  m.2
nrow(m.2) == (nrow(m.2.diff) + nrow(m.2.same))  # [1] TRUE

(nrow(m.2.diff) / 2) == nrow(m.2.diff.dedup)  # [1] TRUE

nrow(m.2.same) == nrow(m.2.same.mq_diff) + nrow(m.2.same.mq_same)  # [1] TRUE

(nrow(m.2.same.mq_diff) / 2) == nrow(m.2.same.mq_diff.dedup)  # [1] TRUE

(nrow(m.2.same.mq_same) / 2) == nrow(m.2.same.mq_same.dedup)  # [1] TRUE

((nrow(m.2.same.mq_diff) + nrow(m.2.same.mq_same)) / 2) == nrow(m.2.same.mq_diff.dedup) + nrow(m.2.same.mq_same.dedup)

#  m.3
nrow(m.3) == (nrow(m.3.diff) + nrow(m.3.same))  # [1] TRUE

nrow(m.3.diff) == (nrow(m.3.diff.n2) + nrow(m.3.diff.n3))  # [1] TRUE

(nrow(m.3.diff.n2) * (2/3)) == nrow(m.3.diff.n2.cut)  # [1] TRUE

nrow(m.3.diff.n2.cut) == (nrow(m.3.diff.n2.cut.same) + nrow(m.3.diff.n2.cut.diff))  # [1] TRUE

nrow(m.3.diff.n2.cut.same) == nrow(m.3.diff.n2.cut.same.mq_diff) + nrow(m.3.diff.n2.cut.same.mq_same)  # [1] TRUE

(nrow(m.3.diff.n2.cut.same.mq_diff) / 2) == nrow(m.3.diff.n2.cut.same.mq_diff.dedup)  # [1] TRUE

(nrow(m.3.diff.n2.cut.same.mq_same) / 2) == nrow(m.3.diff.n2.cut.same.mq_same.dedup)  # [1] TRUE

(nrow(m.3.diff.n2.cut.diff) / 2) == nrow(m.3.diff.n2.cut.diff.dedup)  # [1] TRUE

(nrow(m.3.diff.n3) / 3) == nrow(m.3.diff.n3.dedup)  # [1] TRUE

nrow(m.3.same) == nrow(m.3.same.mq_diff) + nrow(m.3.same.mq_same)  # [1] TRUE

(nrow(m.3.same.mq_diff) / 3) == nrow(m.3.same.mq_diff.dedup)  # [1] TRUE

(nrow(m.3.same.mq_same) / 3) == nrow(m.3.same.mq_same.dedup)  # [1] TRUE


# -----------------------------------------------------------------------------
joint.GB.dedup <- dplyr::bind_rows(
    m.2.diff.dedup,
    m.2.same.mq_diff.dedup,
    m.2.same.mq_same.dedup,
    m.3.diff.n2.cut.diff.dedup,
    m.3.diff.n2.cut.same.mq_diff.dedup,
    m.3.diff.n2.cut.same.mq_same.dedup,
    m.3.diff.n3.dedup,
    m.3.same.mq_diff.dedup,
    m.3.same.mq_same.dedup
)

nrow(joint.GB.all_dups)  # [1] 1186389
nrow(joint.GB.dedup)  # [1] 589593
joint.GB.all_dups[!duplicated(joint.GB.all_dups$coordinate), ] %>% nrow()  # [1] 589593
nrow(joint.GB.dedup) == joint.GB.all_dups[!duplicated(joint.GB.all_dups$coordinate), ] %>% nrow()  # [1] TRUE
#  Nice... :)

joint.GB.dedup[duplicated(joint.GB.dedup$coordinate), ]  # A tibble: 0 x 11


#  Clean up -------------------------------------------------------------------
rm(m, m.2, m.2.diff, m.2.same)
rm(m.2.diff.dedup, dedup.m.2.same)
rm(m.3, m.3.diff, m.3.diff.n2)
rm(m.3.diff.n2.cut, m.3.diff.n2.cut.same, m.3.diff.n2.cut.diff)
rm(dedup.m.3.diff.n2.cut.same, m.3.diff.n2.cut.diff.dedup)
rm(m.3.diff.n2.cut.same.mq_diff.dedup, m.3.diff.n2.cut.same.mq_same.dedup)
rm(m.3.diff.n3, m.3.same)
rm(m.3.diff.n3.dedup, dedup.m.3.same)
rm(m.2.same.mq_diff, m.2.same.mq_same)
rm(m.2.same.mq_diff.dedup, m.2.same.mq_same.dedup)
rm(m.3.diff.n2.cut.same.mq_diff, m.3.diff.n2.cut.same.mq_same)
rm(m.3.same.mq_diff, m.3.same.mq_same)
rm(m.3.same.mq_diff.dedup, m.3.same.mq_same.dedup)
