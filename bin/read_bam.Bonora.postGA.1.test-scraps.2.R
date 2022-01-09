#!/usr/bin/env Rscript

packages <- c(
    "foreach",
    "doParallel",
    "parallel",
    "Rsamtools",
    "tidyverse"
)

for(i in packages){
    suppressPackageStartupMessages(library(i, character.only = TRUE))
}
rm(i, packages)

options(pillar.sigfig = 8, scipen = 10000)

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

loadTibbleFromRDS <- function(variable, file) {
    command <- paste0(variable, " <- loadRDS(file = \"", file, "\")") %>% as.list()
    return(eval(parse(text = command), envir = globalenv()))
}

makeOperation <- function(variable, command) {
    operation <- paste(variable, command) %>% paste(., collapse = "; ")
    return(operation)
}


#  Load tibbles from .rds files -----------------------------------------------
file <- list.files(pattern = "\\.chr1.rds$")
variable <- file %>% gsub("\\.rds$", "", .)
loadTibbleFromRDS(variable = variable, file = file)


#  Append information about the GB assignment ---------------------------------
command <- paste0("<- \"", gsub("\\.chr*(.+)$", "", variable), "\"")
operation <- makeOperation(paste0(variable, "$assignment_GB"), command)
eval(parse(text = operation))


#  Filter out all rows without unique coordinates -----------------------------
#+ 
#+ (Numbers are for data sets comprised of only chr1)
# -------
#  Pre-check: nrow() 
GB.alt.chr1 %>% nrow()  # [1] 2599984
GB.ambig.chr1 %>% nrow()  # [1] 9013879
GB.contra.chr1 %>% nrow()  # [1] 7900
GB.ref.chr1 %>% nrow()  # [1] 2655891

# --
#  Pre-check: Collect duplicate information: Get one of each duplicated element
command <- paste0("<- ", variable, " %>% dplyr::filter(duplicated(.[[\"coordinate\"]]))")
operation <- makeOperation(paste0(variable,".each_dup"), command)
eval(parse(text = operation))
GB.alt.chr1.each_dup %>% nrow()  # [1] 880904
GB.ambig.chr1.each_dup %>% nrow()  # [1] 4090693
GB.contra.chr1.each_dup %>% nrow()  # [1] 2319
GB.ref.chr1.each_dup %>% nrow()  # [1] 901236

# --
#  Pre-check: Collect duplicate information: Get all of the duplicated elements
command <- paste0("<- ", variable, " %>% dplyr::filter(coordinate %in% unique(.[[\"coordinate\"]][duplicated(.[[\"coordinate\"]])]))")
operation <- makeOperation(paste0(variable, ".all_dups"), command)
eval(parse(text = operation))
GB.alt.chr1.all_dups %>% nrow()  # [1] 1344095
GB.ambig.chr1.all_dups %>% nrow()  # [1] 6223797
GB.contra.chr1.all_dups %>% nrow()  # [1] 3554
GB.ref.chr1.all_dups %>% nrow()  # [1] 1374738


#  Deduplicate the tibbles ----------------------------------------------------
command <- paste0("<- ", variable, "[!duplicated(", variable, "$coordinate), ]")
operation <- makeOperation(paste0("dedup.", variable), command)
eval(parse(text = operation))
dedup.GB.alt.chr1 %>% nrow()  # [1] 1719080
dedup.GB.ambig.chr1 %>% nrow()  # [1] 4923186
dedup.GB.contra.chr1 %>% nrow()  # [1] 5581
dedup.GB.ref.chr1 %>% nrow()  # [1] 1754655

((dedup.GB.alt.chr1 %>% nrow()) / (GB.alt.chr1 %>% nrow())) * 100  # [1] 66.11887
((dedup.GB.ambig.chr1 %>% nrow()) / (GB.ambig.chr1 %>% nrow())) * 100  # [1] 54.61784
((dedup.GB.contra.chr1 %>% nrow()) / (GB.contra.chr1 %>% nrow())) * 100  # [1] 70.64557
((dedup.GB.ref.chr1 %>% nrow()) / (GB.ref.chr1 %>% nrow())) * 100  # [1] 66.06653


#  Rename variables for clarity -----------------------------------------------
#+ Remember, "alt" is "CAST", "ref" is "129"
GB.alt.CAST.chr1 <- GB.alt.chr1
GB.ref.129S1.chr1 <- GB.ref.chr1
dedup.GB.alt.CAST.chr1 <- dedup.GB.alt.chr1
dedup.GB.ref.129S1.chr1 <- dedup.GB.ref.chr1

rm(GB.alt.chr1, GB.ref.chr1, dedup.GB.alt.chr1, dedup.GB.ref.chr1)


#  Are duplicates present jointly (i.e., across the four tibbles)? ------------
variable <- c(
    "dedup.GB.alt.CAST.chr1",
    "dedup.GB.ref.129S1.chr1",
    "dedup.GB.ambig.chr1",
    "dedup.GB.contra.chr1"
)
command <- paste0("<- bind_rows(", variable[1], ", ", variable[2], ", ", variable[3], ", ", variable[4], ")")
operation <- makeOperation(paste0("joint.GB.dedup"), command) 
eval(parse(text = operation))

variable <- "joint.GB.dedup"
command <- paste0("<- ", variable, "[order(", variable, "$coordinate, -", variable, "$tag.AS, -", variable, "$mapq), ]")
operation <- makeOperation(variable, command)
eval(parse(text = operation))
joint.GB.dedup %>% nrow()  # [1] 8402502

command <- paste0("<- ", variable, " %>% dplyr::filter(duplicated(.[[\"coordinate\"]]))")
operation <- makeOperation(paste0(variable,".each_dup"), command)
eval(parse(text = operation))
joint.GB.dedup.each_dup %>% nrow()  # [1] 596796

command <- paste0("<- ", variable, " %>% dplyr::filter(coordinate %in% unique(.[[\"coordinate\"]][duplicated(.[[\"coordinate\"]])]))")
operation <- makeOperation(paste0(variable, ".all_dups"), command)
eval(parse(text = operation))
joint.GB.dedup.all_dups %>% nrow()  # [1] 1186389
joint.GB.dedup.all_dups %>%
    dplyr::select(ID, pos, pos_end, tag.AS, mapq, tag.NM, mpos, mpos_end, assignment_GB) %>%
    View()


#  Check if a given element equals the preceding element ----------------------
joint.mint <- joint.GB.dedup.all_dups %>%
    dplyr::select(ID, pos, pos_end, tag.AS, mapq, tag.NM, mpos, mpos_end, assignment_GB, coordinate)

joint.mint <- joint.mint %>% dplyr::group_by(coordinate)  # Groups:   coordinate [589,593]

joint.mint.tally <- joint.mint %>% dplyr::tally(., sort = TRUE)
joint.mint.tally[which(joint.mint.tally$n == 3), ] %>% nrow()  # [1] 7203
joint.mint.tally[which(joint.mint.tally$n == 2), ] %>% nrow()  # [1] 582390
joint.mint.tally[which(joint.mint.tally$n == 1), ] %>% nrow()  # [1] 0

#  Group metadata
joint.mint %>% dplyr::group_keys()
joint.mint %>% dplyr::group_indices()
joint.mint %>% dplyr::group_rows() %>% head(n = 150)
joint.mint %>% dplyr::group_data()  # Gives keys and no. rows
joint.mint %>% dplyr::group_vars()

joint.mint.sample <- joint.mint[1:10000, ]
joint.mint.sample %>% summarise(n = n())  # Another way to get keys and no. rows
joint.mint.sample %>% summarise(id = cur_group_id())
joint.mint.sample %>% summarise(row = cur_group_rows())
joint.mint.sample %>% summarise(data = list(cur_group()))
joint.mint.sample %>% summarise(data = list(cur_data()))  # These seem to be the way to go
joint.mint.sample %>% summarise(data = list(cur_data_all()))  # These seem to be the way to go

#  Working with groups, and trying to create labels for groups
joint.mint.sample %>% mutate(n_members = summarise(n_members = n()))

group.members <- joint.mint %>% summarise(n_members = n()) %>%
    dplyr::arrange(-n_members, .by_group = TRUE) %>%
    dplyr::mutate_at(vars(n_members), dplyr::funs(factor)) %>%
    dplyr::mutate(n_members = recode(n_members, "2" = "two", "3" = "three"))
group.members.sample <- group.members[1:10000, ]

joint.mint.test <- joint.mint
joint.mint.test <- joint.mint.test[1:10000, ]
# joint.mint.test$test <- group.members.sample$coordinate[group.members.sample$coordinate %in% joint.mint$coordinate]
# joint.mint.test$test <- ifelse(
#     group.members.sample$coordinate %in% joint.mint.test$coordinate,
#     group.members.sample$coordinate,
#     NA
# )
joint.mint.test[joint.mint.test$coordinate %in% group.members.sample$coordinate, ]

joint.mint.test$equal <- (joint.mint.test$coordinate == joint.mint.test$test)

#  Check if an element equals the previous element: if no, "0"; if yes, "1"
joint.mint$pos_same <- c(FALSE, diff(joint.mint$pos) == 0) %>% as.integer(as.logical(.))
joint.mint$pos_AS_same <- c(FALSE, diff(joint.mint$pos) == 0 & diff(joint.mint$tag.AS) == 0) %>% as.integer(as.logical(.))
joint.mint$pos_AS_mapq_same <- c(FALSE, diff(joint.mint$pos) == 0 & diff(joint.mint$tag.AS) == 0) %>% as.integer(as.logical(.))
joint.mint$trinary <- paste0(
    joint.mint$pos_same,
    joint.mint$pos_AS_same,
    joint.mint$pos_AS_mapq_same
) %>% as_factor()

joint.mint$trinary %>% head(n = 1)
joint.mint$trinary %>% table()
joint.mint$trinary %>% table() %>% prop.table()

joint.mint %>% dplyr::select(-pos_same, -pos_AS_same, -pos_AS_mapq_same) %>% View()
joint.mint[which(joint.mint$trinary == "111"), ] %>% dplyr::select(-pos_same, -pos_AS_same, -pos_AS_mapq_same) %>% head(n = 1)
joint.mint[which(dplyr::lead(joint.mint$trinary == "111")), ] %>% dplyr::select(-pos_same, -pos_AS_same, -pos_AS_mapq_same) %>% head(n = 1)

a <- joint.mint[which(joint.mint$trinary == "111"), ]
b <- joint.mint[which(dplyr::lag(joint.mint$trinary == "111")), ]
test <- joint.mint[joint.mint$ID %in% a$ID | joint.mint$ID %in% b$ID, ]
test$trinary %>% table()

test.2 <- joint.mint[joint.mint$ID %in% a$ID & joint.mint$ID %in% b$ID, ]
test.2$trinary %>% table()

test.3 <- joint.mint[joint.mint$ID %in% a$ID & joint.mint$ID %notin% b$ID, ]
test.3$trinary %>% table()

test.4 <- joint.mint[joint.mint$ID %notin% a$ID & joint.mint$ID %in% b$ID, ]
test.4$trinary %>% table()

test %>% View()

test %>% arrange(pos, pos_end) %>% View()

test <- joint.mint %>% dplyr::group_by(coordinate)
test %>% View()

# Goal:
# If, for a group of duplicate rows, there are no differences in tag.AS or mapq,
# then pick one row from the group at random; if there are differences, then
# pick the row with the better tag.AS; if tag.AS are equal, then pick the row
# with the better mapq; if both tag.AS and mapq 

#  Clean up
#TODO Determine better spots fot this...
rm(joint.mint)

rm(test, test.2, test.3, test.4)
rm(a, b)

rm(chromosome, file, suffix, variable, variable_assignment)
rm(command, command_pipe, file, i, operation, variable)
