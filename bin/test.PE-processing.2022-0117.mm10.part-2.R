#!/usr/bin/env Rscript

library(stringr)
library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)


#  Set up work directory (locations TB∆) --------------------------------------
directory_user <- "/Users/kalavattam"
directory_base <- "Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do"
directory_work <- "get_unique_fragments/Bonora"

setwd(
    paste0(directory_user, "/", directory_base, "/", directory_work)
)
rm(directory_user, directory_base, directory_work)


#  Having run part 1 first, load part-1 .Rdata into environment ---------------
load("test.PE-processing.2022-0117.mm10.part-1.Rdata")
#  Doing so loads in appropriate functions, variables, etc.


#  Script name ----------------------------------------------------------------
script <- "test.PE-processing.2022-0117.mm10.part-2.R"


#  Add row numbers to sub.* tibbles -------------------------------------------
operation <- paste0(
    variable, "$row_n", " <- ", "seq.int(nrow(", variable, "))"
)
evaluateOperation(operation)

command <- paste0(
    "<- ", variable, " %>% ", "dplyr::relocate(row_n, .before = flag)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Check on mate pairing; generate vectors and tables for mate pairing --------
command <- paste0(
    "<- testMatesPaired(pos = ", variable, "$qpos, mpos = ", variable, "$qmpos)"
)  # VMP: "vector mates paired"
operation <- makeOperation(paste0("VMP.", variable), command)
evaluateOperation(operation)

#  Check: NAs in VMP vector?
command <- paste0("<- ", "is.na(", "VMP.", variable, ")")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)

command <- paste0("<- isNA.VMP.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)  # Yes, 1 NA

command <- paste0("rm(", "isNA.VMP.", variable, ")")
evaluateOperation(command)

#  Convert NAs
#  For "ambiguous" and "mated," assign FALSE to NAs
operation <- paste0(
    "VMP.", variable[!str_detect(variable, "unmated")], "[",
        "is.na(", "VMP.", variable[!str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  For "unmated," assign FALSE to NAs
operation <- paste0(
    "VMP.", variable[str_detect(variable, "unmated")], "[",
        "is.na(", "VMP.", variable[str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  Check: NAs in VMP vector? Should be gone now
command <- paste0("<- ", "is.na(", "VMP.", variable, ")")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)

command <- paste0("<- isNA.VMP.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)  # No NAs

command <- paste0("rm(", "isNA.VMP.", variable, ")")
evaluateOperation(command)

#  Save vectors as tibbles (`TRUE`, `FALSE` columns)
command <- paste0(
    "<- VMP.", variable, " %>% ",
        "table() %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("tbl.VMP.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles"
vector_string_tibbles <- paste0("tbl.VMP.", variable)
z.check.sub.1.VMP <- createTibbleFromTibbles(vector_string_tibbles)
z.check.sub.1.VMP$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.sub.1.VMP <- z.check.sub.1.VMP %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0("rm(", "tbl.VMP.", variable, ")")
evaluateOperation(command)
rm(vector_string_tibbles)


#  Are `TRUE`, `FALSE` *not* alternating? (They should be) --------------------

#  Generate vectors and tables
command <- paste0(
    "<- testLogicalAlternating(logical = ", "VMP.", variable, ")"
)  # VA: "vector alternating"
operation <- makeOperation(paste0("VA.", variable), command)
evaluateOperation(operation)

#  Check: NAs in VA vector?
command <- paste0("<- ", "is.na(", "VA.", variable, ")")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)

command <- paste0("<- isNA.VA.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)  # Yes, 1 NA

command <- paste0("rm(", "isNA.VA.", variable, ")")
evaluateOperation(command)

#  Convert NAs
#  For "ambiguous" and "mated," assign FALSE to NAs
operation <- paste0(
    "VA.", variable[!str_detect(variable, "unmated")], "[",
        "is.na(", "VA.", variable[!str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  For "unmated," assign TRUE to NAs
operation <- paste0(
    "VA.", variable[str_detect(variable, "unmated")], "[",
        "is.na(", "VA.", variable[str_detect(variable, "unmated")], ")",
    "]", " <- ", "TRUE"
)
evaluateOperation(operation)

#  Check: NAs in VA vector? Should be gone now
command <- paste0("<- ", "is.na(", "VA.", variable, ")")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)

command <- paste0("<- isNA.VA.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)  # No NAs

command <- paste0("rm(", "isNA.VA.", variable, ")")
evaluateOperation(command)

#  Save vectors as tibbles (`TRUE`, `FALSE` columns)
command <- paste0(
    "<- VA.", variable, " %>% ",
        "table() %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("tbl.VA.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles"
vector_string_tibbles <- paste0("tbl.VA.", variable)
z.check.sub.2.VA <- createTibbleFromTibbles(vector_string_tibbles)
z.check.sub.2.VA$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.sub.2.VA <- z.check.sub.2.VA %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0("rm(", "tbl.VA.", variable, ")")
evaluateOperation(command)
rm(vector_string_tibbles)


#  Create tbl columns denoting lack of alternation, i.e., $discrepancy --------
command <- paste0("<- ", "VA.", variable)
operation <- makeOperation(paste0(variable, "$discrepancy"), command)
evaluateOperation(operation)

#  In $discrepancy, NA should be FALSE for "mated" and "ambiguous," TRUE for
#+ "unmated;" check via tables
command <- paste0(
    "<- ", "is.na(", variable, "$discrepancy", ")", " %>% ", "table()" 
)
operation <- makeOperation(paste0("isNA.", variable), command)
evaluateOperation(operation)  # No NAs present

#  Remove the temporary tables
command <- paste0("rm(", "isNA.", variable, ")")
evaluateOperation(command)

# #  Convert NAs
#  For "ambiguous" and "mated," assign FALSE to NAs
# operation <- paste0(
#     variable[!str_detect(variable, "unmated")], "$discrepancy[",
#         "is.na(",
#             variable[!str_detect(variable, "unmated")],
#             "$discrepancy",
#         ")",
#     "]", " <- ", "FALSE"
# )  #TODO May delete this b/c no NAs
# evaluateOperation(operation)
#
# #  For "unmated," assign TRUE to NAs
# operation <- paste0(
#     variable[str_detect(variable, "unmated")], "$discrepancy[",
#         "is.na(",
#             variable[str_detect(variable, "unmated")],
#             "$discrepancy",
#         ")",
#     "]", " <- ", "TRUE"
# )  #TODO May delete this b/c no NAs
# evaluateOperation(operation)
#
# #  Check to ensure that the NAs are gone now
# command <- paste0(
#     "<- ", "is.na(", variable, "$discrepancy", ")", " %>% ", "table()" 
# )
# operation <- makeOperation(paste0("isNA.", variable), command)
# evaluateOperation(operation)  # No NAs
#
# #  Remove the temporary tables
# command <- paste0("rm(", "isNA.", variable, ")")
# evaluateOperation(command)


#  Make tibble indices (numeric vectors) for discrepancies ± 1 ----------------
command <- paste0(
    "<- createVectorPlusMinus(",
        "vector = which(", variable, "$discrepancy == TRUE)",
    ")"
)  # VD: "vector discrepancy"
operation <- makeOperation(paste0("VD.", variable), command)
evaluateOperation(operation)


#  Subset tibbles for only discrepancies ± one --------------------------------
command <- paste0(
    "<- ", variable, "[", "VD.", variable, ", ]"
)
operation <- makeOperation(paste0("discrepancy.", variable), command)
evaluateOperation(operation)

#  Sort the "discrepancy tibbles"
command <- paste0(
    "<- ", "discrepancy.", variable, "[",
        "order(",
            "discrepancy.", variable, "$isize_abs, ",
            "discrepancy.", variable, "$qpos, ",
            "discrepancy.", variable, "$qmpos",
        "), ",
    "]"
)
operation <- makeOperation(paste0("discrepancy.", variable), command)
evaluateOperation(operation)


#  From the sub tibbles, filter out duplicate $qpos, $qmpos -------------------

#QUESTION Is this step handled properly?
command <- paste0(
    "<- ", variable, "[",
        "!duplicated(", variable, "$qpos) & ",
        "!duplicated(", variable, "$qmpos), ",
    "]"
)
operation <- makeOperation(
    paste0("filter.", variable) %>% stringr::str_remove("sub."),
    command
)
evaluateOperation(operation)

#  Remove unneeded variables
command <- paste0(
    "rm(",
        "VA.", variable, ", ",
        "VD.", variable, ", ",
        "VMP.", variable,
    ")"
)
evaluateOperation(command)


#  Are mates paired after filtering out duplicate $qpos, $qmpos? --------------
variable <- paste0("filter.", variable) %>% stringr::str_remove("sub.")

#  Add row numbers to filter.* tibbles
operation <- paste0(
    variable, "$row_n", " <- ", "seq.int(nrow(", variable, "))"
)
evaluateOperation(operation)

#  Generate vectors and tibbles 
command <- paste0(
    "<- testMatesPaired(",
        "pos = ", variable, "$qpos, ",
        "mpos = ", variable, "$qmpos",
    ")"
)
operation <- makeOperation(paste0("VMP.", variable), command)
evaluateOperation(operation)

#  Check: NAs in VMP vector?
command <- paste0("<- ", "is.na(", "VMP.", variable, ")")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)

command <- paste0("<- isNA.VMP.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)  # Yes, 1 NA

command <- paste0("rm(", "isNA.VMP.", variable, ")")
evaluateOperation(command)

#  Convert NAs
#  For "ambiguous" and "mated," assign FALSE to NAss
operation <- paste0(
    "VMP.", variable[!str_detect(variable, "unmated")], "[",
        "is.na(", "VMP.", variable[!str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  For "unmated," assign FALSE to NAs
operation <- paste0(
    "VMP.", variable[str_detect(variable, "unmated")], "[",
        "is.na(", "VMP.", variable[str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  Check: NAs in VMP vector?
command <- paste0("<- ", "is.na(", "VMP.", variable, ")")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)

command <- paste0("<- isNA.VMP.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)  # No NAs

command <- paste0("rm(", "isNA.VMP.", variable, ")")
evaluateOperation(command)

command <- paste0(
    "<- VMP.", variable, " %>% ",
        "table() %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("tbl.VMP.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles"
vector_string_tibbles <- paste0("tbl.VMP.", variable)
z.check.filter.1.VMP <- createTibbleFromTibbles(vector_string_tibbles)
z.check.filter.1.VMP$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.filter.1.VMP <- z.check.filter.1.VMP %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0("rm(", "tbl.VMP.", variable, ")")
evaluateOperation(command)
rm(vector_string_tibbles)


#  Are `TRUE`, `FALSE` *not* alternating? (They should be) --------------------

#  That is, `TRUE` means rows *are not* alternating; `FALSE` means they *are*
#+ alternating

#  Create vectors and tibbles
command <- paste0(
    "<- testLogicalAlternating(logical = ", "VMP.", variable, ")"
)  # VA: "vector alternating"
operation <- makeOperation(paste0("VA.", variable), command)
evaluateOperation(operation)

#  Check: NAs in VA vector?
command <- paste0("<- ", "is.na(", "VA.", variable, ")")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)

command <- paste0("<- isNA.VA.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)  # Yes, 1 NA

command <- paste0("rm(", "isNA.VA.", variable, ")")
evaluateOperation(command)

#  Convert NAs
#  For "ambiguous" and "mated," assign FALSE to NAs
operation <- paste0(
    "VA.", variable[!str_detect(variable, "unmated")], "[",
        "is.na(", "VA.", variable[!str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  For "unmated," assign TRUE to NAs
operation <- paste0(
    "VA.", variable[str_detect(variable, "unmated")], "[",
        "is.na(", "VA.", variable[str_detect(variable, "unmated")], ")",
    "]", " <- ", "TRUE"
)
evaluateOperation(operation)

#  Check: NAs in VA vector?
command <- paste0("<- ", "is.na(", "VA.", variable, ")")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)

command <- paste0("<- isNA.VA.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)  # No NAs

command <- paste0("rm(", "isNA.VA.", variable, ")")
evaluateOperation(command)

command <- paste0(
    "<- VA.", variable, " %>% ",
        "table() %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("tbl.VA.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles"
vector_string_tibbles <- paste0("tbl.VA.", variable)
z.check.filter.2.VA <- createTibbleFromTibbles(vector_string_tibbles)
z.check.filter.2.VA$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.filter.2.VA <- z.check.filter.2.VA %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0("rm(", "tbl.VA.", variable, ")")
evaluateOperation(command)
rm(vector_string_tibbles)


#  Create tbl columns denoting lack of alternating, i.e., "discrepancies" -----
command <- paste0(
    "<- ", "VA.", variable
)
operation <- makeOperation(paste0(variable, "$discrepancy"), command)
evaluateOperation(operation)

#  In $discrepancy, should be no NAs? Check via tables
command <- paste0(
    "<- ", "is.na(", variable, "$discrepancy", ")", " %>% ", "table()" 
)
operation <- makeOperation(paste0("isNA.", variable), command)
evaluateOperation(operation)  # No NAs present

#  Remove the temporary tables
command <- paste0("rm(", "isNA.", variable, ")")
evaluateOperation(command)


#  Make numeric vectors (i.e., tibble indices) for discrepancies ± 1 ----------
command <- paste0(
    "<- createVectorPlusMinus(",
        "vector = which(", variable, "$discrepancy == TRUE)",
    ")"
)  # VD: "vector discrepancy"
operation <- makeOperation(paste0("VD.", variable), command)
evaluateOperation(operation)


#  Subset tibbles for only discrepancies ± one --------------------------------
command <- paste0(
    "<- ", variable, "[", "VD.", variable, ", ]"
)
operation <- makeOperation(paste0("discrepancy.", variable), command)
evaluateOperation(operation)

# #  Sort the "discrepancy tibbles"
# command <- paste0(
#     "<- ", "discrepancy.", variable, "[",
#         "order(",
#             "discrepancy.", variable, "$pos, ",
#             "discrepancy.", variable, "$mpos",
#         "), ",
#     "]"
# )
# operation <- makeOperation(paste0("discrepancy.", variable), command)
# evaluateOperation(operation)

#  Sort the "discrepancy tibbles"
command <- paste0(
    "<- ", "discrepancy.", variable, "[",
        "order(",
            "discrepancy.", variable, "$isize_abs, ",
            "discrepancy.", variable, "$qname, ",
            "discrepancy.", variable, "$pos",
        "), ",
    "]"
)
operation <- makeOperation(paste0("discrepancy.", variable), command)
evaluateOperation(operation)


# #TODO Some automated way to do the following
# #  Manually correct filter.a.ambiguous.dedup.CAST -----------------------------
# filter.a.ambiguous.dedup.CAST[paste0("584", c(19, 20, 21, 22, 23, 24)), ] <-
#     discrepancy.filter.a.ambiguous.dedup.CAST
# 
# #  Check
# filter.a.ambiguous.dedup.CAST[c(58419:58429), ]
# filter.a.ambiguous.dedup.CAST$row_n[c(58419:58429)]
# # [1] 58419 58420 58421 58424 58422 58423 58425 58426 58427 58428 58429
# 
# #  Re-label $row_n
# filter.a.ambiguous.dedup.CAST$row_n <- seq.int(
#     nrow(filter.a.ambiguous.dedup.CAST)
# )
# filter.a.ambiguous.dedup.CAST[c(58419:58429), ]
# filter.a.ambiguous.dedup.CAST$row_n[c(58419:58429)]
# # [1] 58419 58420 58421 58422 58423 58424 58425 58426 58427 58428 58429


# #  Correct filter.a.mated.dedup.129S1 -----------------------------------------
# desired_order <- paste0("343", c(79, 80, 81, 83, 82, 84))
# discrepancy.filter.a.mated.dedup.129S1.adjust <-
#     discrepancy.filter.a.mated.dedup.129S1 %>%
#     dplyr::mutate(row_n = factor(row_n, levels = desired_order)) %>%
#     dplyr::arrange(row_n)
# discrepancy.filter.a.mated.dedup.129S1.adjust$row_n
# 
# discrepancy.filter.a.mated.dedup.129S1.adjust$row_n <-
#     discrepancy.filter.a.mated.dedup.129S1.adjust$row_n %>%
#     as.character() %>%
#     as.integer()
# 
# filter.a.mated.dedup.129S1[discrepancy.filter.a.mated.dedup.129S1$row_n, ] <-
#     discrepancy.filter.a.mated.dedup.129S1.adjust
# 
# discrepancy.filter.a.mated.dedup.129S1 <-
#     discrepancy.filter.a.mated.dedup.129S1.adjust
# 
# #  Check
# filter.a.mated.dedup.129S1[c(34379:34389), ]
# filter.a.mated.dedup.129S1$row_n[c(34379:34389)]
# # [1] 34379 34380 34381 34383 34382 34384 34385 34386 34387 34388 34389
# 
# #  Re-label $row_n
# filter.a.mated.dedup.129S1$row_n <- seq.int(
#     nrow(filter.a.mated.dedup.129S1)
# )
# 
# #  Check
# filter.a.mated.dedup.129S1[c(34379:34389), ]
# filter.a.mated.dedup.129S1$row_n[c(34379:34389)]
# # [1] 34379 34380 34381 34382 34383 34384 34385 34386 34387 34388 34389
# 
# #  Remove unneeded variables
# rm(desired_order, discrepancy.filter.a.mated.dedup.129S1.adjust)


#  Clean up unneeded variables before setting up new ones ---------------------
command <- paste0(
    "rm(",
        "VMP.", variable, ", ",
        "VA.", variable, ", ",
        "VD.", variable,
    ")"
)
evaluateOperation(command)


#  Set up new variables -------------------------------------------------------
command <- paste0("<- ", variable)
operation <- makeOperation(
    paste0("adjust.", stringr::str_remove(variable, "filter.")),
    command
)
evaluateOperation(operation)

variable <- paste0("adjust.", stringr::str_remove(variable, "filter."))


#  Generate vectors and tables for mate-pairing readouts ----------------------
command <- paste0(
    "<- testMatesPaired(pos = ", variable, "$qpos, mpos = ", variable, "$qmpos)"
)  # VMP: "vector mates paired"
operation <- makeOperation(paste0("VMP.", variable), command)
evaluateOperation(operation)

#  Check: NAs in VMP vector?
command <- paste0("<- ", "is.na(", "VMP.", variable, ")")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)

command <- paste0("<- isNA.VMP.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)  # Yes, 1 NA

command <- paste0("rm(", "isNA.VMP.", variable, ")")
evaluateOperation(command)

#  Convert NAs
#  For "ambiguous" and "mated," assign FALSE to NAs
operation <- paste0(
    "VMP.", variable[!str_detect(variable, "unmated")], "[",
        "is.na(", "VMP.", variable[!str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  For "unmated," assign FALSE to NAs
operation <- paste0(
    "VMP.", variable[str_detect(variable, "unmated")], "[",
        "is.na(", "VMP.", variable[str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  Check: NAs in VMP vector? Should be gone now
command <- paste0("<- ", "is.na(", "VMP.", variable, ")")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)

command <- paste0("<- isNA.VMP.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VMP.", variable), command)
evaluateOperation(operation)  # No NAs

command <- paste0("rm(", "isNA.VMP.", variable, ")")
evaluateOperation(command)

#  Save vectors as tibbles (`TRUE`, `FALSE` columns)
command <- paste0(
    "<- VMP.", variable, " %>% ",
        "table() %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("tbl.VMP.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles"
vector_string_tibbles <- paste0("tbl.VMP.", variable)
z.check.adjust.1.VMP <- createTibbleFromTibbles(vector_string_tibbles)
z.check.adjust.1.VMP$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.adjust.1.VMP <- z.check.adjust.1.VMP %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0("rm(", "tbl.VMP.", variable, ")")
evaluateOperation(command)
rm(vector_string_tibbles)


#  Are `TRUE`, `FALSE` *not* alternating? (They should be) --------------------

#  Generate vectors and tables
command <- paste0(
    "<- testLogicalAlternating(logical = ", "VMP.", variable, ")"
)  # VA: "vector alternating"
operation <- makeOperation(paste0("VA.", variable), command)
evaluateOperation(operation)

#  Check: NAs in VA vector?
command <- paste0("<- ", "is.na(", "VA.", variable, ")")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)

command <- paste0("<- isNA.VA.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)  # Yes, 1 NA

command <- paste0("rm(", "isNA.VA.", variable, ")")
evaluateOperation(command)

#  Convert NAs
#  For "ambiguous" and "mated," assign FALSE to NAs
operation <- paste0(
    "VA.", variable[!str_detect(variable, "unmated")], "[",
        "is.na(", "VA.", variable[!str_detect(variable, "unmated")], ")",
    "]", " <- ", "FALSE"
)
evaluateOperation(operation)

#  For "unmated," assign TRUE to NAs
operation <- paste0(
    "VA.", variable[str_detect(variable, "unmated")], "[",
        "is.na(", "VA.", variable[str_detect(variable, "unmated")], ")",
    "]", " <- ", "TRUE"
)
evaluateOperation(operation)

#  Check: NAs in VA vector? Should be gone now
command <- paste0("<- ", "is.na(", "VA.", variable, ")")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)

command <- paste0("<- isNA.VA.", variable, " %>% ", "table()")
operation <- makeOperation(paste0("isNA.VA.", variable), command)
evaluateOperation(operation)  # No NAs

command <- paste0("rm(", "isNA.VA.", variable, ")")
evaluateOperation(command)

#  Save vectors as tibbles (`TRUE`, `FALSE` columns)
command <- paste0(
    "<- VA.", variable, " %>% ",
        "table() %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("tbl.VA.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles"
vector_string_tibbles <- paste0("tbl.VA.", variable)
z.check.adjust.2.VA <- createTibbleFromTibbles(vector_string_tibbles)
z.check.adjust.2.VA$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.adjust.2.VA <- z.check.adjust.2.VA %>% dplyr::arrange(rownames)


#  Remove unneeded variables --------------------------------------------------
command <- paste0("rm(", "tbl.VA.", variable, ")")

variable_non_adjust <- variable %>% stringr::str_remove("adjust.")
variable_non_a <- variable_non_adjust %>% stringr::str_remove("a.")
command <- paste0(
    "rm(",
        "VMP.", variable, ", ",
        "VA.", variable, ", ",
        "tbl.VA.", variable, ", ",
        "same.", variable_non_adjust, ", ",
        "sub.", variable_non_adjust, ", ",
        "discrepancy.sub.", variable_non_adjust, ", ",
        "discrepancy.filter.", variable_non_adjust, ", ",
        "filter.", variable_non_adjust, ", ",
        variable_non_adjust, ", ",
        "z.bak.", variable_non_a,
    ")"
)
evaluateOperation(command)

rm(
    vector_string_tibbles,
    variable_non_a,
    variable_non_adjust,
    z.bak.variable,
    z.init.dedup.129S1,
    z.init.dedup.CAST,
    z.init.dedup.mm10
)

command <- paste0("rm(",
                      paste0("z.check.", c(
                          "0.change", "0.equal",
                          "1.MIP", "1.PIM",
                          "2.PeqM", "2.PgtM", "2.PltM",
                          "sub.0.MIP", "sub.0.PIM", "sub.1.VMP", "sub.2.VA",
                          "filter.1.VMP", "filter.2.VA",
                          "adjust.1.VMP", "adjust.2.VA"
                      )),
                  ")")
evaluateOperation(command)


#  Generate $groupid_2 for munged "mated" and "ambiguous" tibbles -------------

#  Work with only "mated" and "ambiguous" tibbles, not "unmated" tibbles
variable <- variable[!str_detect(variable, "unmated")]

for (i in 1:length(variable)) {
    df <- eval(parse(text = variable[i]))
    odd <- seq(1, nrow(df), 2)
    even <- (seq(1, nrow(df), 2)) + 1
    
    odd <- df[odd, ]
    even <- df[even, ]
    
    odd$tibble <- "1"
    even$tibble <- "2"
    
    odd <- odd %>% dplyr::relocate(tibble, .before = row_n)
    even <- even %>% dplyr::relocate(tibble, .before = row_n)
    
    command <- paste0(
        "<- odd %>%",
            "dplyr::mutate(groupid_2 = row_number()) %>%",
            "dplyr::bind_rows(even %>% mutate(groupid_2 = row_number())) %>%",
            "dplyr::arrange(groupid_2, tibble) %>%",
            "dplyr::relocate(groupid_2, .after = groupid) %>%",
            "dplyr::select(-tibble)"
    )
    operation <- makeOperation(variable[i], command)
    evaluateOperation(operation)
    
    n_occur <- data.frame(table(eval(parse(text = variable[i]))$groupid_2))
    print(n_occur[n_occur$Freq > 2, ])
}
rm(i)

rm(df, odd, odd.seq, even, n_occur)


#  Combine "mated" and "ambiguous" tibbles ------------------------------------
bindRows <- function(tbl) {
    dplyr::bind_rows(
        evaluateOperation(tbl[1]),
        evaluateOperation(tbl[2]),
    )
}

tbl <- stringr::str_subset(variable, "mm10")
tbl.mm10 <- bindRows(tbl)

#  Remove unneeded data from environment
rm(tbl, bindRows)


#  Generate groupid.* variables for combination mated/ambiguous tibbles -------
variable_tbl <- paste0("tbl.", "mm10")

for (i in 1:length(variable_tbl)) {
    #  Assign tibble of interest to variable 'df'
    df <- eval(parse(text = variable_tbl[i]))
    
    #  Create tibbles from odd and even rows
    odd.seq <- seq(1, nrow(df), 2)
    even.seq <- (seq(1, nrow(df), 2)) + 1
    
    odd <- df[odd.seq, ]
    even <- df[even.seq, ]
    
    #  Odd tibble is tibble #1, even tibble is tibble #2
    odd$tibble <- "1"
    even$tibble <- "2"
    
    odd <- odd %>% dplyr::relocate(tibble, .before = row_n)
    even <- even %>% dplyr::relocate(tibble, .before = row_n)
    
    #  Interleave the odd/even rows, then arrange them by group ID and tibble
    #+ number
    command <- paste0(
        "<- odd %>%",
            "dplyr::mutate(groupid_3 = row_number()) %>% ",
            "dplyr::bind_rows(even %>% mutate(groupid_3 = row_number())) %>% ",
            "dplyr::arrange(groupid_3, tibble) %>% ",
            "dplyr::relocate(groupid_3, .after = groupid_2) %>% ",
            "dplyr::select(-tibble)"
    )
    operation <- makeOperation(variable_tbl[i], command)
    evaluateOperation(operation)
    
    #  Check to make sure that there are no more than two entries per group ID
    n_occur <- data.frame(table(eval(parse(text = variable_tbl[i]))$groupid_3))
    print(n_occur[n_occur$Freq > 2, ])
    rm(n_occur)
    
    #  Sort the tibble of interest by pos while maintaining proper mate pairs
    df <- eval(parse(text = variable_tbl[i]))
    
    odd <- df[odd.seq, ]
    even <- df[even.seq, ]
    
    #  Create a tibble with two mates per row (joined by "groupid_3"); then,
    #+ sort by pos.x (the odd-row pos value)
    df <- full_join(odd, even, by = "groupid_3") %>%
        dplyr::arrange(pos.x) %>%  # Sort by pos.x
        dplyr::rename(groupid_3.x = groupid_3)  # new_name = old_name
    df$groupid_3.y <- df$groupid_3.x
    df <- df %>% dplyr::relocate(groupid_3.y, .after = groupid_2.y)
    
    #  Rename column names to denote odd/even status
    colnames(df) <- str_replace_all(colnames(df), "\\.x", "\\.odd")
    colnames(df) <- str_replace_all(colnames(df), "\\.y", "\\.even")
    
    #  Split the tibble by odd or even status
    odd <- df[stringr::str_subset(colnames(df), "\\.odd")]
    even <- df[stringr::str_subset(colnames(df), "\\.even")]
    
    #  Strip suffixes from column names
    colnames(odd) <- str_replace_all(colnames(odd), "\\.odd", "")
    colnames(even) <- str_replace_all(colnames(even), "\\.even", "")
    
    #  Odd tibble is tibble #1, even tibble is tibble #2
    odd$tibble <- "1"
    even$tibble <- "2"
    
    #  Again, interleave the odd/even rows, then arrange them by group ID and
    #+ tibble number
    command <- paste0(
        "<- odd %>%",
            "dplyr::mutate(groupid_4 = row_number()) %>% ",
            "dplyr::bind_rows(even %>% mutate(groupid_4 = row_number())) %>% ",
            "dplyr::arrange(groupid_4, tibble) %>% ",
            "dplyr::relocate(groupid_4, .after = groupid_3) %>% ",
            "dplyr::select(-tibble)"
    )
    operation <- makeOperation(variable_tbl[i], command)
    evaluateOperation(operation)
    
    #  Again, check to make sure that there are no more than two entries per
    #+ group ID
    n_occur <- data.frame(table(eval(parse(text = variable_tbl[i]))$groupid_4))
    print(n_occur[n_occur$Freq > 2, ])
    rm(n_occur)
}

rm(df, i, even, even.seq, odd, odd.seq)

#  Remove unneeded columns from the tbl.* tibbles
command <- paste0(
    "<- ", variable_tbl, " %>% ",
        "dplyr::select(-groupid, -groupid_2, -groupid_3, -row_n)", " %>% ",
        "dplyr::rename(groupid = groupid_4)"
)
operation <- makeOperation(variable_tbl, command)
evaluateOperation(operation)


# #  Merging rows ---------------------------------------------------------------
# 
# #  See this URL:
# #  stackoverflow.com/questions/21003237/read-odd-numbered-rows-from-csv-file
# 
# #  Initialize one-line ("uniline") variables
# variable_uniline <- paste0("uniline.", str_remove(variable, "adjust."))
# mapply(
#     assign, variable_uniline, "", MoreArgs = list(envir = parent.frame())
# )
# 
# #  Concatenate even row to preceding odd row in  adjust.* variables
# for (i in 1:length(variable)) {
#     df <- eval(parse(text = variable[i]))
#     odd <- seq(1, nrow(df), 2)
#     even <- (seq(1, nrow(df), 2)) + 1
#     
#     odd <- df[odd, ]
#     even <- df[even, ]
#     
#     df <- full_join(odd, even, by = "groupid_2")
#     colnames(df) <- str_replace_all(colnames(df), "\\.x", "\\.odd")
#     colnames(df) <- str_replace_all(colnames(df), "\\.y", "\\.even")
#     
#     command <- paste0("<- ", "df")
#     operation <- makeOperation(variable_uniline[i], command)
#     evaluateOperation(operation)
# }
# 
# #  Remove unneeded variables
# rm(df, odd, even)
# 
# #  Check flag statuses for adjust.*
# for (i in 1:length(variable)) {
#     print(variable[i])
#     operation <- paste0(
#         variable[i], "$flag %>% table()"
#     )
#     print(evaluateOperation(operation))
# }
# 
# #  $flag.odd for uniline.*: Should be only flags 99 and 163
# #QUESTION Is it a problem if not?
# for (i in 1:length(variable_uniline)) {
#     print(variable_uniline[i])
#     operation <- paste0(
#         variable_uniline[i], "$flag.odd %>% table()"
#     )
#     print(evaluateOperation(operation))
# }
# 
# #  $flag.even for uniline.*: Should be only flags 83 and 147
# #QUESTION Is it a problem if not?
# for (i in 1:length(variable_uniline)) {
#     print(variable_uniline[i])
#     operation <- paste0(
#         variable_uniline[i], "$flag.even %>% table()"
#     )
#     print(evaluateOperation(operation))
# }
# rm(i)


# #  Final munging --------------------------------------------------------------
# 
# #  #NOTE
# #+ Whether merged or not (probably easier if I don't merge), I need to see if
# #+ there's a SNP present in at least one of each group's $seq; if so, then as-
# #+ sign that group appropriately
# #+ 
# #+ Ensure that the environment is clean and save an .Rdata image for a new
# #+ subsequent script that will determine if SNPs are/aren't present per
# #+ 
# #+ Also, make sure there are separate .odd/.even $groupid_2 columns
# 
# command <- paste0(
#     "<- ", variable_uniline, " %>% ",
#         "dplyr::rename(groupid_2.odd = groupid_2)"
# )
# operation <- makeOperation(variable_uniline, command)
# evaluateOperation(operation)
# 
# command <- paste0("<- ", variable_uniline, "$groupid_2.odd")
# operation <- makeOperation(
#     paste0(variable_uniline, "$groupid_2.even"), command
# )
# evaluateOperation(operation)
# 
# command <- paste0(
#     "<- ", variable_uniline, " %>% ",
#         "dplyr::relocate(groupid_2.even, .after = groupid.even)"
# )
# operation <- makeOperation(variable_uniline, command)
# evaluateOperation(operation)


#  Save environment image -----------------------------------------------------
save.image(stringr::str_remove(script, ".R") %>% paste0(., ".Rdata"))


#  Script is completed; clean up environment ----------------------------------
rm(list = ls())

#  Now, go to part 3
