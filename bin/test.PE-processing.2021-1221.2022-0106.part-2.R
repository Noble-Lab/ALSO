#!/usr/bin/env Rscript

library(stringr)
library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)


###############################################################################

#  Run part 1 first; load .Rdata into environment
# -----------------------------------------------------------------------------
#  Temporary: Set up work directory (location TB∆)
directory_work <- "/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora"
setwd(directory_work)
rm(directory_work)

load("test.PE-processing.2021-1221.2022-0106.part-1.Rdata")
script <- "test.PE-processing.2021-1221.2022-0106.part-2.R"

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

# # #  Convert NAs
# #  For "ambiguous" and "mated," assign FALSE to NAs
# # operation <- paste0(
# #     variable[!str_detect(variable, "unmated")], "$discrepancy[",
# #         "is.na(", variable[!str_detect(variable, "unmated")], "$discrepancy", ")",
# #     "]", " <- ", "FALSE"
# # )  #TODO May delete this b/c no NAs
# # evaluateOperation(operation)
# 
# # #  For "unmated," assign TRUE to NAs
# # operation <- paste0(
# #     variable[str_detect(variable, "unmated")], "$discrepancy[",
# #         "is.na(", variable[str_detect(variable, "unmated")], "$discrepancy", ")",
# #     "]", " <- ", "TRUE"
# # )  #TODO May delete this b/c no NAs
# # evaluateOperation(operation)
# 
# #  Check to ensure that the NAs are gone now
# # command <- paste0(
# #     "<- ", "is.na(", variable, "$discrepancy", ")", " %>% ", "table()" 
# # )
# # operation <- makeOperation(paste0("isNA.", variable), command)
# # evaluateOperation(operation)  # No NAs
# 
# #  Remove the temporary tables
# # command <- paste0("rm(", "isNA.", variable, ")")
# # evaluateOperation(command)


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


#  From the sub tibbles, filter out duplicate $qpos, $qmpos -----------
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


###############################################################################

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


#  Correct filter.a.ambiguous.dedup.CAST --------------------------------------
filter.a.ambiguous.dedup.CAST[paste0("584", c(19, 20, 21, 22, 23, 24)), ] <-
    discrepancy.filter.a.ambiguous.dedup.CAST

#  Check
filter.a.ambiguous.dedup.CAST[c(58419:58429), ]
filter.a.ambiguous.dedup.CAST$row_n[c(58419:58429)]
# [1] 58419 58420 58421 58424 58422 58423 58425 58426 58427 58428 58429

#  Re-label $row_n
filter.a.ambiguous.dedup.CAST$row_n <- seq.int(
    nrow(filter.a.ambiguous.dedup.CAST)
)
filter.a.ambiguous.dedup.CAST[c(58419:58429), ]
filter.a.ambiguous.dedup.CAST$row_n[c(58419:58429)]
# [1] 58419 58420 58421 58422 58423 58424 58425 58426 58427 58428 58429


#  Correct filter.a.mated.dedup.129S1 -----------------------------------------
desired_order <- paste0("343", c(79, 80, 81, 83, 82, 84))
discrepancy.filter.a.mated.dedup.129S1.adjust <-
    discrepancy.filter.a.mated.dedup.129S1 %>%
        dplyr::mutate(row_n = factor(row_n, levels = desired_order)) %>%
        dplyr::arrange(row_n)
discrepancy.filter.a.mated.dedup.129S1.adjust$row_n

discrepancy.filter.a.mated.dedup.129S1.adjust$row_n <-
    discrepancy.filter.a.mated.dedup.129S1.adjust$row_n %>%
        as.character() %>%
        as.integer()

filter.a.mated.dedup.129S1[discrepancy.filter.a.mated.dedup.129S1$row_n, ] <-
    discrepancy.filter.a.mated.dedup.129S1.adjust

discrepancy.filter.a.mated.dedup.129S1 <-
    discrepancy.filter.a.mated.dedup.129S1.adjust

#  Check
filter.a.mated.dedup.129S1[c(34379:34389), ]
filter.a.mated.dedup.129S1$row_n[c(34379:34389)]
# [1] 34379 34380 34381 34383 34382 34384 34385 34386 34387 34388 34389

#  Re-label $row_n
filter.a.mated.dedup.129S1$row_n <- seq.int(
    nrow(filter.a.mated.dedup.129S1)
)

#  Check
filter.a.mated.dedup.129S1[c(34379:34389), ]
filter.a.mated.dedup.129S1$row_n[c(34379:34389)]
# [1] 34379 34380 34381 34382 34383 34384 34385 34386 34387 34388 34389

#  Remove unneeded variables
rm(desired_order, discrepancy.filter.a.mated.dedup.129S1.adjust)


###############################################################################

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


###############################################################################

#  Remove unneeded variables --------------------------------------------------
command <- paste0("rm(", "tbl.VA.", variable, ")")

variable_non_adjust <- variable %>% stringr::str_remove("adjust.")
variable_non_a <- variable_non_adjust %>% stringr::str_remove("a.")
command <- paste0(
    "rm(",
        "VMP.", variable, ", ",
        "VA.", variable, ", ",
        "tbl.VA.", variable, ", ",
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


###############################################################################

#  Generate $groupid_2 for munged "mated" and "ambiguous" tibbles -------------
generateVectorGroupID <- function(tibble) {
    output <- rep("", nrow(tibble))
    for (i in 1:nrow(tibble)) {
        if (i == nrow(tibble)) {
            paste0(
                deparse(substitute(tibble)), ": Completed"
            ) %>%
                print()
        } else if (
            tibble$pos[i] == tibble$mpos[i + 1] &
            tibble$mpos[i + 1] == tibble$pos[i]
        ) {
            output[i] <- i
            output[i + 1] <- i
        }
    }
    
    output <- output %>% as.integer()
    output <- output - ((rep(1:(length(output) / 2), 2) %>% sort()) - 1)
    return(output)
}

#  Work with only "mated" and "ambiguous" tibbles, not "unmated" tibbles
variable <- variable[!str_detect(variable, "unmated")]

#  Generate groupid.* variables for "mated" and "ambiguous" tibbles
variable_groupid <- paste0("groupid.", str_remove(variable, "adjust."))
command <- paste0("<- ", "generateVectorGroupID(", variable, ")")
operation <- makeOperation(variable_groupid, command)
evaluateOperation(operation)

#  Generate empty $groupid_2 columns for "mated" and "ambiguous" tibbles
command <- paste0("<- ", "rep(\"\", ", "nrow(", variable, "))")
operation <- makeOperation(paste0(variable, "$groupid_2"), command)
evaluateOperation(operation)

command <- paste0(
    "<- ", variable, " %>% ",
        "dplyr::relocate(groupid_2, .after = groupid)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  Assign group IDs to $groupid_2
for (i in 1:length(variable)) {
    if (
        length(eval(parse(text = variable_groupid[i]))) ==
        nrow(eval(parse(text = variable[i])))
    ) {
        print(paste0(variable[i], "$groupid_2 <- ", variable_groupid[i]))
        operation <- paste0(
            variable[i], "$groupid_2", " <- ", variable_groupid[i]
        )
        evaluateOperation(operation)
    } else {
        print(paste0(
                "length(", variable_groupid[i], ") != ",
                "nrow(", variable[i], ")", "\n",
                "Not assigning group IDs to $groupid_2"
        ))
    }
}

#  Remove unneeded variables
command <- paste0("rm(", "groupid.", str_remove(variable, "adjust."), ")")
evaluateOperation(command)
rm(variable_groupid)


###############################################################################

#  Merging rows ---------------------------------------------------------------
#  stackoverflow.com/questions/21003237/read-odd-numbered-rows-from-csv-file

#  Initialize one-line ("uniline") variables
variable_uniline <- paste0("uniline.", str_remove(variable, "adjust."))
mapply(
    assign, variable_uniline, "", MoreArgs = list(envir = parent.frame())
)

#  Concatenate even row to preceding odd row in  adjust.* variables
for (i in 1:length(variable)) {
    df <- eval(parse(text = variable[i]))
    odd <- seq(1, nrow(df), 2)
    even <- (seq(1, nrow(df), 2)) + 1
    
    odd <- df[odd, ]
    even <- df[even, ]
    
    df <- full_join(odd, even, by = "groupid_2")
    colnames(df) <- str_replace_all(colnames(df), "\\.x", "\\.odd")
    colnames(df) <- str_replace_all(colnames(df), "\\.y", "\\.even")
    
    command <- paste0("<- ", "df")
    operation <- makeOperation(variable_uniline[i], command)
    evaluateOperation(operation)
}

#  Remove unneeded variables
rm(df, odd, even)

#  Check flag statuses for adjust.*
for (i in 1:length(variable)) {
    print(variable[i])
    operation <- paste0(
        variable[i], "$flag %>% table()"
    )
    print(evaluateOperation(operation))
}

#  $flag.odd for uniline.*: Should be only flags 99 and 163
#QUESTION Is it a problem if not?
for (i in 1:length(variable_uniline)) {
    print(variable_uniline[i])
    operation <- paste0(
        variable_uniline[i], "$flag.odd %>% table()"
    )
    print(evaluateOperation(operation))
}

#  $flag.even for uniline.*: Should be only flags 83 and 147
#QUESTION Is it a problem if not?
for (i in 1:length(variable_uniline)) {
    print(variable_uniline[i])
    operation <- paste0(
        variable_uniline[i], "$flag.even %>% table()"
    )
    print(evaluateOperation(operation))
}
rm(i)


#  Final munging --------------------------------------------------------------
#NOTE
#  Whether merged or not (probably easier if I don't merge), I need to see if
#+ there's a SNP present in at least one of each group's $seq; if so, then as-
#+ sign that group appropriately
#+ 
#+ Ensure that the environment is clean and save an .Rdata image for a new
#+ subsequent script that will determine if SNPs are/aren't present per
#+ 
#+ Also, make sure there are separate .odd/.even $groupid_2 columns

command <- paste0(
    "<- ", variable_uniline, " %>% ",
        "dplyr::rename(groupid_2.odd = groupid_2)"
)
operation <- makeOperation(variable_uniline, command)
evaluateOperation(operation)

command <- paste0("<- ", variable_uniline, "$groupid_2.odd")
operation <- makeOperation(
    paste0(variable_uniline, "$groupid_2.even"), command
)
evaluateOperation(operation)

command <- paste0(
    "<- ", variable_uniline, " %>% ",
        "dplyr::relocate(groupid_2.even, .after = groupid.even)"
)
operation <- makeOperation(variable_uniline, command)
evaluateOperation(operation)


#  Save environment image -----------------------------------------------------
save.image(stringr::str_remove(script, ".R") %>% paste0(., ".Rdata"))


###############################################################################

#  Got to part 3
