
#  Is $qpos in $qmpos? Vice versa? Create tables and tibbles of answers -------
command <- paste0(
    "<- testPosInMpos(",
        "pos = ", variable, "$qpos, ",
        "mpos = ", variable, "$qmpos",
    ") %>% ",
        "table()"
)
operation <- makeOperation(paste0("PIM.", variable), command)
evaluateOperation(operation)

command <- paste0(
    "<- testMposInPos(",
        "pos = ", variable, "$qpos, ",
        "mpos = ", variable, "$qmpos",
    ") %>% ",
        "table()"
)
operation <- makeOperation(paste0("MIP.", variable), command)
evaluateOperation(operation)

#  Convert the tables to tibbles: PIM
command <- paste0(
    "<- PIM.", variable, " %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("PIM.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles:" PIM
command <- paste0(
    "<- ", "tibble::tibble(", "PIM.", variable, ")"
)
operation <- makeOperation(paste0("PIM.", variable), command)
evaluateOperation(operation)

vector_string_tibbles <- paste0("PIM.", variable)
z.check.1.PIM <- createTibbleFromTibbles(vector_string_tibbles)
z.check.1.PIM$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.1.PIM <- z.check.1.PIM %>% dplyr::arrange(rownames)

#  Convert the tables to tibbles: MIP
command <- paste0(
    "<- MIP.", variable, " %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("MIP.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles:" MIP
command <- paste0(
    "<- ", "tibble::tibble(", "MIP.", variable, ")"
)
operation <- makeOperation(paste0("MIP.", variable), command)
evaluateOperation(operation)

vector_string_tibbles <- paste0("MIP.", variable)
z.check.1.MIP <- createTibbleFromTibbles(vector_string_tibbles)
z.check.1.MIP$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.1.MIP <- z.check.1.MIP %>% dplyr::arrange(rownames)


#  Is PIM.*$`FALSE` == MIP.*$`FALSE`? -----------------------------------------
command <- paste0(
    "<- PIM.", variable, "$`FALSE` == MIP.", variable, "$`FALSE`"
)
operation <- makeOperation(paste0("PeqM.", variable), command)
evaluateOperation(operation)

#  Create a "tibble of tibbles"
command <- paste0(
    "<- ", "tibble::tibble(", "PeqM.", variable, ")"
)
operation <- makeOperation(paste0("PeqM.", variable), command)
evaluateOperation(operation)

operation <- paste0("colnames(", "PeqM.", variable, ") ", "<- ", "\"n\"")
evaluateOperation(operation)

vector_string_tibbles <- paste0("PeqM.", variable)
z.check.2.PeqM <- createTibbleFromTibbles(vector_string_tibbles)
z.check.2.PeqM$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.2.PeqM <- z.check.2.PeqM %>% dplyr::arrange(rownames)


#  Is PIM.*$`FALSE` > MIP.*$`FALSE`? ------------------------------------------
command <- paste0(
    "<- PIM.", variable, "$`FALSE` > MIP.", variable, "$`FALSE`"
)
operation <- makeOperation(paste0("PgtM.", variable), command)
evaluateOperation(operation)

#  Create a tibble
command <- paste0(
    "<- ", "tibble::tibble(", "PgtM.", variable, ")"
)
operation <- makeOperation(paste0("PgtM.", variable), command)
evaluateOperation(operation)

operation <- paste0("colnames(", "PgtM.", variable, ") ", "<- ", "\"n\"")
evaluateOperation(operation)

vector_string_tibbles <- paste0("PgtM.", variable)
z.check.2.PgtM <- createTibbleFromTibbles(vector_string_tibbles)
z.check.2.PgtM$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.2.PgtM <- z.check.2.PgtM %>% dplyr::arrange(rownames)


#  Is PIM.*$`FALSE` < MIP.*$`FALSE`? ------------------------------------------
command <- paste0(
    "<- PIM.", variable, "$`FALSE` < MIP.", variable, "$`FALSE`"
)
operation <- makeOperation(paste0("PltM.", variable), command)
evaluateOperation(operation)

#  Create a tibble
command <- paste0(
    "<- ", "tibble::tibble(", "PltM.", variable, ")"
)
operation <- makeOperation(paste0("PltM.", variable), command)
evaluateOperation(operation)

operation <- paste0("colnames(", "PltM.", variable, ") ", "<- ", "\"n\"")
evaluateOperation(operation)

vector_string_tibbles <- paste0("PltM.", variable)
z.check.2.PltM <- createTibbleFromTibbles(vector_string_tibbles)
z.check.2.PltM$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.2.PltM <- z.check.2.PltM %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0(
    "rm(",
        "PIM.", variable, ", ",
        "MIP.", variable, ", ",
        "PeqM.", variable, ", ",
        "PltM.", variable, ", ",
        "PgtM.", variable,
    ")"
)
evaluateOperation(command)
rm(vector_string_tibbles)


#  Make logical vectors for $qpos in $qmpos and vice versa --------------------

#  PIM
command <- paste0(
    "<- testPosInMpos(pos = ", variable, "$qpos, mpos = ", variable, "$qmpos)"
)
operation <- makeOperation(paste0("vPIM.", variable), command)
evaluateOperation(operation)

# #  Check: NAs in vPIM vector?
# command <- paste0("<- ", "is.na(", "vPIM.", variable, ")")
# operation <- makeOperation(paste0("isNA.vPIM.", variable), command)
# evaluateOperation(operation)
# 
# command <- paste0(
#     paste0("<- isNA.vPIM.", variable, " %>% ", "table()")
# )
# operation <- makeOperation(paste0("isNA.vPIM.", variable), command)
# evaluateOperation(operation)  # No NAs
# 
# command <- paste0("rm(", "isNA.vPIM.", variable, ")")
# evaluateOperation(command)

#  MIP
command <- paste0(
    "<- testMposInPos(pos = ", variable, "$qpos, mpos = ", variable, "$qmpos)"
)
operation <- makeOperation(paste0("vMIP.", variable), command)
evaluateOperation(operation)

#  Check: NAs in vMIP vector?
# command <- paste0("<- ", "is.na(", "vMIP.", variable, ")")
# operation <- makeOperation(paste0("isNA.vMIP.", variable), command)
# evaluateOperation(operation)
# 
# command <- paste0(
#     paste0("<- isNA.vMIP.", variable, " %>% ", "table()")
# )
# operation <- makeOperation(paste0("isNA.vMIP.", variable), command)
# evaluateOperation(operation)  # No NAs
# 
# command <- paste0("rm(", "isNA.vMIP.", variable, ")")
# evaluateOperation(command)

#NOTE
#  Based on PeqM.*, PgtM.*, and PltM.* results, use the following vectors:
#+ All vPIM.* are fine, apparently

#  Remove unneeded variables
command <- paste0(
    "rm(", "vMIP.", variable, ")"
)
evaluateOperation(command)


#  Subset tibbles based on logical vectors for $qpos in $qmpos ----------------

#NOTE vPIM: logical vectors for $qpos %in% $qmpos
command <- paste0("<- ", variable, "[", "vPIM.", variable, ", ]")
operation <- makeOperation(paste0("sub.", variable), command)
evaluateOperation(operation)

#  Remove unneeded vPIM vectors
command <- paste0("rm(", "vPIM.", variable, ")")
evaluateOperation(command)


#  For the subsetted tibbles, save qpos %in$ qmpos, vice versa in vectors -----
variable <- paste0("sub.", variable)
command <- paste0(  # Table qpos %in$ qmpos
    "<- testPosInMpos(",
        "pos = ", variable, "$qpos, ",
        "mpos = ", variable, "$qmpos",
    ") %>% ",
        "table()"
)
operation <- makeOperation(paste0("PIM.", variable), command)
evaluateOperation(operation)

command <- paste0(  # Table qmpos %in$ qpos
    "<- testMposInPos(",
        "pos = ", variable, "$qpos, ",
        "mpos = ", variable, "$qmpos",
    ") %>% ",
        "table()"
)
operation <- makeOperation(paste0("MIP.", variable), command)
evaluateOperation(operation)

#  Convert the tables to tibbles
command <- paste0(
    "<- PIM.", variable, " %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("PIM.", variable), command)
evaluateOperation(operation)

command <- paste0(
    "<- MIP.", variable, " %>% ",
        "t() %>% ",
        "cbind() %>% ",
        "as_tibble()"
)
operation <- makeOperation(paste0("MIP.", variable), command)
evaluateOperation(operation)

#  Collect pertinent readouts for the subsetted tibbles in one tibble; i.e.,
#+ create a "tibble of tibbles"

#  PIM
vector_string_tibbles <- paste0("PIM.", variable)
z.check.sub.0.PIM <- createTibbleFromTibbles(vector_string_tibbles)
z.check.sub.0.PIM$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.sub.0.PIM <- z.check.sub.0.PIM %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0("rm(", "PIM.", variable, ")")
evaluateOperation(command)
rm(vector_string_tibbles)

#  MIP
vector_string_tibbles <- paste0("MIP.", variable)
z.check.sub.0.MIP <- createTibbleFromTibbles(vector_string_tibbles)
z.check.sub.0.MIP$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
z.check.sub.0.MIP <- z.check.sub.0.MIP %>% dplyr::arrange(rownames)

#  Remove unneeded variables
command <- paste0("rm(", "MIP.", variable, ")")
evaluateOperation(command)
rm(vector_string_tibbles)


#  Order the tibbles based on $isize_abs, $qpos, and $qmpos -------------------
command <- paste0(
    "<- ", variable, "[",
        "order(",
            variable, "$isize_abs, ",
            variable, "$qpos, ",
            variable, "$qmpos",
        "), ",
    "]"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Save environment image -----------------------------------------------------
save.image(stringr::str_remove(script, ".R") %>% paste0(., ".Rdata"))


#  Script is completed; clean up environment ----------------------------------
rm(list = ls())

#  Now, go to part 2


###############################################################################
#  Check .bam $flag and $mate_status ------------------------------------------

#  Check $flag
variable_flag <- paste0("flag.", variable)
command <- paste0("<- ", variable, "$flag %>% table()")
operation <- makeOperation(variable_flag, command)
evaluateOperation(operation)

flag.dedup.129S1
flag.liftedToMm10.CAST

#  Remove variables no longer needed
command <- paste0("rm(", variable_flag, ")")
evaluateOperation(command)

rm(variable_flag)

#  Check $mate_status
variable_status <- paste0("status.", variable)
command <- paste0("<- ", variable, "$mate_status %>% table()")
operation <- makeOperation(variable_status, command)
evaluateOperation(operation)

#  Print out $mate_status
for (i in 1:length(variable_status)) {
    cat(variable_status[i], "\n")
    cat("mated ambiguous unmated\n")
    command <- paste0(variable_status[i])
    evaluateOperation(command) %>% cat("\n\n")
}

#  Remove variables no longer needed
command <- paste0("rm(", variable_status, ")")
evaluateOperation(command)

rm(i, variable_status)

# #  Remove flags not equal to...
# #+ 65, 129, 67, 131, 113, 177, 81, 161, 83, 163, 97, 145, 99, 147
# #+ e.g., see seqanswers.com/forums/showthread.php?t=17314
# command <- paste0(
#     "<- ", variable, "[",
#         variable, "$flag %in% c(65, 129, 67, 131, 113, 177, 81, 161, 83, 163, 97, 145, 99, 147)", ", ",
#     "]"
# )
# operation <- makeOperation(variable, command)
# evaluateOperation(operation)
# 
# #  Check $flag again
# variable_flag <- paste0("flag.", variable)
# command <- paste0("<- ", variable, "$flag %>% table()")
# operation <- makeOperation(variable_flag, command)
# evaluateOperation(operation)
# 
# flag.liftedToMm10.129S1
# flag.liftedToMm10.CAST
# 
# #  Remove variables no longer needed
# command <- paste0("rm(", variable_flag, ")")
# evaluateOperation(command)
# 
# rm(variable_flag)


# #  Split tibbles based on $mate_status ----------------------------------------
# split <- c("ambiguous", "mated", "unmated")
# for (i in 1:length(split)) {
#     variable_split <- paste0(split[i], ".", variable)
#     command <- paste0(
#         "<- ", variable, "[", variable, "$mate_status == '", split[i], "', ]"
#     )
#     operation <- makeOperation(variable_split, command)
#     evaluateOperation(operation)
# }
# rm(i)
# 
# # #  Save initial variables as 'z.init.*'
# # command <- paste0("<- ", variable)
# # operation <- makeOperation(paste0("z.init.", variable), command)
# # evaluateOperation(operation)
# 
# #  Remove initial variables
# command <- paste0("rm(", variable, ")")
# evaluateOperation(command)
# 
# #  Create a vector of names of "split variables"
# split <- paste0(split, ".")
# variable <- purrr::cross(list(split, variable)) %>%
#     purrr::map_chr(paste0, collapse = "")
# 
# #  Clean up variables
# rm(split, variable_split)

#  For each variable, print out duplicated status based on $pos
for (i in 1:length(variable)) {
    paste0(variable[i], "$pos") %>% cat("\n")
    cat("FALSE TRUE\n")
    command <- paste0(variable[i], "$pos %>% duplicated() %>% table()")
    evaluateOperation(command) %>% cat()
    cat("\n\n")
}
rm(i)

# #  Create backup of variable 'variable'
# variable <- z.bak.variable
# z.bak.variable <- variable


#  Reorder rname factor levels ------------------------------------------------
command <- paste0(
    "<- forcats::fct_relevel(", variable, "$rname, chromosome)"
)
operation <- makeOperation(paste0(variable, "$rname"), command)
evaluateOperation(operation)


#  Drop unused rname factor levels --------------------------------------------
command <- paste0(
    "<- ", variable, "$rname %>% forcats::fct_drop()"
)
operation <- makeOperation(paste0(variable, "$rname"), command)
evaluateOperation(operation)


#  Drop rows that are not chr1-19, chrX ---------------------------------------
chromosomes <- c(paste0("chr", 1:19), "chrX")
command <- paste0(
    "<- ", variable, " %>% ", "filter(., rname %in% chromosomes)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  Remove unneeded variables
rm(chromosome, chromosomes)


#  Create and append $pos_end, $mpos_end --------------------------------------
command <- paste0("<- ", variable, "$pos + 49")
operation <- makeOperation(paste0(variable, "$pos_end"), command)
evaluateOperation(operation)

command <- paste0("<- ", variable, "$mpos + 49")
operation <- makeOperation(paste0(variable, "$mpos_end"), command)
evaluateOperation(operation)

#  Move *_end to appropriate locations
command <- paste0(
    "<- ", variable, " %>% ",
        "dplyr::relocate(pos_end, .after = pos) %>% ",
        "dplyr::relocate(mpos_end, .after = mpos)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Filter out rows with mapq less than 30 -------------------------------------
command <- paste0("<- ", variable, "[", variable, "$mapq >= 30, ]")
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Prior to sorting, deduplicate and mate-label the tibbles -------------------

#  Set up $criteria, a variable needed for sorting: qname, flag, pos, mpos
command <- paste0(
    "<- paste0(",
        variable, "$qname, ", "'_', ",
        variable, "$flag, ", "'_', ",
        variable, "$pos, ", "'_', ",
        variable, "$mpos",
    ")"
)
operation <- makeOperation(paste0(variable, "$criteria"), command)
evaluateOperation(operation)

#  Set up $qpos, a variable needed for sorting
command <- paste0(
    "<- paste0(",
        variable, "$qname, ", "'_', ",
        variable, "$pos",
    ")"
)
operation <- makeOperation(paste0(variable, "$qpos"), command)
evaluateOperation(operation)

#  Set up $qmpos, a variable needed for sorting
command <- paste0(
    "<- paste0(",
        variable, "$qname, ", "'_', ",
        variable, "$mpos",
    ")"
)
operation <- makeOperation(paste0(variable, "$qmpos"), command)
evaluateOperation(operation)

#  Set up $isize_abs, i.e., absolute insert-size values, for sorting
command <- paste0(
    "<- ", variable, "$isize %>% abs() %>% as.numeric()"
)
operation <- makeOperation(paste0(variable, "$isize_abs"), command)
evaluateOperation(operation)

# #  Create backups of variables
# command <- paste0("<- ", variable)
# operation <- makeOperation(paste0("z.bak.", variable), command)
# evaluateOperation(operation)
# 
# #  Load back-ups of variables
# command <- paste0("<- z.bak.", variable)
# operation <- makeOperation(variable, command)
# evaluateOperation(operation)

# #  To keep the original munged tibbles together, rename the tibbles associated
# #+ with 'variable'
# command <- paste0("<- ", variable)
# operation <- makeOperation(paste0("a.", variable), command)
# evaluateOperation(operation)
# 
# #  Remove unneeded variables
# command <- paste0("rm(", variable, ")")
# evaluateOperation(command)
# 
# #  Adjust values associated 'variable'
# variable <- paste0("a.", variable)

#  Rearrange the order of tibble columns for quick, easy surveying
command <- paste0(
    "<- ", variable, " %>% ",
        "dplyr::relocate(mpos, .after = pos) %>% ",
        "dplyr::relocate(mpos_end, .after = pos_end) %>% ",
        "dplyr::relocate(mrnm, .after = rname) %>% ",
        "dplyr::relocate(isize, .after = mpos_end) %>% ",
        "dplyr::relocate(criteria, .before = groupid) %>% ",
        "dplyr::relocate(AS, .after = cigar) %>% ",
        "dplyr::relocate(qpos, .after = criteria) %>% ",
        "dplyr::relocate(qmpos, .after = qpos) %>% ",
        "dplyr::relocate(qname, .after = qual) %>% ",
        "dplyr::relocate(isize_abs, .before = isize)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  To survey duplicates, split out entries in which $qpos is the same as
#+ $qmpos, i.e., entries in which each member of the pair maps to the same
#+ location
command <- paste0(
    "<- ", variable, "[", variable, "$qpos == ", variable, "$qmpos, ]"
)
operation <- makeOperation(paste0("same.", variable), command)
evaluateOperation(operation)

#  Regarding the main tibbles, remove entries in which $qpos is the same as
#+ $qmpos, i.e., entries in which each member of the pair maps to the same
#+ location
command <- paste0(
    "<- ", variable, "[!(", variable, "$qpos == ", variable, "$qmpos), ]"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#QUESTION?  Does the below leave one of the entries post deduplication? Yes.
#+
#+ Deduplicate the main tibbles based on $criteria: There should be no more
#+ than one entry with a given $criteria value (combination of $qname, $flag,
#+ $pos, and $mpos)
command <- paste0(
    "<- ", variable, "[!duplicated(", variable, "$criteria), ]"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  Remove the unneeded "same" variables
operation <- paste0("rm(", "same.", variable, ")")
evaluateOperation(operation)


#  Save environment image -----------------------------------------------------
save.image(stringr::str_remove(script, ".R") %>% paste0(., ".Rdata"))

#  Go to test.PE-processing.part-7.R
