#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(reshape2)
library(scales)
library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)


#  Set up work directory (locations TB∆) --------------------------------------
directory_user <- "/Users/kalavattam"
directory_base <- "Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do"
directory_work <- "get_unique_fragments/Bonora"

path.1 <- paste0(directory_user, "/", directory_base, "/", directory_work)
path.2 <- "segregatedReads.SNPTHRESH1.Q30"
rm(directory_user, directory_base, directory_work)


#  Set up functions -----------------------------------------------------------
`%notin%` <- Negate(`%in%`)


convertPercent <- function(x) {
    if(is.numeric(x)) {
        ifelse(is.na(x), x, paste0(round(x * 100L, 2), "%")) 
    } else {
        x
    } 
}


evaluateOperation <- function(operation = operation) {
    return(eval(parse(text = operation), envir = .GlobalEnv))
}


loadRDS <- readRDS


loadTibbleFromRDS <- function(variable, file) {
    command <- paste0(variable, " <- loadRDS(file = \"", file, "\")") %>%
        as.list()
    return(eval(parse(text = command), envir = globalenv()))
}


makeOperation <- function(variable, command) {
    operation <- paste(variable, command) #%>% paste(., collapse = "; ")
    return(operation)
}


mungeTmpGB <- function(tbl) {
    #  Rename and select columns of interest
    tbl %>%
        plyr::rename(c(
                "rname" = "GB_rname",
                "strand" = "GB_strand",
                "pos" = "GB_pos",
                "mrnm" = "GB_mrnm",
                "mpos" = "GB_mpos"
        )) %>% 
        dplyr::select(
            coordinate,
            GB_rname,
            GB_strand,
            GB_pos,
            GB_mrnm,
            GB_mpos
        )
}


#  Load tibbles (already sorted) from .rds files ------------------------------

#  Load KA assignments
path.1 %>% setwd()

# chromosome <- "chr1"
chromosome <- "chrX"

search <- paste0("\\.", chromosome, ".rds$")
file <- c(list.files(path = ".", pattern = search))
variable <- file %>% gsub("\\.rds$", "", .)
loadTibbleFromRDS(variable = variable, file = file)
#NOTE Remember, "alt" is "CAST", "ref" is "129"

variable_KA <- c(
    "KA.129S1",
    "KA.CAST",
    "KA.NA",
    "KA.Ambiguous"
)
command <- paste0("<- ", variable)
operation <- makeOperation(variable_KA, command)
evaluateOperation(operation)

command <- paste0("rm(", variable, ")")
eval(parse(text = command))

#  Load GB assignments
paste0(path.1, "/", path.2) %>% setwd()

search <- paste0("^(dedup.joint.*", chromosome, ".rds$)")
file <- c(list.files(path = ".", pattern = search))
variable <- file %>% gsub("\\.rds$", "", .)
loadTibbleFromRDS(variable = variable, file = file)

variable_GB <- c(
    "GB.alt.CAST",
    "GB.ambig",
    "GB.contra",
    "GB.ref.129S1"
)
command <- paste0("<- ", variable)
operation <- makeOperation(variable_GB, command)
evaluateOperation(operation)

command <- paste0("rm(", variable, ")")
eval(parse(text = command))

rm(variable)

#  Back-ups for testing
# bak.KA.129S1 <- KA.129S1
# bak.KA.CAST <- KA.CAST
# bak.KA.Ambiguous <- KA.Ambiguous
# bak.KA.NA <- KA.NA
# 
# KA.129S1 <- bak.KA.129S1
# KA.CAST <- bak.KA.CAST
# KA.Ambiguous <- bak.KA.Ambiguous
# KA.NA <- bak.KA.NA


#  "Smart join" of GB to KA tibbles, then munging -----------------------------

#  Join the GB data; then, munge the column names
joint.GB <- dplyr::bind_rows(GB.alt.CAST, GB.ref.129S1, GB.ambig, GB.contra)
colnames(joint.GB) <- paste0(colnames(joint.GB), ".GB")
colnames(joint.GB)[3] <- "coordinate"
colnames(joint.GB)[colnames(joint.GB) %>% length()] <- "assignment_GB"

#  Change NA in KA.NA$assignment to character "NA"
KA.NA$assignment[is.na(KA.NA$assignment)] <- "NA"

#  Come up with names for tibbles after they're joined; then, do the joins
variable_joint <- variable_KA %>% strsplit(., "\\.") %>% lapply(., `[[`, 2) %>% unlist() %>% paste0("KA_GB.", .)
command <- paste0(
    "<- full_join(",
        variable_KA, ", ",
        "joint.GB %>% select(-ID.GB, -read.GB, -tag.XS.GB, -tag.NM.GB), ",
        "by = \"coordinate\"",
    ") %>% ",
        "tidyr::drop_na(\"assignment\")"
)
operation <- makeOperation(variable_joint, command)
evaluateOperation(operation)

#  Change NA in KA_GB.*$assignment_GB to character "GB.not_present"
command <- paste0(
    variable_joint, "$assignment_GB[",
        "is.na(", variable_joint, "$assignment_GB)",
    "] <- \"GB.not_present\""
)
eval(parse(text = command))

#  Rename $assignment_GB to $GB_assignment.1
command <- paste0(
    "<- ", variable_joint, " %>% ",
        "dplyr::rename(., GB_assignment.1 = assignment_GB)"
)  # new_name = old_name
operation <- makeOperation(variable_joint, command)
evaluateOperation(operation)

#  Clean up: Remove the prejoined KA.* tibbles
command <- paste0("rm(", variable_KA, ")")
eval(parse(text = command))


#  Identify intersections between "assignment" and "dedup.joint" --------------
#TODO 1/2 Because of the "smart joins," this section is no longer needed;
#TODO 2/2 however, keep it for now
for (i in c(1:4)) {
    #  "GB.ref.129S1"
    command <- paste0(
        "<- ifelse(",
            variable_joint[i], "$coordinate %in% ", variable_GB[4], "$coordinate,
            '1',
            '0'",
        ")"
    )
    operation <- makeOperation(
        paste0(variable_joint[i], "$in_ref_129S1"), command
    )
    # print(operation)
    evaluateOperation(operation)
    
    #  "GB.alt.CAST"
    command <- paste0(
        "<- ifelse(",
            variable_joint[i], "$coordinate %in% ", variable_GB[1], "$coordinate,
            '1',
            '0'",
        ")"
    )
    operation <- makeOperation(
        paste0(variable_joint[i], "$in_alt_CAST"), command
    )
    # print(operation)
    evaluateOperation(operation)
    
    #  GB.ambig
    command <- paste0(
        "<- ifelse(",
            variable_joint[i], "$coordinate %in% ", variable_GB[2], "$coordinate,
            '1',
            '0'",
        ")"
    )
    operation <- makeOperation(paste0(variable_joint[i], "$in_ambig"), command)
    # print(operation)
    evaluateOperation(operation)
    
    #  GB.contra
    command <- paste0(
        "<- ifelse(",
            variable_joint[i], "$coordinate %in% ", variable_GB[3], "$coordinate,
            '1',
            '0'",
        ")"
    )
    operation <- makeOperation(
        paste0(variable_joint[i], "$in_contra"), command
    )
    # print(operation)
    evaluateOperation(operation)
    
    command <- paste0(
        "<- paste0(",
            variable_joint[i], "$in_ref_129S1, ",
            variable_joint[i], "$in_alt_CAST, ",
            variable_joint[i], "$in_ambig, ",
            variable_joint[i], "$in_contra",
        ")"
    )
    operation <- makeOperation(
        paste0(variable_joint[i], "$GB_intersection"), command
    )
    # print(operation)
    evaluateOperation(operation)
    
    command <- paste0(
        "<- ", variable_joint[i], "$GB_intersection %>% as.factor()"
    )
    operation <- makeOperation(
        paste0(variable_joint[i], "$GB_intersection"), command
    )
    # print(operation)
    evaluateOperation(operation)
}

#  Reorder the factor levels for column/variable GB_intersection
order_intersection <- c("1000", "0100", "0010", "0001", "0000")
command <- paste0(
    "<- forcats::fct_relevel(",
        variable_joint, "$GB_intersection, order_intersection",
    ")"
)
operation <- makeOperation(
    paste0(variable_joint, "$GB_intersection"), command
)
evaluateOperation(operation)

#  Copy the vector of binary names to a new variable; then, give them easier-
#+ to-read names 
command <- paste0("<- ", variable_joint, "$GB_intersection")
operation <- makeOperation(
    paste0(variable_joint, "$GB_assignment.2"), command
)
evaluateOperation(operation)

command <- paste0(
    "<- ", variable_joint, "$GB_assignment.2 %>% plyr::revalue(c(",
        "'0000' = 'GB.Not_present', ",
        "'0001' = 'GB.Contra', ",
        "'0010' = 'GB.Ambiguous', ",
        "'0100' = 'GB.CAST', ",
        "'1000' = 'GB.129S1'",
    "))"
)
operation <- makeOperation(
    paste0(variable_joint, "$GB_assignment.2"), command
)
evaluateOperation(operation)

command <- paste0(
    "<- ", variable_joint, "$assignment %>% plyr::revalue(c(",
        "'129S1-SvImJ' = 'KA.129S1', ",
        "'CAST-EiJ' = 'KA.CAST', ",
        "'Neutral' = 'KA.Ambiguous', ",
        "'NA' = 'KA.NA' ",
    "))"
)
operation <- makeOperation(paste0(variable_joint, "$KA_assignment"), command)
evaluateOperation(operation)  # The errors reported here are fine...


#  Set up KA × GB assignment factors, KA × GB assignment factors --------------
KA_GB.129S1$KA_GB_assignment <- paste0(
    KA_GB.129S1$KA_assignment, ' × ', KA_GB.129S1$GB_assignment.2
)
KA_GB.CAST$KA_GB_assignment <- paste0(
    KA_GB.CAST$KA_assignment, ' × ', KA_GB.CAST$GB_assignment.2
)
KA_GB.NA$KA_GB_assignment <- paste0(
    KA_GB.NA$KA_assignment, ' × ', KA_GB.NA$GB_assignment.2
)
KA_GB.Ambiguous$KA_GB_assignment <- paste0(
    KA_GB.Ambiguous$KA_assignment, ' × ', KA_GB.Ambiguous$GB_assignment.2
)

KA_GB.129S1$GB_KA_assignment <- paste0(
    KA_GB.129S1$GB_assignment.2, ' × ', KA_GB.129S1$KA_assignment
)
KA_GB.CAST$GB_KA_assignment <- paste0(
    KA_GB.CAST$GB_assignment.2, ' × ', KA_GB.CAST$KA_assignment
)
KA_GB.NA$GB_KA_assignment <- paste0(
    KA_GB.NA$GB_assignment.2, ' × ', KA_GB.NA$KA_assignment
)
KA_GB.Ambiguous$GB_KA_assignment <- paste0(
    KA_GB.Ambiguous$GB_assignment.2, ' × ', KA_GB.Ambiguous$KA_assignment
)

command <- paste0("<- ", variable_joint, "$KA_GB_assignment %>% as_factor()")
operation <- makeOperation(
    paste0(variable_joint, "$KA_GB_assignment"), command
)
evaluateOperation(operation)

command <- paste0("<- ", variable_joint, "$GB_KA_assignment %>% as_factor()")
operation <- makeOperation(
    paste0(variable_joint, "$GB_KA_assignment"), command
)
evaluateOperation(operation)

#TODO Reorder these factor levels


#  Row-bind the KA tibbles ----------------------------------------------------
joint <- bind_rows(KA_GB.129S1, KA_GB.CAST, KA_GB.Ambiguous, KA_GB.NA)

setwd(path.1)
saveRDS(joint, file = paste0("KA_GB.joint.", chromosome,".rds"))


#  Script is over; clean up the environment -----------------------------------
rm(list = ls())
