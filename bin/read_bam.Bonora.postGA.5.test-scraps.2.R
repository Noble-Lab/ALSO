
library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)


# -----------------------------------------------------------------------------
#  Set up functions
`%notin%` <- Negate(`%in%`)


convertPercent <- function(x) {
    if(is.numeric(x)) {
        ifelse(is.na(x), x, paste0(round(x * 100L, 2), "%")) 
    } else x 
}


loadRDS <- readRDS


loadTibbleFromRDS <- function(variable, file) {
    command <- paste0(variable, " <- loadRDS(file = \"", file, "\")") %>% as.list()
    return(eval(parse(text = command), envir = globalenv()))
}


makeOperation <- function(variable, command) {
    operation <- paste(variable, command) %>% paste(., collapse = "; ")
    return(operation)
}


mungeTmpGB <- function(tbl) {
    #  Rename and select columns of interest
    tbl %>%
        plyr::rename(  # with plyr::rename(), old_name = new_name
            c(
                "rname" = "GB_rname",
                "strand" = "GB_strand",
                "pos" = "GB_pos",
                "mrnm" = "GB_mrnm",
                "mpos" = "GB_mpos"
            )
        ) %>% 
        dplyr::select(
            coordinate,
            GB_rname,
            GB_strand,
            GB_pos,
            GB_mrnm,
            GB_mpos
        )
}


#  Set up work directories (locations to be changed) --------------------------
path.1 <- "/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora"
path.2 <- "/segregatedReads.SNPTHRESH1.Q30"


#  Load tibbles (already sorted) from .rds files ------------------------------

#  KA assignments ---------------------
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
eval(parse(text = operation))

command <- paste0("rm(", variable, ")")
eval(parse(text = command))

#  GB assignments ---------------------
paste0(path.1, path.2) %>% setwd()
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
eval(parse(text = operation))

command <- paste0("rm(", variable, ")")
eval(parse(text = command))

rm(variable)


#  Identify intersections between "assignment" and "dedup.joint" --------------
for (i in c(1:4)) {
# for (i in c(4)) {
    command <- paste0(
        "<- ifelse(", variable_KA[i], "$coordinate %in% ", variable_GB[4], "$coordinate, '1', '0')"
    )
    operation <- makeOperation(paste0(variable_KA[i], "$in_ref_129S1"), command)
    print(operation)
    eval(parse(text = operation))
    
    command <- paste0(
        "<- ifelse(", variable_KA[i], "$coordinate %in% ", variable_GB[1], "$coordinate, '1', '0')"
    )
    operation <- makeOperation(paste0(variable_KA[i], "$in_alt_CAST"), command)
    print(operation)
    eval(parse(text = operation))
    
    command <- paste0(
        "<- ifelse(", variable_KA[i], "$coordinate %in% ", variable_GB[2], "$coordinate, '1', '0')"
    )
    operation <- makeOperation(paste0(variable_KA[i], "$in_ambig"), command)
    print(operation)
    eval(parse(text = operation))
    
    command <- paste0(
        "<- ifelse(", variable_KA[i], "$coordinate %in% ", variable_GB[3], "$coordinate, '1', '0')"
    )
    operation <- makeOperation(paste0(variable_KA[i], "$in_contra"), command)
    print(operation)
    eval(parse(text = operation))
    
    command <- paste0(
        "<- paste0(",
        variable_KA[i], "$in_ref_129S1, ",
        variable_KA[i], "$in_alt_CAST, ",
        variable_KA[i], "$in_ambig, ",
        variable_KA[i], "$in_contra)"
    )
    operation <- makeOperation(paste0(variable_KA[i], "$GB_intersection"), command)
    print(operation)
    eval(parse(text = operation))
    
    command <- paste0("<- ", variable_KA[i], "$GB_intersection %>% as.factor()")
    operation <- makeOperation(paste0(variable_KA[i], "$GB_intersection"), command)
    print(operation)
    eval(parse(text = operation))
}

#  Change the NAs of KA.NA$assignment with character-type "NA"
KA.NA$assignment[is.na(KA.NA$assignment)] <- "NA"
KA.NA$assignment <- KA.NA$assignment %>% as_factor()

#  Row-bind the KA tibbles
joint <- bind_rows(KA.129S1, KA.CAST, KA.Ambiguous, KA.NA)

#  Reorder the factor levels for column/variable GB_intersection
order_intersection <- c("1000", "0100", "0010", "0001", "0000")
joint$GB_intersection <- forcats::fct_relevel(joint$GB_intersection, order_intersection)

#  Copy the vector of binary names to a new variable; then, give them easier-
#+ to-read names 
joint$GB_assignment <- joint$GB_intersection %>% plyr::revalue(
    c(
        '0000' = 'GB.Not_present',
        '0001' = 'GB.Contra',
        '0010' = 'GB.Ambiguous',
        '0100' = 'GB.CAST',
        '1000' = 'GB.129S1'
    )
)

#  Give KA_assignment levels easy-to-read names
joint$KA_assignment <- joint$assignment %>% plyr::revalue(
    c(
        '129S1-SvImJ' = 'KA.129S1',
        'CAST-EiJ' = 'KA.CAST',
        'Neutral' = 'KA.Ambiguous',
        'NA' = 'KA.NA'
    )
)

#  Set up KA × GB assignment factors, KA × GB assignment factors
joint$KA_GB_assignment <- paste0(joint$KA_assignment, ' × ', joint$GB_assignment)
joint$GB_KA_assignment <- paste0(joint$GB_assignment, ' × ', joint$KA_assignment)

joint$KA_GB_assignment <- joint$KA_GB_assignment %>% as_factor()
joint$GB_KA_assignment <- joint$GB_KA_assignment %>% as_factor()

#TODO Reorder these factor levels

setwd(path.1)
saveRDS(joint, file = "KA_GB.joint.rds")
