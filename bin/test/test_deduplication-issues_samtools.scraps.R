
library(tidyverse)
library(Rsamtools)


evaluateOperation <- function(operation = operation) {
    return(eval(parse(text = operation), envir = .GlobalEnv))
}


makeOperation <- function(variable = variable, command = command) {
    operation <- paste0(variable, " <- ", command)
    return(operation)
}


#  Load in .bam information, including mate information -----------------------
path <- "/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora/test_deduplication-issues"
setwd(path)

# chromosome <- "chr1"
chromosome <- "chrX"

#  Load .bam files
file <- list.files(pattern = paste0("\\.rmdup.", chromosome, ".bam$")) %>%
    stringr::str_subset("mm10\\.", negate = TRUE)
variable <- file %>%
    strsplit(., "\\.") %>%
    lapply(., `[[`, 1) %>% unlist() %>%
    paste0("dedup.", .) %>%
    stringr::str_remove("S1")
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

#  Load .bai files
file <- list.files(pattern = paste0("\\.rmdup.", chromosome, ".bam.bai$")) %>%
    stringr::str_subset("mm10\\.", negate = TRUE)
index <- file %>%
    strsplit(., "\\.") %>%
    lapply(., `[[`, 1) %>% unlist() %>%
    paste0("index.", .) %>%
    stringr::str_remove("S1")
mapply(
    assign, index, file, MoreArgs = list(envir = parent.frame())
)

#  Using Rsamtools, load in standard .bam fields
command <- paste0(
    variable, " %>% ",
        "Rsamtools::BamFile(., index = ", index, ", asMates = TRUE)", " %>% ",
        "Rsamtools::scanBam()", " %>% ",
        "as.data.frame()", " %>% ",
        "tibble::as_tibble()"
)
operation <- makeOperation(paste0(variable, ".full"), command)
evaluateOperation(operation)

#  Using Rsamtools, load in .bam AS and MD fields
map_params <- Rsamtools::ScanBamParam(tag = c("AS", "MD"))
command <- paste0(
    variable, " %>% ",
        "Rsamtools::BamFile(., index = ", index, ", asMates = TRUE)", " %>% ",
        "Rsamtools::scanBam(param = map_params)", " %>% ",
        "as.data.frame()", " %>% ",
        "tibble::as_tibble()"
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

#  Remove variables no longer needed
operation <- paste0("rm(", variable, ".full, ", index, ")")
evaluateOperation(operation)

rm(file, index)