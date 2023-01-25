#!/usr/bin/env Rscript

#  auxiliary.R
#  KA

#  Collection of auxiliary functions


#  Functions ------------------------------------------------------------------
importLibrary <- function(x) {
    # Suppress messages when loading libraries into the R session
    # 
    # :param x: a vector of libraries <chr>
    invisible(
        suppressWarnings(
            suppressMessages(
                lapply(x, library, character.only = TRUE, quietly = TRUE)
            )
        )
    )
}


#TODO `findReplace()` is no longer used: Check on this...
findReplace <- function(x, pattern, replacement = pattern, fill = NA, ...) {
    # Iterate over a vector to search for matched patterns and, whenever there is
    # a match, progressively replace the matches while reducing the amount of
    # comparisons
    # 
    # :param x: a vector
    # :param pattern: pattern(s) to search among elements
    # :param replacement: values with which to replace pattern-matched elements
    # :param fill: NA if not specified
    # :param ...: arguments to pass to grepl
    # :return: ...
    stopifnot(length(pattern) == length(replacement))
    answer <- rep_len(as.character(fill), length(x))    
    empty <- seq_along(x)
    for(i in seq_along(pattern)) {
        greps <- grepl(pattern[[i]], x[empty], ...)
        answer[empty[greps]] <- replacement[[i]]  
        empty <- empty[!greps]
    }
    return(answer)
}


printScriptBeginning <- function(script) {
    #TODO Documentation
    cat("\n# -----------------------------------------------\n")
    paste0("Beginning of script: ", script, ", ") %>% cat()
    Sys.time() %>% paste0() %>% cat(., "\n\n")
    cat("The following arguments have been called:\n")
    arguments[4:length(arguments)] %>% print(., quote = TRUE)
}


printScriptEnding <- function(script) {
    #TODO Documentation
    cat("End of script:", script, "\n\n")
}


generateExperimentName <- function() {
    #TODO Documentation
    library(argparser)

    experiment_name <- paste0(
        arguments$experiment_date, "_", arguments$experiment_description
    )
    return(experiment_name)
}


generateAnalysisName <- function(prefix = "experiment") {
    #TODO Documentation
    require(argparser)

    analysis_name <- paste0(
        prefix,
        "_", arguments$series,
        "_", arguments$alignment_strategy,
        "_", arguments$first_cell_type,
        "_", arguments$first_condition,
        "_", arguments$second_cell_type,
        "_", arguments$second_condition
    )
    return(analysis_name)
}


generateAnalysisNameTest <- function(prefix = "experiment") {
    #TODO Documentation
    require(argparser)

    analysis_name <- paste0(
        prefix,
        "_", arguments$in_string
    )
    return(analysis_name)
}


generateDataPath <- function(branch = "data") {
    # Create a file path that branches off from the parent of the current working
    # directory
    # 
    # :param string: String for the name of the directory to branch from parent
    # :return: String for file path
    library(argparser)

    path <- getwd() %>%
        paste0("/", branch, "/") %>% 
        file.path()
    return(path)
}


generateExperimentPath <- function(branch = "results", prefix = "experiment") {
    #TODO Documentation
    library(argparser)

    path <- getwd() %>%
    paste0(
        "/", branch, "/",
        arguments$username, "/",
        generateExperimentName(), "/",
        generateAnalysisName(prefix), "/"
    ) %>%
    file.path()
    return(path)
}


generateOutputDirectory <- function(branch = "results", prefix = "experiment") {
    #TODO Documentation
    library(argparser)  #TODO Not require() but some check, Is argparse loaded?

    path <- ifelse(
        !dir.exists(generateExperimentPath(branch, prefix)),
        dir.create(generateExperimentPath(branch, prefix), recursive = TRUE),
        FALSE
    )
    return(path)
}


makeOperation <- function(variable, command) {
    #TODO Documentation
    operation <- paste(variable, command) %>% paste(., collapse = "; ")
    return(operation)
}


loadCategoryRDS <- function(directory, string, file_suffix) {
    #TODO Documentation
    load <- readRDS(
        paste0(directory, "/", string, file_suffix, ".rds")
    )
    return(load)
}


reportCategoryRDS <- function(directory, string, file_suffix) {
    #TODO Documentation
    report <- paste0(
        "Have loaded ", directory, "/", string, file_suffix, ".rds"
    )
    return(report %>% cat(., "\n"))
}


# Set options and load packages -----------------------------------------------
options(
    # error = recover,
    scipen = 999
)
preload <- c(
    "argparser",
    "BiocParallel",
    "ggtext",
    "magrittr",
    "scales",
    "tidyverse"
)
importLibrary(preload)
