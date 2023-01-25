#!/usr/bin/env Rscript

#  auxiliary_experiment-metadata.R
#  KA

#  Auxiliary function to record exprimenter, experiment data, and experiment
#+ description

callMetadataArguments <- function(parser = ap) {
    #TODO Documentation
    library(argparser)

    #  Metadata arguments
    ap <- add_argument(
        ap,
        short = "-u",
        arg = "--username",
        type = "character",
        default = "kga0",
        nargs = 1,
        help = "
            string for username directory under, e.g., '../data/', '../results/' <chr>
        "
    )
    ap <- add_argument(
        ap,
        short = "-e1",
        arg = "--experiment_date",
        type = "character",
        default = NA,
        nargs = 1,
        help = "string for date of experiment; use '<YYYY>_<MMDD>' format <chr>"
    )
    ap <- add_argument(
        ap,
        short = "-e2",
        arg = "--experiment_description",
        type = "character",
        default = NULL,
        nargs = 1,
        help = "string for description of experiment <chr>"
    )
    return(parser)
}
