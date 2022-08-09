#!/bin/Rscript

#  munge_STAR_log_files.R
#  KA


#  Load libraries, set up options ---------------------------------------------
library(ggplot2)
library(tidyverse)
options(pillar.sigfig = 8, scipen = 10000)


#  Functions ------------------------------------------------------------------
extract_samples <- function(files) {
    #  rdrr.io/github/HuntsmanCancerInstitute/hciR/src/R/extract_samples.R
    
    #' Parse sample from file names
    #'
    #' Sample names are parsed from file names without the extension.  If the file name is not unique,
    #'    the parent directory is used.
    #'
    #' @param files file name with path
    #'
    #' @return A vector
    #'
    #' @author Chris Stubben
    #'
    #' @examples
    #'     # File name or parent directory should be unique
    #'     extract_samples(c("align1/1355X1.counts", "align2/1355X2.counts"))
    #'     extract_samples(c("align1/1355X1/Log.out", "align2/1355X2/Log.out"))
    #' @export
    
    # Capture sample in filename or path
    x <- lapply(strsplit(files, "/"), rev)
    samples <- sapply(x, "[", 1)
    
    # Remove file extension
    samples <- tools::file_path_sans_ext(samples)
    if(any(duplicated(samples))) {
        # Use parent directory, 13555X2/Log.final.out
        samples <- sapply(x, "[", 2)
    }
    
    if(any(duplicated(samples))) {
        stop(
            "Sample names are not unique:  \n  ",
            paste(files, collapse = "\n  "),
            call. = FALSE
        )
    }
    
    return(samples)
}


read_sample_files <- function(
        path = ".", pattern = "\\.counts$", delim = "\t", ...
    ) {
    #  rdrr.io/github/HuntsmanCancerInstitute/hciR/src/R/read_sample_files.R
    
    #' Read and combine common output files
    #'
    #' Read output files and parse the file name to add sample IDs in the first column
    #'
    #' @param path the path to output files, the default corresponds to the working directory.
    #' @param pattern regular expression for file name matching
    #' @param delim separator in output files, default table
    #' @param \dots additional options such as col_names passed to \code{read_delim}.
    #'
    #' @note Sample names are parsed from file names without extensions.  If the file name is not unique,
    #'     the parent directory is used.   Requires tibble > version 1.2 to avoid error in add_column
    #'
    #' @return A list with coverage and stats data.frames
    #'
    #' @author Chris Stubben
    #'
    #' @examples
    #' \dontrun{
    #'     #FeatureCounts summary (second column name with *.bam is always unique, so skip and assign)
    #'     fc <- read_sample_files(".summary$", skip = 1, col_names = c("status", "count"))
    #'     filter(fc, count! = 0) %>%
    #'     hchart("bar", x = sample, y = count, group = status) %>%
    #'     hc_plotOptions(bar = list(stacking = "normal"))
    #' }
    #' @export
    
    outfiles <- list.files(path, pattern, recursive = TRUE, full.names = TRUE)
    if(length(outfiles) == 0) {
        stop("No ", pattern, " files found in ", path, call. = FALSE)
    }
    samples <- extract_samples(outfiles)
    
    out1 <- vector("list", length(outfiles))
    for(i in seq_along(outfiles)){
        message("Reading ", outfiles[i])
            x <- suppressMessages(
               readr::read_delim(outfiles[i], delim = delim, ...)
            )
            # Requires tibble > 1.2 (to add single name into column)
            out1[[i]] <- tibble::add_column (
                x, sample = samples[i], .before = 1
            )
        }
        dplyr::bind_rows(out1)
}


read_STAR <- function(path = ".", pattern, reshape = FALSE) {
    #  rdrr.io/github/HuntsmanCancerInstitute/hciR/src/R/read_STAR.R
    
    #' Read STAR log files
    #'
    #' Read STAR Log.final.out files and optionally reshape into wide format.
    #'
    #' @param path the path to STAR log files, the default corresponds to the working directory.
    #' @param pattern regular expression for file name matching, default .final.out
    #' @param reshape reshape percent mapping into wide format with samples in rows
    #'
    #' @note Reading output files requires a unique sample identifier in either the file name or parent directory
    #'
    #' @return A tibble
    #'
    #' @author Chris Stubben
    #'
    #' @examples
    #' \dontrun{
    #' library(hciRdata)
    #' res <- results_all(pasilla$dds, fly98)
    #' res
    #' }
    #' @export
    
    #  For testing...
    # path <- "."
    # pattern <- "\\.final.out$"
    # reshape <- FALSE
    
    if(missing(pattern)) pattern <- "\\.final.out$"
    #  Suppress warnings about subheaders... Warning: 4 parsing failures.
    x <- suppressWarnings(
        read_sample_files(
            path, pattern, col_names = c("stat", "value"),
            skip = 5,
            trim_ws = TRUE
        )
    )

    x <- dplyr::filter(x, !is.na(value)) %>% # drop subheader with NA value, e.g., 'UNIQUE READS:'
        dplyr::mutate(
            stat = gsub(" \\|$", "", stat),  # remove pipe from end of statistic
            value = gsub("%", "", value),  # drop % from value
            value = as.numeric(value)  # change value to numeric
        )

    #  Add unmapped reads
    columns <- c(
        "Number of input reads",
        "Average input read length",
        "Uniquely mapped reads number",
        "Uniquely mapped reads %",
        "Average mapped length",
        "Number of splices: Total",
        "Number of splices: Annotated (sjdb)",
        "Number of splices: GT/AG",
        "Number of splices: GC/AG",
        "Number of splices: AT/AC",
        "Number of splices: Non-canonical",
        "Mismatch rate per base, %",
        "Deletion rate per base",
        "Deletion average length",
        "Insertion rate per base",
        "Insertion average length",
        "Number of reads mapped to multiple loci",
        "% of reads mapped to multiple loci",
        "Number of reads mapped to too many loci",
        "% of reads mapped to too many loci",
        "Number of reads unmapped: too many mismatches",
        "% of reads unmapped: too many mismatches",
        "Number of reads unmapped: too short",
        "% of reads unmapped: too short",
        "Number of reads unmapped: other",
        "% of reads unmapped: other"
    )

    # #  See stackoverflow.com/questions/40749742/add-missing-subtotals-to-each-group-using-dplyr
    # x <- dplyr::group_by(x, sample) %>%
    #     dplyr::summarize(
    #         value = value[stat == 'Number of input reads'] -
    #             sum(value[stat %in% mapped]),
    #         stat = 'Unmapped reads'
    #     ) %>%
    #     dplyr::bind_rows(x) %>%
    #     dplyr::select(sample, stat, value) %>%  # back to original column order
    #     dplyr::arrange(sample)

    #  Split unmapped reads into too many mismatches, too short and other?
    #  Divide %unmapped too short by total %unmapped and muliply by Unmapped reads to get estimate
    
    if(isTRUE(reshape)){
        # n <- c(mapped, "Unmapped reads", "Number of input reads")
        n <- columns
        x <- dplyr::filter(x, stat %in% n) %>%
            dplyr::mutate(stat = factor(stat, levels = n)) %>%
            tidyr::spread(stat, value) %>%
            dplyr::mutate_each(dplyr::funs(as.integer), -1)
    }
    
    return(x)
}


#  Set up custom ggplot2 plot themes ------------------------------------------
theme_slick <- theme_classic() +
    theme(
        panel.grid.major = element_line(size = 0.4),
        panel.grid.minor = element_line(size = 0.2),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.4),
        axis.text = element_text(color = "black"),
        axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown(),
        plot.title = ggtext::element_markdown(),
        legend.title = element_text(size = 0)
    )

theme_slick_no_legend <- theme_slick + theme(legend.position = "none")

theme_slick_no_boundary <- theme_classic() +
    theme(
        panel.grid.major = element_line(size = 0.4),
        panel.grid.minor = element_line(size = 0.2),
        axis.line = element_line(size = 0),
        axis.ticks = element_line(size = 0),
        axis.text = element_text(color = "black"),
        axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown(),
        plot.title = ggtext::element_markdown(),
        legend.title = element_text(size = 0)
    )


# #  Declare options, flags, types, default values, etc. ------------------------
# script <- "munge_STAR_log_files.R"
# 
# #  Create a parser
# ap <- arg_parser(
#     name = script,
#     description = "
#         This script takes ... and ... .
#     ",
#     hide.opts = TRUE
# )
# 
# #  Add command line arguments
# ap <- callMetadataArguments()
# ap <- add_argument(
#     ap,
#     short = "-t",
#     arg = "--TbyC",
#     type = "character",
#     default = NA,
#     nargs = 1,
#     help = "
#         Path and .rds for the table of TAD boundary overlaps stratified by
#         compartments transitions between two samples <chr>
#     "
# )
# ap <- add_argument(
#     ap,
#     short = "-s",
#     arg = "--deseq2_TAD",
#     type = "character",
#     default = NA,
#     nargs = 1,
#     help = "
#         Path and .rds for a table of DESeq2 information with respect to TAD
#         boundary overlap information <chr>
#     "
# )
# ap <- add_argument(
#     ap,
#     short = "-j",
#     arg = "--joint",
#     type = "character",
#     default = 4,
#     nargs = 1,
#     help = "
#         Handle boundaries, including 'joint' (i.e., 'shared') boundaries
#         through either three or four categories; if '3', then
#         'shared.sample_1', 'shared.sample_2', and/or 'shared.both' are combined
#         into 'shared': the categories are 'sample_1', 'sample_2', and 'shared';
#         if '4', then the categories are 'sample_1', 'sample_2',
#         'shared.sample_1', and 'shared.sample_2'; if '5', then the categories
#         are 'sample_1', 'sample_2', 'shared.both', 'shared.sample_1', and
#         'shared.sample_2' <int = 3 | int = 4 | int = 5>
#     "
# )
# ap <- add_argument(
#     ap,
#     short = "-a",
#     arg = "--deseq2_numerator",
#     type = "character",
#     default = NA,
#     nargs = 1,
#     help = "String for the 'numerator sample' in the DESeq2 analysis <chr>"
# )
# ap <- add_argument(
#     ap,
#     short = "-b",
#     arg = "--deseq2_denominator",
#     type = "character",
#     default = NA,
#     nargs = 1,
#     help = "String for the 'denominator sample' in the DESeq2 analysis <chr>"
# )
# ap <- add_argument(
#     ap,
#     short = "-c",
#     arg = "--color",
#     type = "character",
#     default = NA,
#     help = "Hex code color for ggplot2 geom_col() plots <chr>"
# )
# ap <- add_argument(
#     ap,
#     short = "-y",
#     arg = "--ylim",
#     type = "integer",
#     default = 4,
#     help = "ylim max and min for ggplot2 geom_col() plots <int>"
# )
# ap <- add_argument(
#     ap,
#     short = "-r",
#     arg = "--out_compress_rds",
#     type = "logical",
#     default = FALSE,
#     nargs = 1,
#     help = "Compress output .rds or not <logical>"
# )


#  TBD
setwd("/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/Berletch_work")
list.files()

x.reshape.no <- read_STAR(pattern = ".final.out$", reshape = FALSE)
x.reshape.yes <- read_STAR(pattern = ".final.out$", reshape = TRUE)
# View(x.reshape.no)

y <- x.reshape.no$sample %>%
    stringr::str_split("\\.", simplify = TRUE) %>%
    tibble::as_tibble(.name_repair = "unique") %>%
    dplyr::select(c(`...1`, `...3`)) %>%
    dplyr::rename(prefix = `...1`, suffix = `...3`) %>%
    tidyr::unite(sample, c(prefix, suffix), sep = ".")
x.reshape.no <- dplyr::bind_cols(x.reshape.no, y) %>%
    dplyr::select(-`sample...1`) %>%
    dplyr::rename(sample = `sample...4`) %>%
    dplyr::relocate(sample, .before = stat)
rm(y)

# entries <- c(
#     "Uniquely mapped reads number",
#     "Number of reads mapped to multiple loci",
#     "Number of reads mapped to too many loci",
#     "Number of reads unmapped: too many mismatches",
#     "Number of reads unmapped: too short",
#     "Number of reads unmapped: other"
# )
entries <- c(
    "Number of reads unmapped: other",
    "Number of reads unmapped: too short",
    "Number of reads unmapped: too many mismatches",
    "Number of reads mapped to too many loci",
    "Number of reads mapped to multiple loci",
    "Uniquely mapped reads number"
)
x.reshape.no.filt <- x.reshape.no[x.reshape.no$stat %in% entries, ]
x.reshape.no.filt$stat <- x.reshape.no.filt$stat %>%
    factor(., levels = entries)

ggplot2::ggplot(x.reshape.no.filt, aes(x = sample, y = value, fill = stat)) +
    geom_bar(position = "stack", stat = "identity") +
    theme_slick +
    ggplot2::theme(
        axis.text.x = element_text(angle = 30, hjust = 1)
    )

# ggplot2::ggplot(x, aes(x = sample, y = `Number of input reads`)) +
#     geom_bar(position = "dodge", stat = "identity") +
#     theme_slick_no_boundary +
#     ggplot2::theme(
#         axis.text.x = element_text(angle = 30, hjust = 1)
#     )
