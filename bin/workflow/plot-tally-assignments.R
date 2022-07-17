#!/bin/Rscript

#  plot-tally-assignments.R
#  KA

library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)

theme_slick <- theme_classic() +
    theme(
        panel.grid.major = ggplot2::element_line(size = 0.4),
        panel.grid.minor = ggplot2::element_line(size = 0.2),
        axis.line = ggplot2::element_line(size = 0.2),
        axis.ticks = ggplot2::element_line(size = 0.4),
        axis.text = ggplot2::element_text(color = "black"),
        # axis.title.x = ggtext::element_markdown(),
        # axis.title.y = ggtext::element_markdown(),
        # plot.title = ggtext::element_markdown(),
        axis.title.x = ggplot2::element_text(),
        axis.title.y = ggplot2::element_text(),
        plot.title = ggplot2::element_text(),
        legend.title = ggplot2::element_text(size = 0)
    )

theme_slick_no_legend <- theme_slick + theme(legend.position = "none")


dir_projects <- "/Users/kalavattam/Dropbox/UW/projects-etc"
dir_base <- "2021_kga0_4dn-mouse-cross"
dir_results <- "results/kga0/2022-0627_tally-assignments"
dir_path <- paste(dir_projects, dir_base, dir_results, sep = "/")
tally <- paste(dir_path, list.files(dir_path)[2], sep = "/")

tally <- readr::read_delim(tally)
tally$pct_sample_1 <- tally$sample_1/tally$total
tally$pct_sample_2 <- tally$sample_2/tally$total
tally$pct_ambiguous <- tally$ambiguous/tally$total

tally$exp <- tally$exp %>%
    gsub("Disteche_", "", .) %>%
    gsub("_ALSO", "", .)

tally$exp <- tally$exp %>%
    factor(., levels = gtools::mixedsort(tally$exp))


bar_count <- tally %>%
    dplyr::select(-c(pct_sample_1, pct_sample_2, pct_ambiguous, total)) %>%
    dplyr::rename(c(
        "mm10" = "sample_1",
        "CAST" = "sample_2",
        "ambiguous" = "ambiguous",
        "experiment" = "exp"
    )) %>%
    tidyr::pivot_longer(!experiment, names_to = "assignment", values_to = "read pairs") %>%
    ggplot2::ggplot(aes(x = experiment, y = `read pairs`, fill = assignment)) +
    geom_bar(position = "stack", stat = "identity") +
    ggplot2::ggtitle("Read pairs stratified by assignment") +
    theme_slick +
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

bar_pct <- tally %>%
    dplyr::select(-c(sample_1, sample_2, ambiguous, total)) %>%
    dplyr::rename(c(
        "mm10" = "pct_sample_1",
        "CAST" = "pct_sample_2",
        "ambiguous" = "pct_ambiguous",
        "experiment" = "exp"
    )) %>%
    tidyr::pivot_longer(!experiment, names_to = "assignment", values_to = "fraction") %>%
    ggplot2::ggplot(aes(x = experiment, y = fraction, fill = assignment)) +
    ggplot2::geom_bar(position = "fill", stat = "identity") +
    ggplot2::ggtitle("Fractions of read pairs stratified by assignment") +
    theme_slick +
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )


bar_count
bar_pct

ggplot2::ggsave(
    paste(
        dir_projects, dir_base, dir_results, "tally_bar_count.pdf",
        sep = "/"
    ),
    plot = bar_count
)
ggplot2::ggsave(
    paste(
        dir_projects, dir_base, dir_results, "tally_bar_pct.pdf",
        sep = "/"
    ),
    plot = bar_pct
)

tally %>% print(n = 25)
