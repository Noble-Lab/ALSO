#!/bin/Rscript

#  plot-tally-assignments_bulk-RNA-seq.R
#  KA

library(plyr)
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
dir_data <- "data/Berletch_work"
dir_path <- paste(dir_projects, dir_base, dir_data, sep = "/")
tally <- paste(dir_path, "scoring-ALSO-SNPsplit_page-1.csv", sep = "/")

tally <- readr::read_csv(tally, show_col_types = FALSE)
tally.init <- tally
# tally <- tally.init
tally <- tally %>%
    dplyr::select(-contains(
        c("gte", "complement", "intersection", "method", "rep")
    ))

tally.detailed <- tally
tally.detailed[["mm11.only_in_mm11"]] <- tally.detailed[["mm11.mm11"]] - tally.detailed[["mm11.SPRET"]]
tally.detailed[["SPRET.only_in_SPRET"]] <- tally.detailed[["SPRET.SPRET"]] - tally.detailed[["SPRET.mm11"]]
tally.detailed <- tally.detailed %>%
    dplyr::rename(
        "mm11.also_in_SPRET" = "mm11.SPRET",
        "SPRET.also_in_mm11" = "SPRET.mm11",
        "ambiguous" = "mm11.ambiguous"
    ) %>%
    dplyr::select(-c("mm11.mm11", "SPRET.SPRET", "SPRET.ambiguous"))

#  Pivot to format appropriate for bar charts: "Tally detailed"
tally.detailed <- tally.detailed %>%
    tidyr::pivot_longer(!sample, names_to = "assignment", values_to = "reads")

tally.detailed$sample <- factor(tally.detailed$sample)
tally.detailed$assignment <- factor(tally.detailed$assignment)
tally.detailed$assignment <- factor(
    tally.detailed$assignment, levels = c(
        "SPRET.only_in_SPRET", "SPRET.also_in_mm11", "mm11.only_in_mm11",
        "mm11.also_in_SPRET", "ambiguous"
    )
)

tally.detailed <- plyr::ddply(
    tally.detailed, .(sample), transform, percent = reads/sum(reads) * 100
)
tally.detailed <- ddply(tally.detailed, .(sample), transform, pos = (cumsum(reads) - 0.5 * reads))
tally.detailed$label <- paste0(sprintf("%.0f", tally.detailed$percent), "%")

tally <- tally %>%
    dplyr::rename("ambiguous" = "mm11.ambiguous") %>%
    dplyr::select(-c("SPRET.ambiguous", "mm11.SPRET", "SPRET.mm11")) %>%
    dplyr::rename("mm11" = "mm11.mm11", "SPRET" = "SPRET.SPRET") %>%
    tidyr::pivot_longer(!sample, names_to = "assignment", values_to = "reads")

tally$sample <- factor(tally$sample)
tally$assignment <- factor(tally$assignment)
tally$assignment <- factor(
    tally$assignment, levels = c("SPRET", "mm11", "ambiguous")
)

tally <- plyr::ddply(
    tally, .(sample), transform, percent = reads/sum(reads) * 100
)
tally <- ddply(tally, .(sample), transform, pos = (cumsum(reads) - 0.33 * reads))
tally$label <- paste0(sprintf("%.0f", tally$percent), "%")


#  "Tally"
bar_count <- tally %>%
    ggplot2::ggplot(aes(x = sample, y = reads, fill = assignment)) +
    geom_bar(position = "stack", stat = "identity") +
    # geom_text(aes(y = pos, label = label), size = 3, position = "stack") +
    ggplot2::ggtitle("Reads stratified by assignment") +
    theme_slick +
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
bar_count

bar_pct <- tally %>%
    ggplot2::ggplot(aes(x = sample, y = reads, fill = assignment)) +
    geom_bar(position = "fill", stat = "identity") +
    geom_text(aes(y = pos, label = label), size = 3, position = "fill") +
    ggplot2::ggtitle("Reads stratified by assignment") +
    theme_slick +
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
bar_pct

#  "Tally detailed"
bar_count_detailed <- tally.detailed %>%
    ggplot2::ggplot(aes(x = sample, y = reads, fill = assignment)) +
    geom_bar(position = "stack", stat = "identity") +
    # geom_text(aes(y = pos, label = label), size = 3, position = "stack") +
    ggplot2::ggtitle("Reads stratified by assignment") +
    theme_slick +
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
bar_count_detailed

bar_pct_detailed <- tally.detailed %>%
    ggplot2::ggplot(aes(x = sample, y = reads, fill = assignment)) +
    geom_bar(position = "fill", stat = "identity") +
    geom_text(aes(y = pos, label = label), size = 3, position = "fill") +
    ggplot2::ggtitle("Reads stratified by assignment") +
    theme_slick +
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
bar_pct_detailed

ggplot2::ggsave(
    paste(
        dir_projects, dir_base, dir_data, "tally_bar_count.pdf",
        sep = "/"
    ),
    plot = bar_count
)
ggplot2::ggsave(
    paste(
        dir_projects, dir_base, dir_data, "tally_bar_pct.pdf",
        sep = "/"
    ),
    plot = bar_pct
)
ggplot2::ggsave(
    paste(
        dir_projects, dir_base, dir_data, "tally_bar_count_detailed.pdf",
        sep = "/"
    ),
    plot = bar_count_detailed
)
ggplot2::ggsave(
    paste(
        dir_projects, dir_base, dir_data, "tally_bar_pct_detailed.pdf",
        sep = "/"
    ),
    plot = bar_pct_detailed
)
