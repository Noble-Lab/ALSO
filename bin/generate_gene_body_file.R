#!/usr/bin/env Rscript

#  Chengxiang (CX) Qiu
#  Adapted, edited by KA

library(rtracklayer)
library(argparse)

#  Avoid scientific notation
options(scipen = 999)
options(stringsAsFactors = FALSE)

parser <- argparse::ArgumentParser(description = 'Makes a .bed file of gene bodies provided a GENCODE/Ensembl .gtf file.')
parser$add_argument('input_file', help = 'GENCODE/Ensembl .gtf file.')
parser$add_argument('output_file', help = 'Output gene body .bed file.')
args <- parser$parse_args()

get_gene_body <- function(gene_ann) {
    gene_ann <- subset(gene_ann, type == "gene")
    gene_ann <- gene_ann[!duplicated(gene_ann$gene_id), ]

    #  Output .bed format
    gene_ann[, "start"] <- gene_ann[, "start"] - 1

    return(gene_ann)
}

gene_bodies <- get_gene_body(rtracklayer::readGFF(args$input_file))
gene_bodies$score <- '.'
write.table(
    gene_bodies[,
    c(
        'seqid',
        'start',
        'end',
        'gene_id',
        'score',
        'strand',
        'gene_name',
        'gene_type'
    )],
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = '\t',
    file = args$output_file
)
