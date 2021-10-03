#!/usr/bin/env Rscript

#  Chengxiang (CX) Qiu
#  Adapted, edited by KA

library(rtracklayer)
library(argparser)

#  Avoid scientific notation
options(scipen = 999)
options(stringsAsFactors = FALSE)

#  Create a parser, add command-line arguments
ap <- arg_parser(
    "This script makes a gene-body .bed file from a GENCODE/Ensembl .gtf file."
)
ap <- add_argument(
    ap,
    "--infile",
    type = "character",
    help = "GENCODE/Ensembl .gtf file, including path <chr>"
)
ap <- add_argument(
    ap,
    "--outfile",
    type = "character",
    help = "Output gene-body .bed file, including path <chr>"
)

#  Parse the arguments
cl <- c(  # Use for interactive testing
    "--infile", "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/doc/Mus_musculus_129s1svimj.129S1_SvImJ_v1.104.chr.gtf.gz",
    "--outfile", "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/doc/129S1-SvImJ.gene-body.bed"
)
arguments <- parse_args(ap, cl)  # Use for interactive testing
# arguments <- parse_args(ap)  # Use for command-line calls

get_gene_body <- function(gene_annotation) {
    gene_annotation <- subset(gene_annotation, type == "gene")
    gene_annotation <- gene_annotation[!duplicated(gene_annotation$gene_id), ]

    #  Output .bed format
    gene_annotation[, "start"] <- gene_annotation[, "start"] - 1

    return(gene_annotation)
}

gene_bodies <- get_gene_body(rtracklayer::readGFF(arguments$infile))
gene_bodies$score <- '.'
write.table(
    gene_bodies[, c(
        'seqid',
        'start',
        'end',
        'gene_id',
        'score',
        'strand',
        'gene_name',
        'gene_biotype'
    )],
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = '\t',
    file = arguments$outfile
)
