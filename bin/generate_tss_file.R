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
    "This script makes a transcription-start-site .bed file from a
    GENCODE/Ensembl .gtf file."
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
    help = "Output TSS .bed file, including path <chr>"
)

#  Parse the arguments
# cl <- c(  # Use for interactive testing
#     "--infile", "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/doc/Mus_musculus_129s1svimj.129S1_SvImJ_v1.104.chr.gtf.gz",
#     "--outfile", "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/doc/129S1-SvImJ.TSS.bed"
# )
# arguments <- parse_args(ap, cl)  # Use for interactive testing
arguments <- parse_args(ap)  # Use for command-line calls

get_TSS <- function(gene_annotation) {
    #TODO Documentation, etc.
    gene_annotation <- subset(gene_annotation, type == "transcript")
    pos <- subset(gene_annotation, strand == "+")
    neg <- subset(gene_annotation, strand == "-")
    pos <- pos[order(pos[, "start"]), ]
    neg <- neg[order(neg[, "start"], decreasing = TRUE), ]
    pos <- pos[!duplicated(pos$transcript_id), ]
    neg <- neg[!duplicated(neg$transcript_id), ]

    #  Output .bed format
    pos[, "end"] <- pos[, "start"]
    pos[, "start"] <- pos["start"] - 1

    neg[, "start"] <- neg[, "end"] - 1
    neg[, "end"] <- neg[, "start"] + 1
    return(rbind(pos, neg))
}

TSS <- get_TSS(rtracklayer::readGFF(arguments$infile))

if ('transcript_support_level' %in% colnames(TSS)) {
    TSS <- subset(
        TSS,
        is.na(transcript_support_level) | transcript_support_level >= 2
    )
}
TSS$score <- '.'
write.table(
    TSS[, c(
        'seqid',
        'start',
        'end',
        'gene_id',
        'score',
        'strand',
        'gene_name',
        'transcript_id',
        'transcript_biotype'
    )],
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = '\t',
    file = arguments$outfile
)
