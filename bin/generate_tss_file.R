library(rtracklayer)
library(argparse)

# Avoid scientific notation
options(scipen=999)
options(stringsAsFactors = FALSE)

parser = argparse::ArgumentParser(description='Makes a BED file of TSSs provided a GENCODE GTF file.')
parser$add_argument('input_file', help='GENCODE GTF file.')
parser$add_argument('output_file', help='Output TSS BED file.')
args = parser$parse_args()

get_tss <- function(gene_ann) {
  gene_ann = subset(gene_ann, type == "transcript")
  pos <- subset(gene_ann, strand == "+")
  neg <- subset(gene_ann, strand == "-")
  pos <- pos[order(pos[,"start"]),]
  neg <- neg[order(neg[,"start"], decreasing = TRUE),]
  pos <- pos[!duplicated(pos$transcript_id),]
  neg <- neg[!duplicated(neg$transcript_id),]

  # Output BED format
  pos[,"end"] <- pos[,"start"]
  pos[,"start"] <- pos["start"] - 1

  neg[,"start"] <- neg[,"end"] - 1
  neg[,"end"] <- neg[,"start"] + 1
  return(rbind(pos, neg))
}

tss = get_tss(rtracklayer::readGFF(args$input_file))

if ('transcript_support_level' %in% colnames(tss)) {
  tss = subset(tss, is.na(transcript_support_level) | transcript_support_level >= 2)
}
tss$score = '.'
write.table(tss[, c('seqid', 'start', 'end', 'gene_id', 'score', 'strand', 'gene_name', 'transcript_id', 'transcript_type')], quote=F, row.names=F, col.names=F, sep='\t', file=args$output_file)
