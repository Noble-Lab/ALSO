library(rtracklayer)
library(argparse)

# Avoid scientific notation
options(scipen=999)
options(stringsAsFactors = FALSE)

parser = argparse::ArgumentParser(description='Makes a BED file of TSSs provided a GENCODE GTF file.')
parser$add_argument('input_file', help='GENCODE GTF file.')
parser$add_argument('output_file', help='Output gene body BED file.')
args = parser$parse_args()

get_gene_body <- function(gene_ann) {
  gene_ann = subset(gene_ann, type == "gene")
  gene_ann <- gene_ann[!duplicated(gene_ann$gene_id),]

  # Output BED format
  gene_ann[,"start"] <- gene_ann[,"start"] - 1

  return(gene_ann)
}

gene_bodies = get_gene_body(rtracklayer::readGFF(args$input_file))
gene_bodies$score = '.'
write.table(gene_bodies[, c('seqid', 'start', 'end', 'gene_id', 'score', 'strand', 'gene_name', 'gene_type')], quote=F, row.names=F, col.names=F, sep='\t', file=args$output_file)
