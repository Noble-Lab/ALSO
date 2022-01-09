
library(AnnotationDbi)
library(BSgenome)
library(GenomicFeatures)
library(Mmusculus129S1SvImJ)
library(MmusculusCASTEiJ)
library(Mmusculus129Inserted)
library(MmusculusCASTInserted)
library(MmusculusCAST129Inserted)
library(Mmusculus129Nmasked)
library(MmusculusCASTNmasked)
library(MmusculusCAST129Nmasked)
library(tidyverse)
library(reticulate)


# -----------------
genome_129_BS <- Mmusculus129S1SvImJ
genome_CAST_BS <- MmusculusCASTEiJ
genome_129_inserted_BS <- Mmusculus129Inserted
genome_CAST_inserted_BS <- MmusculusCASTInserted
genome_CAST_129_inserted_BS <- MmusculusCAST129Inserted
genome_129_Nmasked_BS <- Mmusculus129Nmasked
genome_CAST_Nmasked_BS <- MmusculusCASTNmasked
genome_CAST_129_Nmasked_BS <- MmusculusCAST129Nmasked


# -----------------
# genome_CAST_GR <- GenomeInfoDb::seqinfo(MmusculusCASTEiJ) %>%
#     GenomicRanges::GRanges()
# genome_CAST_GR

# genome_CAST_BS <- MmusculusCASTEiJ
# genome_CAST_BS
# genome_CAST_BS %>% seqnames()
# genome_CAST_BS %>% organism()
# genome_CAST_BS %>% provider()


# genome_CASTNmasked_GR <- GenomeInfoDb::seqinfo(MmusculusCASTNmasked) %>%
#     GenomicRanges::GRanges()
# genome_CASTNmasked_GR

# genome_CASTNmasked_BS <- MmusculusCASTNmasked
# genome_CASTNmasked_BS
# genome_CASTNmasked_BS %>% seqnames()
# genome_CASTNmasked_BS %>% organism()
# genome_CASTNmasked_BS %>% provider()


# genome_CASTinserted_GR <- GenomeInfoDb::seqinfo(MmusculusCASTInserted) %>%
#     GenomicRanges::GRanges()
# genome_CASTinserted_GR

# genome_CASTinserted_BS <- MmusculusCASTInserted
# genome_CASTinserted_BS
# genome_CASTinserted_BS %>% seqnames()
# genome_CASTinserted_BS %>% organism()
# genome_CASTinserted_BS %>% provider()


# -----------------
chr = "chr1"
start = 1
end = 60
genome_129_BS %>% getSeq(., chr, start = start, end = end)
genome_CAST_BS %>% getSeq(., chr, start = start, end = end)
genome_129_inserted_BS %>% getSeq(., chr, start = start, end = end)
genome_CAST_inserted_BS %>% getSeq(., chr, start = start, end = end)
genome_CAST_129_inserted_BS %>% getSeq(., chr, start = start, end = end)
genome_129_Nmasked_BS %>% getSeq(., chr, start = start, end = end)
genome_CAST_Nmasked_BS %>% getSeq(., chr, start = start, end = end)
genome_CAST_129_Nmasked_BS %>% getSeq(., chr, start = start, end = end)

chr = "chr1"
start = 501
end = 600
genome_129_BS %>% getSeq(., chr, start = start, end = end)
genome_CAST_BS %>% getSeq(., chr, start = start, end = end)
genome_129_inserted_BS %>% getSeq(., chr, start = start, end = end)
genome_CAST_inserted_BS %>% getSeq(., chr, start = start, end = end)
genome_CAST_129_inserted_BS %>% getSeq(., chr, start = start, end = end)
genome_129_Nmasked_BS %>% getSeq(., chr, start = start, end = end)
genome_CAST_Nmasked_BS %>% getSeq(., chr, start = start, end = end)
genome_CAST_129_Nmasked_BS %>% getSeq(., chr, start = start, end = end)

chr = "chr4"
start = 24999981
end = 25000050
genome_129_BS %>% getSeq(., chr, start = start, end = end)
genome_CAST_BS %>% getSeq(., chr, start = start, end = end)
genome_129_inserted_BS %>% getSeq(., chr, start = start, end = end)
genome_CAST_inserted_BS %>% getSeq(., chr, start = start, end = end)
genome_CAST_129_inserted_BS %>% getSeq(., chr, start = start, end = end)
genome_129_Nmasked_BS %>% getSeq(., chr, start = start, end = end)
genome_CAST_Nmasked_BS %>% getSeq(., chr, start = start, end = end)
genome_CAST_129_Nmasked_BS %>% getSeq(., chr, start = start, end = end)


# -----------------
chain_129_to_mm10 <- paste0("./liftOver/", "129S1-SvImJ-to-mm10.over.chain.munged") %>%
    rtracklayer::import.chain()

chain_mm10_to_129 <- paste0("./liftOver/", "mm10-to-129S1-SvImJ.over.chain.munged") %>%
    rtracklayer::import.chain()

chain_CAST_to_mm10 <- paste0("./liftOver/", "CAST-EiJ-to-mm10.over.chain.munged") %>%
    rtracklayer::import.chain()

chain_mm10_to_CAST <- paste0("./liftOver/", "mm10-to-CAST-EiJ.over.chain.munged") %>%
    rtracklayer::import.chain()


# -----------------
chr = "chr11"
start = 101442663
end = start + 19
gr <- GRanges(
    seqnames = Rle(chr),
    ranges = IRanges(start = start, end = end)
)
gr_liftOver_129_to_mm10 <- rtracklayer::liftOver(gr, chain_129_to_mm10) %>% unlist()
gr_liftOver_CAST_to_mm10 <- rtracklayer::liftOver(gr, chain_CAST_to_mm10) %>% unlist()
gr_liftOver_mm10_to_129 <- rtracklayer::liftOver(gr, chain_mm10_to_129) %>% unlist()
gr_liftOver_mm10_to_CAST <- rtracklayer::liftOver(gr, chain_mm10_to_CAST) %>% unlist()

# -----------------
genome_129_BS %>% getSeq(., names = gr_liftOver_129_to_mm10) %>% unlist()
genome_129_Nmasked_BS %>% getSeq(., names = gr) %>% unlist()
genome_129_inserted_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_129_inserted_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_129_Nmasked_BS %>% getSeq(., names = gr) %>% unlist()

#  Go with this.
genome_129_BS %>% getSeq(., names = gr) %>% unlist()
genome_129_Nmasked_BS %>% getSeq(., names = gr_liftOver_129_to_mm10) %>% unlist()
genome_129_inserted_BS %>% getSeq(., names = gr_liftOver_129_to_mm10) %>% unlist()
genome_CAST_129_inserted_BS %>% getSeq(., names = gr_liftOver_129_to_mm10) %>% unlist()
genome_CAST_129_Nmasked_BS %>% getSeq(., names = gr_liftOver_129_to_mm10) %>% unlist()

# -----------------
genome_CAST_BS %>% getSeq(., names = gr_liftOver_CAST_to_mm10) %>% unlist()
genome_CAST_Nmasked_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_inserted_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_129_Nmasked_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_129_inserted_BS %>% getSeq(., names = gr) %>% unlist()

#  Go with this.
genome_CAST_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_Nmasked_BS %>% getSeq(., names = gr_liftOver_CAST_to_mm10) %>% unlist()
genome_CAST_inserted_BS %>% getSeq(., names = gr_liftOver_CAST_to_mm10) %>% unlist()
genome_CAST_129_Nmasked_BS %>% getSeq(., names = gr_liftOver_CAST_to_mm10) %>% unlist()
genome_CAST_129_inserted_BS %>% getSeq(., names = gr_liftOver_CAST_to_mm10) %>% unlist()


# -----------------
chr = "chr11"
start = 101442663
end = start + 19
gr <- GRanges(
    seqnames = Rle(chr),
    ranges = IRanges(start = start, end = end)
)
gr_liftOver_129_to_mm10 <- rtracklayer::liftOver(gr, chain_129_to_mm10) %>% unlist()
gr_liftOver_CAST_to_mm10 <- rtracklayer::liftOver(gr, chain_CAST_to_mm10) %>% unlist()
gr_liftOver_mm10_to_129 <- rtracklayer::liftOver(gr, chain_mm10_to_129) %>% unlist()
gr_liftOver_mm10_to_CAST <- rtracklayer::liftOver(gr, chain_mm10_to_CAST) %>% unlist()

# -----------------
genome_129_BS %>% getSeq(., names = gr_liftOver_129_to_mm10) %>% unlist()
genome_129_Nmasked_BS %>% getSeq(., names = gr) %>% unlist()
genome_129_inserted_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_129_inserted_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_129_Nmasked_BS %>% getSeq(., names = gr) %>% unlist()

#  Go with this.
genome_129_BS %>% getSeq(., names = gr) %>% unlist()
genome_129_Nmasked_BS %>% getSeq(., names = gr_liftOver_129_to_mm10) %>% unlist()
genome_129_inserted_BS %>% getSeq(., names = gr_liftOver_129_to_mm10) %>% unlist()
genome_CAST_129_inserted_BS %>% getSeq(., names = gr_liftOver_129_to_mm10) %>% unlist()
genome_CAST_129_Nmasked_BS %>% getSeq(., names = gr_liftOver_129_to_mm10) %>% unlist()

# -----------------
genome_CAST_BS %>% getSeq(., names = gr_liftOver_CAST_to_mm10) %>% unlist()
genome_CAST_Nmasked_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_inserted_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_129_Nmasked_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_129_inserted_BS %>% getSeq(., names = gr) %>% unlist()

#  Go with this.
genome_CAST_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_Nmasked_BS %>% getSeq(., names = gr_liftOver_CAST_to_mm10) %>% unlist()
genome_CAST_inserted_BS %>% getSeq(., names = gr_liftOver_CAST_to_mm10) %>% unlist()
genome_CAST_129_Nmasked_BS %>% getSeq(., names = gr_liftOver_CAST_to_mm10) %>% unlist()
genome_CAST_129_inserted_BS %>% getSeq(., names = gr_liftOver_CAST_to_mm10) %>% unlist()


# -----------------
chr = "chr19"
start = 10000000
end = start + 59

gr <- GRanges(
    seqnames = Rle(chr),
    ranges = IRanges(start = start, end = end)
)
gr_liftOver_129_to_mm10 <- rtracklayer::liftOver(gr, chain_129_to_mm10) %>% unlist()
gr_liftOver_CAST_to_mm10 <- rtracklayer::liftOver(gr, chain_CAST_to_mm10) %>% unlist()
gr_liftOver_mm10_to_129 <- rtracklayer::liftOver(gr, chain_mm10_to_129) %>% unlist()
gr_liftOver_mm10_to_CAST <- rtracklayer::liftOver(gr, chain_mm10_to_CAST) %>% unlist()

# -----------------
genome_129_BS %>% getSeq(., names = gr_liftOver_129_to_mm10) %>% unlist()
genome_129_Nmasked_BS %>% getSeq(., names = gr) %>% unlist()
genome_129_inserted_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_129_inserted_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_129_Nmasked_BS %>% getSeq(., names = gr) %>% unlist()

#  Go with this.
genome_129_BS %>% getSeq(., names = gr) %>% unlist()
genome_129_Nmasked_BS %>% getSeq(., names = gr_liftOver_129_to_mm10) %>% unlist()
genome_129_inserted_BS %>% getSeq(., names = gr_liftOver_129_to_mm10) %>% unlist()
genome_CAST_129_inserted_BS %>% getSeq(., names = gr_liftOver_129_to_mm10) %>% unlist()
genome_CAST_129_Nmasked_BS %>% getSeq(., names = gr_liftOver_129_to_mm10) %>% unlist()

# -----------------
genome_CAST_BS %>% getSeq(., names = gr_liftOver_CAST_to_mm10) %>% unlist()
genome_CAST_Nmasked_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_inserted_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_129_Nmasked_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_129_inserted_BS %>% getSeq(., names = gr) %>% unlist()

#  Go with this.
genome_CAST_BS %>% getSeq(., names = gr) %>% unlist()
genome_CAST_Nmasked_BS %>% getSeq(., names = gr_liftOver_CAST_to_mm10) %>% unlist()
genome_CAST_inserted_BS %>% getSeq(., names = gr_liftOver_CAST_to_mm10) %>% unlist()
genome_CAST_129_Nmasked_BS %>% getSeq(., names = gr_liftOver_CAST_to_mm10) %>% unlist()
genome_CAST_129_inserted_BS %>% getSeq(., names = gr_liftOver_CAST_to_mm10) %>% unlist()


# -----------------------------------------------------------------------------
string_1 <- genome_CAST_BS %>%
    BSgenome::getSeq(., names = gr) %>%
    as.character()
string_2 <- genome_CAST_Nmasked_BS %>%
    BSgenome::getSeq(., names = gr_liftOver_CAST_to_mm10) %>%
    as.character()


reticulate::use_condaenv(
    "/Users/kalavattam/miniconda3/envs/r41_env",
    required = TRUE
)


python_code <- 
"
def compare(string_1, string_2, no_match_c=' ', match_c='|'):
    if len(string_2) < len(string_1):
        string_1, string_2 = string_2, string_1
    result = ''
    n_diff = 0
    for c1, c2 in zip(string_1, string_2):
        if c1 == c2:
            result += match_c
        else:
            result += no_match_c
            n_diff += 1
    delta = len(string_2) - len(string_1)
    result += delta * no_match_c
    n_diff += delta
    return (result, n_diff)


def main():
    string_1 = r.string_1
    string_2 = r.string_2
    result, n_diff = compare(string_1, string_2, no_match_c='*')
    print(string_1)
    print(result)
    print(string_2)
    print('''%d difference(s).\n''' % n_diff)
    

main()
"


python_code %>% reticulate::py_run_string()
