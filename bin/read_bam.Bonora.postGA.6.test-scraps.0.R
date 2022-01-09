c_CAST129Nmasked <- sapply(
    row, compSeqCAST129Nmasked, bam = bam, genome = MmusculusCAST129Nmasked
)
c_CAST129Inserted <- sapply(
    row, compSeqCAST129Inserted, bam = bam, genome = MmusculusCAST129Inserted
)
c_129 <- sapply(
    row, compSeq129S1, bam = bam
)


paste0(
    "\n\n# ----------------------------------------------------\n",
    "mm10 assembly, N-masked at sites of 129S1-SvImJ and CAST-EiJ SNPs;",
    "\ncoordinates selected by GB's approach\n",
    "# ----------------------\n",
    c_CAST129Nmasked,
    "\n\nmm10 assembly, 129S1-SvImJ and CAST-EiJ SNPs inserted;",
    "\ncoordinates selected by GB's approach\n",
    "# ----------------------\n",
    c_CAST129Inserted,
    "\n\n129S1-SvImJ assembly;",
    "\ncoordinates selected by KA's approach\n",
    "# ----------------------\n",
    c_129,
    "\n\n"
    ) %>% 
    cat()


#TODO  Functions in the process of being rewritten, debugged...
compSeqCAST129Inserted <- function(bam, row, genome = MmusculusCAST129Inserted) {
    chr <- bam$rname[row] %>% as.character()  # rname = liftOver_rname (new_name = old_name), except not in the case of KA.ambiguous
    start <- bam$pos[row]  # pos = liftOver_pos (new_name = old_name), except not in the case of KA.ambiguous
    end <- bam$pos_end[row]  # pos_end = liftOver_pos_end (new_name = old_name), except not in the case of KA.ambiguous
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq[row]  #TODO Fix this...
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this...
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    comparison <- paste0(
        p$string_1, "\n",
        p$result, "\n",
        p$string_2, "\n",
        p$n_diff, " difference(s)."
    )
    return(comparison)
}


compSeqCAST129Nmasked <- function(bam, row, genome = MmusculusCAST129Nmasked) {
    chr <- bam$rname[row] %>% as.character()  # rname = liftOver_rname (new_name = old_name), except not in the case of KA.ambiguous
    start <- bam$pos[row]  # pos = liftOver_pos (new_name = old_name), except not in the case of KA.ambiguous
    end <- bam$pos_end[row]  # pos_end = liftOver_pos_end (new_name = old_name), except not in the case of KA.ambiguous
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq[row]  #TODO Fix this...
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this...
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    comparison <- paste0(
        p$string_1, "\n",
        p$result, "\n",
        p$string_2, "\n",
        p$n_diff, " difference(s)."
    )
    return(comparison)
}


compSeq129S1 <- function(bam, row, genome = Mmusculus129S1SvImJ) {
    chr <- bam$old_rname[row] %>% as.character()  # old_rname = rname (new_name = old_name), except not in the case of KA.ambiguous
    start <- bam$old_pos[row]  # old_pos = pos (new_name = old_name), except not in the case of KA.ambiguous
    end <- bam$old_pos_end[row]  # old_end = end (new_name = old_name), except not in the case of KA.ambiguous
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- bam$seq[row]  #TODO Fix this...
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()  #TODO Fix this...
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    comparison <- paste0(
        p$string_1, "\n",
        p$result, "\n",
        p$string_2, "\n",
        p$n_diff, " difference(s)."
    )
    return(comparison)
}

