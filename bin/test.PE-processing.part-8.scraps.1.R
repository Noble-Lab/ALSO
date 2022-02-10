
#  Perform analyses of SNPs, etc. for [[5]] KA.129 × GB.NA --------------------
bam <- categories[[5]]

#  Add in matches to *.odd
bind <- bind.129S1 %>%
    dplyr::rename(coordinate.KA.129S1.odd = coordinate.standard)
colnames(bind) <- colnames(bind) %>%
    gsub("coordinate.standard", "coordinate.KA.129S1.odd", .) %>%
    gsub("standard", "standard.odd", .) %>%
    gsub("sam2pairwise", "sam2pairwise.odd", .)

bam <- bam %>% dplyr::left_join(
    bind[!duplicated(bind$coordinate.KA.129S1.odd), ]
)
rm(bind)

#  Add in matches to *.even
bind <- bind.129S1 %>%
    dplyr::rename(coordinate.KA.129S1.even = coordinate.standard)
colnames(bind) <- colnames(bind) %>%
    gsub("coordinate.standard", "coordinate.KA.129S1.even", .) %>%
    gsub("standard", "standard.even", .) %>%
    gsub("sam2pairwise", "sam2pairwise.even", .)

bam <- bam %>% dplyr::left_join(
    bind[!duplicated(bind$coordinate.KA.129S1.even), ]
)
rm(bind)

bam$matches.GB.odd %>% head()
bam$matches.KA.odd %>% head(n = 20)
bam$matches.sam2pairwise.odd %>% head(n = 20)
bam$cigar.KA.odd %>% head(20)
bam$cigar.sam2pairwise.odd %>% head(20)


bam$matches.GB.even %>% head()
bam$matches.KA.even %>% head(n = 20)
bam$matches.sam2pairwise.even %>% head(n = 20)
bam$cigar.KA.even %>% head(20)
bam$cigar.sam2pairwise.even %>% head(20)

(bam$matches.KA.odd == bam$matches.sam2pairwise.odd) %>% table()
(bam$matches.KA.even == bam$matches.sam2pairwise.even) %>% table()

(bam$cigar.KA.odd == bam$cigar.sam2pairwise.odd) %>% table()
(bam$cigar.KA.even == bam$cigar.sam2pairwise.even) %>% table()

#NOTE This does not work: These never aligned to the Nmasked mm10 genome
#TODO 1/3 Need to associate the initial coordinate information with the lifted-
#TODO 2/3 over information, then use the initial information to join the
#TODO 3/3 lifted-over information to tibble "bam"

#  Randomly sample 'int' number of rows from tibble 'bam'; if commented out,
#+ then process all observations in the tibble
# int <- 5L
# int <- 100L
# int <- 500L
int <- 1000L
bam <- bam[
    sample(
        1:nrow(bam),
        if (nrow(bam) < int) {
            nrow(bam)
        } else {
            int
        },
        replace = FALSE
    ), 
]
bam <- bam %>%
    tibble::rowid_to_column() %>%
    dplyr::rename(ID.post_sample = rowid)

row <- bam$ID.post_sample %>% as.double()
rm(int)


# -----------------------------------------------------------------------------
category <- bam %>% dplyr::select(
    ID.post_sample,
    read_sequence.KA.odd, read_sequence.KA.even,
    matches.KA.odd, matches.KA.even,
    reference_sequence.KA.odd, reference_sequence.KA.even,
    # read_sequence.GB.odd, read_sequence.GB.even,
    # matches.GB.odd, matches.GB.even,
    # reference_sequence.GB.odd, reference_sequence.GB.even,
    read_sequence.sam2pairwise.odd, read_sequence.sam2pairwise.even,
    matches.sam2pairwise.odd, matches.sam2pairwise.even,
    reference_sequence.sam2pairwise.odd, reference_sequence.sam2pairwise.even,
    cigar.KA.odd, cigar.KA.even,
    # cigar.GB.odd.2, cigar.GB.even.2,
    cigar.sam2pairwise.odd, cigar.sam2pairwise.even,
    rname.129S1.odd, pos.129S1.odd, pos_end.129S1.odd,
    rname.129S1.even, pos.129S1.even, pos_end.129S1.even,
    # rname.GB.odd, pos.GB.odd, pos_end.GB.odd,
    # rname.GB.even, pos.GB.even, pos_end.GB.even,
    rname.standard.odd, pos.standard.odd,
    rname.standard.even, pos.standard.even 
)

category$cigar.KA.both <- paste0(
    category$cigar.KA.odd, " × ", category$cigar.KA.even
)
category$read_sequence.KA.both <- paste0(
    category$read_sequence.KA.odd, " × ", category$read_sequence.KA.even
)
category$matches.KA.both <- paste0(
    category$matches.KA.odd, " × ", category$matches.KA.even
)
category$reference_sequence.KA.both <- paste0(
    category$reference_sequence.KA.odd, " × ", category$reference_sequence.KA.even
)
# category$cigar.GB.both <- paste0(
#     category$cigar.GB.odd.2, " × ", category$cigar.GB.even.2
# )
# category$read_sequence.GB.both <- paste0(
#     category$read_sequence.GB.odd, " × ", category$read_sequence.GB.even
# )
# category$matches.GB.both <- paste0(
#     category$matches.GB.odd, " × ", category$matches.GB.even
# )
# category$reference_sequence.GB.both <- paste0(
#     category$reference_sequence.GB.odd, " × ", category$reference_sequence.GB.even
# )
category$cigar.sam2pairwise.both <- paste0(
    category$cigar.sam2pairwise.odd, " × ", category$cigar.sam2pairwise.even
)
category$read_sequence.sam2pairwise.both <- paste0(
    category$read_sequence.sam2pairwise.odd, " × ", category$read_sequence.sam2pairwise.even
)
category$matches.sam2pairwise.both <- paste0(
    category$matches.sam2pairwise.odd, " × ", category$matches.sam2pairwise.even
)
category$reference_sequence.sam2pairwise.both <- paste0(
    category$reference_sequence.sam2pairwise.odd, " × ", category$reference_sequence.sam2pairwise.even
)


# -----------------------------------------------------------------------------
category$no_differences.KA <- stringi::stri_count(
    category$matches.KA.both, fixed = "*"
)
category$no_differences.sam2pairwise <- stringi::stri_count(
    category$matches.sam2pairwise.both, fixed = "*"
)

category$cigar.KA.both %>% table(., useNA = "ifany")
category$cigar.sam2pairwise.both %>% table(., useNA = "ifany")


# -----------------------------------------------------------------------------
out.KA <- paste0(
    "\n",
    "read pair:", " ",
    category$rname.129S1.odd, ":",
    category$pos.129S1.odd, "-",
    category$pos_end.129S1.odd, " × ",
    category$rname.129S1.even, ":",
    category$pos.129S1.even, "-",
    category$pos_end.129S1.even, "\n",
    "    CIGAR:", " ",
    category$cigar.KA.both, "\n",
    "     read:", " ",
    category$read_sequence.KA.both, "\n",
    "  matches:", " ",
    stringr::str_replace(category$matches.KA.both, "×", " "), "\n",
    "reference:", " ",
    category$reference_sequence.KA.both, "\n",
    "          ", " ",
    category$no_differences.KA,
    " general difference(s)", "\n"
)
out.KA %>% cat()

out.CO <- paste0(
    "\n",
    "read pair:", " ",
    category$rname.standard.odd, ":",
    category$pos.standard.odd, "-",
    (category$pos.standard.odd + 49), " × ",
    category$rname.standard.even, ":",
    category$pos.standard.even, "-",
    (category$pos.standard.even + 49), "\n",
    "    CIGAR:", " ",
    category$cigar.sam2pairwise.both, "\n",
    "     read:", " ",
    category$read_sequence.sam2pairwise.both, "\n",
    "  matches:", " ",
    stringr::str_replace(category$matches.sam2pairwise.both, "×", " "), "\n",
    "reference:", " ",
    category$reference_sequence.sam2pairwise.both, "\n",
    "          ", " ",
    category$no_differences.sam2pairwise,
    " general difference(s)", "\n"
)
out.CO %>% cat()

#  Histograms for '*.both'
hist(
    category$no_differences.KA,
    breaks = max(category$no_differences.KA),
    main = "129 alignments compared to 129 assembly",
    xlab = "General differences",
    right = FALSE
)
hist(
    category$no_differences.sam2pairwise,
    breaks = max(category$no_differences.sam2pairwise),
    main = "129 alignments lifted over to mm10 assembly",
    xlab = "General differences",
    right = FALSE
)

write.table(
    out,
    "2022-0201.sample-500.129-alignments-compared-to-129-assembly.txt",
    row.names = FALSE,
    col.names = FALSE
)

category_not50M <- category[category$cigar.KA.both != "50M × 50M", ]
out <- paste0(
    "\n",
    "read pair:", " ",
    category_not50M$rname.129S1.odd, ":",
    category_not50M$pos.129S1.odd, "-",
    category_not50M$pos_end.129S1.odd, " × ",
    category_not50M$rname.129S1.even, ":",
    category_not50M$pos.129S1.even, "-",
    category_not50M$pos_end.129S1.even, "\n",
    "    CIGAR:", " ",
    category_not50M$cigar.KA.both, "\n",
    "     read:", " ",
    category_not50M$read_sequence.KA.both, "\n",
    "  matches:", " ",
    stringr::str_replace(category_not50M$matches.KA.both, "×", " "), "\n",
    "reference:", " ",
    category_not50M$reference_sequence.KA.both, "\n",
    "          ", " ",
    category_not50M$no_differences,
    " mismatch(es)", "\n"
)
out %>% cat()

write.table(
    out,
    "2022-0201.sample-500.CIGAR-no-50M.129-alignments-compared-to-129-assembly.txt",
    row.names = FALSE,
    col.names = FALSE
)
# rm(out)



###############################################################################
# #  Set up liftOver chains -----------------------------------------------------
# liftover_directory <- "liftOver"
#
# #  129S1
# chain_129S1_to_mm10 <- paste0(
#     liftover_directory, "/", "129S1-SvImJ-to-mm10.over.chain.munged"
# ) %>%
#     rtracklayer::import.chain()
# chain_mm10_to_129S1 <- paste0(
#     liftover_directory, "/", "mm10-to-129S1-SvImJ.over.chain.munged"
# ) %>%
#     rtracklayer::import.chain()
#
# #  CAST
# chain_CAST_to_mm10 <- paste0(
#     liftover_directory, "/", "CAST-EiJ-to-mm10.over.chain.munged"
# ) %>%
#     rtracklayer::import.chain()
# chain_mm10_to_CAST <- paste0(
#     liftover_directory, "/", "mm10-to-CAST-EiJ.over.chain.munged"
# ) %>%
#     rtracklayer::import.chain()


#  Set up Python code ---------------------------------------------------------
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


string_1 = r.string_1
string_2 = r.string_2
result, n_diff = compare(string_1, string_2, no_match_c='*')
"


#  Generate tibble of comparisons, metrics, etc. for [[5]] KA.129 × GB.NA -----
sample.1000.Mmusculus129S1SvImJ <- joinAndSelectWorst(
    row = row,
    genome = Mmusculus129S1SvImJ,
    rname.odd = bam$rname.129S1.odd,
    pos.odd = bam$pos.129S1.odd,
    pos_end.odd = bam$pos_end.129S1.odd,
    seq.odd = bam$seq.129S1.odd,
    rname.even = bam$rname.129S1.even,
    pos.even = bam$pos.129S1.even,
    pos_end.even = bam$pos_end.129S1.even,
    seq.even = bam$seq.129S1.even
)

sample.1000.Mmusculus129Nmasked <- joinAndSelectWorst(
    row = row,
    genome = Mmusculus129Nmasked,
    rname.odd = bam$lO_rname.129S1.odd,
    pos.odd = bam$lO_pos.129S1.odd,
    pos_end.odd = bam$lO_pos.129S1.odd + 49,
    seq.odd = bam$seq.129S1.odd,
    rname.even = bam$lO_rname.129S1.even,
    pos.even = bam$lO_pos.129S1.even,
    pos_end.even = bam$lO_pos.129S1.even + 49,
    seq.even = bam$seq.129S1.even
    
    # pos_end.odd = bam$lO_pos_end.129S1.odd,
    # pos_end.even = bam$lO_pos_end.129S1.even,
)

#  Histograms for '*.both'
hist(
    sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both,
    # breaks = max(sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both),
    main = "129 alignments compared to 129 assembly",
    xlab = "Mismatches",
    right = FALSE
)

hist(
    sample.1000.Mmusculus129Nmasked$no_difference.sum.both,
    # breaks = max(sample.1000.Mmusculus129Nmasked$no_difference.sum.both),
    main = "129 alignments lifted to mm10 coordinates; then,\nlifted alignments compared to mm10 assembly\nN-masked at 129 SNPs",
    xlab = "Mismatches",
    right = FALSE
)

hist(
    sample.1000.Mmusculus129Nmasked$no_diff_not_N.both,
    # breaks = max(sample.1000.Mmusculus129Nmasked$no_difference.sum.both),
    main = "129 alignments lifted to mm10 coordinates; then,\nlifted alignments compared to mm10 assembly\nN-masked at 129 SNPs",
    xlab = "Mismatches (not including SNPs)",
    right = FALSE
)

hist(
    sample.1000.Mmusculus129Nmasked$no_N.both,
    breaks = max(sample.1000.Mmusculus129Nmasked$no_N.both),
    main = "129 alignments lifted to mm10 coordinates; then,\nlifted alignments compared to mm10 assembly\nN-masked at 129 SNPs",
    xlab = "Numbers of N-masked SNPs",
    right = FALSE
)


out <- paste0(
    "\n#######################################################################################################\n",
    "129 alignments compared to the 129 assembly", "\n",
    "-------------------------------------------------------------------------------------------------------", "\n",
    sample.1000.Mmusculus129S1SvImJ$seq_query.both, "\n",
    sample.1000.Mmusculus129S1SvImJ$seq_match.both, "\n",
    sample.1000.Mmusculus129S1SvImJ$seq_genome.both, "\n",
    sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both, " total difference(s)", "\n",
    sample.1000.Mmusculus129S1SvImJ$no_difference.both, " total difference(s)", "\n",
    sample.1000.Mmusculus129S1SvImJ$no_diff_not_N.both, " difference(s) that are not N-related", "\n",
    sample.1000.Mmusculus129S1SvImJ$no_N.both, " difference(s) that *are* N-related", "\n",
    "\n",
    "129 alignments lifted to mm10 coordinates; then, lifted alignments compared to mm10 assembly N-masked at 129 SNPs", "\n",
    "-------------------------------------------------------------------------------------------------------", "\n",
    sample.1000.Mmusculus129Nmasked$seq_query.both, "\n",
    sample.1000.Mmusculus129Nmasked$seq_match.both, "\n",
    sample.1000.Mmusculus129Nmasked$seq_genome.both, "\n",
    sample.1000.Mmusculus129Nmasked$no_difference.sum.both, " total difference(s)", "\n",
    sample.1000.Mmusculus129Nmasked$no_difference.both, " total difference(s)", "\n",
    sample.1000.Mmusculus129Nmasked$no_diff_not_N.both, " difference(s) that are not N-related", "\n",
    sample.1000.Mmusculus129Nmasked$no_N.both, " difference(s) that *are* N-related", "\n",
    "\n"
)
write.table(
    out,
    "2022-0124.sample.1000.sample.1000.Mmusculus129S1SvImJ.sample.1000.Mmusculus129Nmasked.txt",
    row.names = FALSE,
    col.names = FALSE
)
rm(out)

#  sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10
out <- paste0(
    "\n#######################################################################################################\n",
    "129 alignments compared to the 129 assembly", "\n",
    "-------------------------------------------------------------------------------------------------------", "\n",
    sample.1000.Mmusculus129S1SvImJ$seq_query.both[sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10], "\n",
    sample.1000.Mmusculus129S1SvImJ$seq_match.both[sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10], "\n",
    sample.1000.Mmusculus129S1SvImJ$seq_genome.both[sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10], "\n",
    sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both[sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10], " total difference(s)", "\n",
    sample.1000.Mmusculus129S1SvImJ$no_difference.both[sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10], " total difference(s)", "\n",
    sample.1000.Mmusculus129S1SvImJ$no_diff_not_N.both[sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10], " difference(s) that are not N-related", "\n",
    sample.1000.Mmusculus129S1SvImJ$no_N.both[sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10], " difference(s) that *are* N-related", "\n",
    "\n",
    "129 alignments lifted to mm10 coordinates; then, lifted alignments compared to mm10 assembly N-masked at 129 SNPs", "\n",
    "-------------------------------------------------------------------------------------------------------", "\n",
    sample.1000.Mmusculus129Nmasked$seq_query.both[sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10], "\n",
    sample.1000.Mmusculus129Nmasked$seq_match.both[sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10], "\n",
    sample.1000.Mmusculus129Nmasked$seq_genome.both[sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10], "\n",
    sample.1000.Mmusculus129Nmasked$no_difference.sum.both[sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10], " total difference(s)", "\n",
    sample.1000.Mmusculus129Nmasked$no_difference.both[sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10], " total difference(s)", "\n",
    sample.1000.Mmusculus129Nmasked$no_diff_not_N.both[sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10], " difference(s) that are not N-related", "\n",
    sample.1000.Mmusculus129Nmasked$no_N.both[sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10], " difference(s) that *are* N-related", "\n",
    "\n"
)
write.table(
    out,
    "2022-0124.sample.1000.sample.1000.Mmusculus129S1SvImJ.sample.1000.Mmusculus129Nmasked.Mmusculus129Nmasked-gt-10-differences.txt",
    row.names = FALSE,
    col.names = FALSE
)
rm(out)


#  sample.1000.Mmusculus129Nmasked$no_difference.sum.both > 10
out <- paste0(
    "\n#######################################################################################################\n",
    "129 alignments compared to the 129 assembly", "\n",
    "-------------------------------------------------------------------------------------------------------", "\n",
    sample.1000.Mmusculus129S1SvImJ$seq_query.both[sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both > 10], "\n",
    sample.1000.Mmusculus129S1SvImJ$seq_match.both[sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both > 10], "\n",
    sample.1000.Mmusculus129S1SvImJ$seq_genome.both[sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both > 10], "\n",
    sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both[sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both > 10], " total difference(s)", "\n",
    sample.1000.Mmusculus129S1SvImJ$no_difference.both[sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both > 10], " total difference(s)", "\n",
    sample.1000.Mmusculus129S1SvImJ$no_diff_not_N.both[sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both > 10], " difference(s) that are not N-related", "\n",
    sample.1000.Mmusculus129S1SvImJ$no_N.both[sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both > 10], " difference(s) that *are* N-related", "\n",
    "\n",
    "129 alignments lifted to mm10 coordinates; then, lifted alignments compared to mm10 assembly N-masked at 129 SNPs", "\n",
    "-------------------------------------------------------------------------------------------------------", "\n",
    sample.1000.Mmusculus129Nmasked$seq_query.both[sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both > 10], "\n",
    sample.1000.Mmusculus129Nmasked$seq_match.both[sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both > 10], "\n",
    sample.1000.Mmusculus129Nmasked$seq_genome.both[sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both > 10], "\n",
    sample.1000.Mmusculus129Nmasked$no_difference.sum.both[sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both > 10], " total difference(s)", "\n",
    sample.1000.Mmusculus129Nmasked$no_difference.both[sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both > 10], " total difference(s)", "\n",
    sample.1000.Mmusculus129Nmasked$no_diff_not_N.both[sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both > 10], " difference(s) that are not N-related", "\n",
    sample.1000.Mmusculus129Nmasked$no_N.both[sample.1000.Mmusculus129S1SvImJ$no_difference.sum.both > 10], " difference(s) that *are* N-related", "\n",
    "\n"
)
write.table(
    out,
    "2022-0124.sample.1000.sample.1000.Mmusculus129S1SvImJ.sample.1000.Mmusculus129Nmasked.Mmusculus129S1SvImJ-gt-10-differences.txt",
    row.names = FALSE,
    col.names = FALSE
)
rm(out)


#  Cell-wise analyses of SNPs, etc. for [[17]] KA.NA × GB.Ambiguous -----------
bam <- joint.assign[[17]] %>% dplyr::select(sort(dplyr::current_vars()))

#  Randomly sample 'int' number of rows from tibble 'bam'; if commented out,
#+ then process all observations in the tibble
# int <- 5L
# int <- 100L
# int <- 500L
int <- 1000L
bam <- bam[
    sample(
        1:nrow(bam),
        if (nrow(bam) < int) {
            nrow(bam)
        } else {
            int
        },
        replace = FALSE
    ), 
]
bam <- bam %>%
    tibble::rowid_to_column() %>%
    dplyr::rename(ID.post_sample = rowid)

row <- bam$ID.post_sample %>% as.double()
rm(int)

# colnames(bam)


#  Generate tibble of comparisons, metrics, etc. for [[17]] KA.NA × GB.Ambiguous
sample.1000.MmusculusCAST129Nmasked <- joinAndSelectWorst(
    row = row,
    genome = MmusculusCAST129Nmasked,
    rname.odd = bam$rname.GB.odd,
    pos.odd = bam$pos.GB.odd,
    pos_end.odd = bam$pos_end.GB.odd,
    seq.odd = bam$seq.GB.odd,
    rname.even = bam$rname.GB.even,
    pos.even = bam$pos.GB.even,
    pos_end.even = bam$pos_end.GB.even,
    seq.even = bam$seq.GB.even
)

sample.1000.MmusculusCAST129Inserted <- joinAndSelectWorst(
    row = row,
    genome = MmusculusCAST129Inserted,
    rname.odd = bam$rname.GB.odd,
    pos.odd = bam$pos.GB.odd,
    pos_end.odd = bam$pos_end.GB.odd,
    seq.odd = bam$seq.GB.odd,
    rname.even = bam$rname.GB.even,
    pos.even = bam$pos.GB.even,
    pos_end.even = bam$pos_end.GB.even,
    seq.even = bam$seq.GB.even
)

#  Histograms for '*.both'
hist(
    sample.1000.MmusculusCAST129Nmasked$no_difference.sum.both,
    breaks = max(sample.1000.MmusculusCAST129Nmasked$no_difference.sum.both),
    main = "Alignments to mm10 reference N-masked at 129 and CAST SNPs\ncompared to mm10 assembly N-masked at 129 and CAST SNPs",
    xlab = "Mismatches",
    right = FALSE
)

hist(
    sample.1000.MmusculusCAST129Inserted$no_difference.sum.both,
    breaks = max(sample.1000.MmusculusCAST129Inserted$no_difference.sum.both),
    main = "Alignments to mm10 reference N-masked at 129 and CAST SNPs\ncompared to mm10 assembly that includes 129 and CAST SNPs",
    xlab = "Mismatches",
    right = FALSE
)


out <- paste0(
    "\n#######################################################################################################\n",
    "Alignments to mm10 reference N-masked at 129 and CAST SNPs compared to mm10 assembly N-masked at 129 and CAST SNPs", "\n",
    "-------------------------------------------------------------------------------------------------------", "\n",
    sample.1000.MmusculusCAST129Nmasked$seq_query.both, "\n",
    sample.1000.MmusculusCAST129Nmasked$seq_match.both, "\n",
    sample.1000.MmusculusCAST129Nmasked$seq_genome.both, "\n",
    sample.1000.MmusculusCAST129Nmasked$no_difference.sum.both, " total difference(s)", "\n",
    sample.1000.MmusculusCAST129Nmasked$no_difference.both, " total difference(s)", "\n",
    sample.1000.MmusculusCAST129Nmasked$no_diff_not_N.both, " difference(s) that are not N-related", "\n",
    sample.1000.MmusculusCAST129Nmasked$no_N.both, " difference(s) that *are* N-related", "\n",
    "\n",
    "Alignments to mm10 reference that includes 129 and CAST SNPs compared to mm10 assembly that includes 129 and CAST SNPs", "\n",
    "-------------------------------------------------------------------------------------------------------", "\n",
    sample.1000.MmusculusCAST129Inserted$seq_query.both, "\n",
    sample.1000.MmusculusCAST129Inserted$seq_match.both, "\n",
    sample.1000.MmusculusCAST129Inserted$seq_genome.both, "\n",
    sample.1000.MmusculusCAST129Inserted$no_difference.sum.both, " total difference(s)", "\n",
    sample.1000.MmusculusCAST129Inserted$no_difference.both, " total difference(s)", "\n",
    sample.1000.MmusculusCAST129Inserted$no_diff_not_N.both, " difference(s) that are not N-related", "\n",
    sample.1000.MmusculusCAST129Inserted$no_N.both, " difference(s) that *are* N-related", "\n",
    "\n"
)

write.table(
    out,
    "2022-0124.sample.1000.MmusculusCAST129Nmasked.MmusculusCAST129Inserted.txt",
    row.names = FALSE,
    col.names = FALSE
)
rm(out)

#  Scraps, to-dos, etc. -------------------------------------------------------
#TODO 1/2 Get the 'seq' associated with the pmin(), i.e., the worse alignment;
#TODO 2/2 put them in a separate column; then use them for the above analyses

#TODO 1/3 Test seq.even and seq.odd for most mismatches; select the one with;
#TODO 2/3 the most mismatches; put it in its own column; then use them for the
#TODO 3/3 above analyses
