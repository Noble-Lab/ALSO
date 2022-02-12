
library(GenomicAlignments)
library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)

#  Change to appropriate working directory
project <- "2021_kga0_4dn-mouse-cross"
default <- "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
current <- stringr::str_split(getwd(), "/")[[1]][
    length(stringr::str_split(getwd(), "/")[[1]])
]
if(current != project) {
    if(dir.exists(default)) {
        setwd(default)
    } else {
        setwd(
            readline(
                prompt = paste0(
                    "Enter path to and including the project directory, ",
                    project, ":"
                )
            )
        )
    }
}
rm(project, default, current)


#  Set up work directory, etc. (locations TB∆) --------------------------------
directory_user <- "/Users/kalavattam"
directory_project <- "Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
directory_bam <- "results/kga0/2022-0121_segregated_reads.thresh_SNP_1.thresh_Q_30"
directory_data <- "data/kga0"
path.1 <- paste0(directory_user, "/", directory_project)
path.2 <- paste0(directory_bam)
path.3 <- paste0(directory_data)

rm(directory_bam, directory_data, directory_project, directory_user)

setwd(paste0(path.1, "/", path.3))


#  Script name ----------------------------------------------------------------
script <- "test.PE-processing.part-8.R"


#  Load in datasets -----------------------------------------------------------

#  Having run part 7, load *.Rdata file
load(file = "test.PE-processing.part-7.Rdata")

#  Load in .Rdata from 'test.liftedOver-PE-processing.part-1.R'
load(file = "test.liftedOver-PE-processing.part-1.Rdata")


#  Perform analyses of SNPs, etc. for [[1]] KA.129 × GB.Ambiguous -------------
bam <- categories[[1]]

bam <- bam %>%
    dplyr::select(-c(cigar.GB.odd.1, cigar.GB.even.1)) %>%
    dplyr::rename(cigar.GB.odd = cigar.GB.odd.2) %>%
    dplyr::rename(cigar.GB.even = cigar.GB.even.2)

# (bam$matches.KA.odd == bam$matches.GB.odd) %>% table(useNA = "ifany")
# (bam$matches.KA.even == bam$matches.GB.even) %>% table(useNA = "ifany")
# names(bam)

#  Randomly sample 'int' number of rows from tibble 'bam'
# int <- 5L
# int <- 100L
# int <- 500L
# int <- 1000L
# bam <- bam[
#     sample(
#         1:nrow(bam),
#         if (nrow(bam) < int) {
#             nrow(bam)
#         } else {
#             int
#         },
#         replace = FALSE
#     ), 
# ]
bam <- bam %>%
    tibble::rowid_to_column() %>%
    dplyr::rename(ID.post_sample = rowid)

row <- bam$ID.post_sample %>% as.double()
rm(int)


# -----------------------------------------------------------------------------
cigar.KA.odd <- cigarOpTable(bam$cigar.KA.odd) %>% as_tibble()
cigar.KA.even <- cigarOpTable(bam$cigar.KA.even) %>% as_tibble()
cigar.GB.odd <- cigarOpTable(bam$cigar.GB.odd) %>% as_tibble()
cigar.GB.even <- cigarOpTable(bam$cigar.GB.even) %>% as_tibble()
cigar.KA.both <- cigar.KA.odd + cigar.KA.odd
cigar.GB.both <- cigar.GB.odd + cigar.GB.odd

variable <- c(
    "cigar.KA.odd",
    "cigar.KA.even",
    "cigar.KA.both",
    "cigar.GB.odd",
    "cigar.GB.even",
    "cigar.GB.both"
)
operation <- paste0(
    "colnames(", variable, ")", " <- ",
        "paste0(",
            "colnames(", variable, "), \"",
            stringr::str_remove(variable, "cigar"), "\"",
        ")"
)
evaluateOperation(operation)

bam <- dplyr::bind_cols(
    bam,
    cigar.KA.odd,
    cigar.KA.even,
    cigar.GB.odd,
    cigar.GB.even,
    cigar.KA.both,
    cigar.GB.both
)


# -----------------------------------------------------------------------------
names(bam)
category <- bam %>% dplyr::select(
    ID.post_sample,
    read_sequence.KA.odd, read_sequence.KA.even,
    matches.KA.odd, matches.KA.even,
    reference_sequence.KA.odd, reference_sequence.KA.even,
    read_sequence.GB.odd, read_sequence.GB.even,
    matches.GB.odd, matches.GB.even,
    reference_sequence.GB.odd, reference_sequence.GB.even,
    cigar.KA.odd, cigar.KA.even,
    cigar.GB.odd, cigar.GB.even,
    rname.129S1.odd, pos.129S1.odd, pos_end.129S1.odd,
    rname.129S1.even, pos.129S1.even, pos_end.129S1.even,
    rname.GB.odd, pos.GB.odd, pos_end.GB.odd,
    rname.GB.even, pos.GB.even, pos_end.GB.even,
    M.KA.both, M.GB.both, I.KA.both, I.GB.both, D.KA.both, D.GB.both,
    N.KA.both, N.GB.both, S.KA.both, S.GB.both, H.KA.both, H.GB.both,
    P.KA.both, P.GB.both, `=.KA.both`, `=.GB.both`, X.KA.both, X.GB.both
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
category$cigar.GB.both <- paste0(
    category$cigar.GB.odd, " × ", category$cigar.GB.even
)
category$read_sequence.GB.both <- paste0(
    category$read_sequence.GB.odd, " × ", category$read_sequence.GB.even
)
category$matches.GB.both <- paste0(
    category$matches.GB.odd, " × ", category$matches.GB.even
)
category$reference_sequence.GB.both <- paste0(
    category$reference_sequence.GB.odd, " × ", category$reference_sequence.GB.even
)

category <- category %>% dplyr::relocate(
    c(
        M.KA.both, M.GB.both, I.KA.both, I.GB.both, D.KA.both, D.GB.both,
        N.KA.both, N.GB.both, S.KA.both, S.GB.both, H.KA.both, H.GB.both,
        P.KA.both, P.GB.both, `=.KA.both`, `=.GB.both`, X.KA.both, X.GB.both
    ),
    .after = reference_sequence.GB.both
)


# -----------------------------------------------------------------------------
category$no_differences.KA <- stringi::stri_count(
    category$matches.KA.both, fixed = "*"
)
category$no_differences.GB <- stringi::stri_count(
    category$matches.GB.both, fixed = "*"
)


#  Looking at all entries in the cell -----------------------------------------
# category <- category.init
out <- paste0(
    "\n",
    "KA", "\n",
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
    " general difference(s)", "\n",
    "GB", "\n",
    "read pair:", " ",
    category$rname.GB.odd, ":",
    category$pos.GB.odd, "-",
    category$pos_end.GB.odd, " × ",
    category$rname.GB.even, ":",
    category$pos.GB.even, "-",
    category$pos_end.GB.even, "\n",
    "    CIGAR:", " ",
    category$cigar.GB.both, "\n",
    "     read:", " ",
    category$read_sequence.GB.both, "\n",
    "  matches:", " ",
    stringr::str_replace(category$matches.GB.both, "×", " "), "\n",
    "reference:", " ",
    category$reference_sequence.GB.both, "\n",
    "          ", " ",
    category$no_differences.GB,
    " general difference(s)", "\n\n"
)
out %>% cat()
# out %>% head() %>% cat()

#  Histograms for '*.both'
hist(
    category$no_differences.KA,
    breaks = max(category$no_differences.KA),
    main = paste0("Alignments assigned '129' by KA pipeline\n(best AS score from alignment to\n129 assembly, not CAST assembly)"),
    xlab = "General differences",
    right = FALSE
)
hist(
    category$no_differences.GB,
    breaks = max(category$no_differences.GB),
    main = paste0("Alignments assigned 'Ambiguous' by GB pipeline\n(reads do not contain '129' or 'CAST' SNPs)"),
    xlab = "General differences",
    right = FALSE
)

hist(
    category$`=.KA.both`,
    # breaks = max(category$no_differences.GB),
    main = paste0("Alignments assigned '129' by KA pipeline:\nSequence matches"),
    xlab = "No. sequence macthes",
    right = FALSE
)
hist(
    category$`=.GB.both`,
    # breaks = max(category$no_differences.GB),
    main = paste0("Alignments assigned 'Ambiguous' by GB pipeline:\nSequence matches"),
    xlab = "No. sequence macthes",
    right = FALSE
)

hist(
    category$X.KA.both,
    breaks = max(category$no_differences.GB),
    main = paste0("Alignments assigned '129' by KA pipeline:\nSequence mismatches"),
    xlab = "No. sequence mismacthes",
    right = FALSE
)
hist(
    category$X.GB.both,
    breaks = max(category$no_differences.GB),
    main = paste0("Alignments assigned 'Ambiguous' by GB pipeline:\nSequence mismatches"),
    xlab = "No. sequence mismacthes",
    right = FALSE
)

hist(
    category$I.KA.both,
    breaks = max(category$I.KA.both),
    main = paste0("Alignments assigned '129' by KA pipeline:\nInsertions to the reference"),
    xlab = "No. insertions to the reference",
    right = FALSE
)
hist(
    category$I.GB.both,
    breaks = max(category$I.GB.both),
    main = paste0("Alignments assigned 'Ambiguous' by GB pipeline:\nInsertions to the reference"),
    xlab = "No. insertions to the reference",
    right = FALSE
)

hist(
    category$D.KA.both,
    breaks = max(category$D.KA.both),
    main = paste0("Alignments assigned '129' by KA pipeline:\nDeletions from the reference"),
    xlab = "No. deletions from the reference",
    right = FALSE
)
hist(
    category$D.GB.both,
    breaks = max(category$D.GB.both),
    main = paste0("Alignments assigned 'Ambiguous' by GB pipeline:\nDeletions from the reference"),
    xlab = "No. deletions from the reference",
    right = FALSE
)

hist(
    category$M.KA.both,
    breaks = max(category$M.KA.both),
    main = paste0("Alignments assigned '129' by KA pipeline:\nAlignments to N in reference"),
    xlab = "No. of alignments to N in reference (from soft-clipping)",
    right = FALSE
)
hist(
    category$M.GB.both,
    breaks = max(category$M.GB.both),
    main = paste0("Alignments assigned 'Ambiguous' by GB pipeline:\nAlignments to N in reference"),
    xlab = "No. of alignments to N in reference (from N-insertion or soft-clipping)",
    right = FALSE
)

write.table(
    out,
    "2022-0207.KA.129-by-GB.Ambiguous.sample-all.txt",
    row.names = FALSE,
    col.names = FALSE
)

category.init <- category


#  Looking at all entries with M > 0 ------------------------------------------
category <- category.init[category.init$M.GB.both > 0 | category.init$M.KA.both > 0, ]
# category <- category.init[category.init$M.KA.both > 0, ]
out <- paste0(
    "\n",
    "KA", "\n",
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
    " general difference(s)", "\n",
    "GB", "\n",
    "read pair:", " ",
    category$rname.GB.odd, ":",
    category$pos.GB.odd, "-",
    category$pos_end.GB.odd, " × ",
    category$rname.GB.even, ":",
    category$pos.GB.even, "-",
    category$pos_end.GB.even, "\n",
    "    CIGAR:", " ",
    category$cigar.GB.both, "\n",
    "     read:", " ",
    category$read_sequence.GB.both, "\n",
    "  matches:", " ",
    stringr::str_replace(category$matches.GB.both, "×", " "), "\n",
    "reference:", " ",
    category$reference_sequence.GB.both, "\n",
    "          ", " ",
    category$no_differences.GB,
    " general difference(s)", "\n\n"
)
out %>% cat()
# out %>% head() %>% cat()

write.table(
    out,
    "2022-0207.KA.129-by-GB.Ambiguous.sample-all.M-gt-0.txt",
    row.names = FALSE,
    col.names = FALSE
)
# rm(out)

#  I dont know what i was expecting; i guess i was expecting to see more differences
#+ between KA.129 and GB.Ambiguous, esp in terms of I and D, but these are essentially the same;
#+ instead, the biggest differences are seen in terms of M (Ns from s-c and N-insertion) and
#+ "general differences" (again from N-insertion)

#+ So, that leads to the question, why 129 from KA and Ambiguous from GB? I think's
#+ b/c my assignments come from only comparing CAST AS to 129 AS; I should also compare
#+ to mm10 or N-masked and if, for example, 129 AS is gt CAST but the same as mm10-N-masked or mm10,
#+ should that really be assigned 129?
