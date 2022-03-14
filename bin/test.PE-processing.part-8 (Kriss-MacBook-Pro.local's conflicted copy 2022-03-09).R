#!/usr/bin/env Rscript

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
directory_bam <- "results/kga0/2022-0214_segregated_reads.thresh_SNP_1.thresh_Q_30"
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

# #  Load in .Rdata from 'test.liftedOver-PE-processing.part-1.R'
# load(file = "test.liftedOver-PE-processing.part-1.Rdata")


#  Perform analyses of SNPs, etc. for [[1]] KA.129 × GB.Ambiguous -------------
bam <- categories[[1]]

bam <- bam %>%
    dplyr::select(-c(cigar.GB.odd.1, cigar.GB.even.1)) %>%
    dplyr::rename(cigar.GB.odd = cigar.GB.odd.2) %>%
    dplyr::rename(cigar.GB.even = cigar.GB.even.2)

bam <- bam %>%
    tibble::rowid_to_column() %>%
    dplyr::rename(ID.post_sample = rowid)

row <- bam$ID.post_sample %>% as.double()
rm(int)


#  Create cigarOpTable for each CIGAR entry in tibble 'bam' -------------------

#  cigarOpTable: Integer matrix with number of rows equal to the length of
#+ CIGAR and nine columns, one for each extended CIGAR operation
cigar.GB.odd <- GenomicAlignments::cigarOpTable(bam$cigar.GB.odd) %>%
    as_tibble()
cigar.GB.even <- GenomicAlignments::cigarOpTable(bam$cigar.GB.even) %>%
    as_tibble()
cigar.KA.129.odd <- GenomicAlignments::cigarOpTable(bam$cigar.KA.129.odd) %>%
    as_tibble()
cigar.KA.129.even <- GenomicAlignments::cigarOpTable(bam$cigar.KA.129.even) %>%
    as_tibble()
cigar.KA.CAST.odd <- GenomicAlignments::cigarOpTable(bam$cigar.KA.CAST.odd) %>%
    as_tibble()
cigar.KA.CAST.even <- GenomicAlignments::cigarOpTable(bam$cigar.KA.CAST.even) %>%
    as_tibble()

variable <- c(
    "cigar.GB.odd",
    "cigar.GB.even",
    "cigar.KA.129.odd",
    "cigar.KA.129.even",
    "cigar.KA.CAST.odd",
    "cigar.KA.CAST.even"
)

#  Make all row values NA if one value is NA
operation <- paste0(
    variable, "[!!rowSums(is.na(", variable, ")), ]", " <- ", "NA"
)
evaluateOperation(operation)

#  Sum the rows
command <- paste0("rowSums(", variable, ")")
operation <- makeOperation(paste0(variable, "$sum"), command)
evaluateOperation(operation)
#NOTE Any sum > 50 is due to presence of D, ∴ subtract D from =

#  Subtract D from =
command <- paste0(variable, "$`=`", " - ", variable, "$D")
operation <- makeOperation(paste0(variable, "$`=`"), command)
evaluateOperation(operation)

#  Sum the rows again
command <- paste0(variable, " %>% ", "dplyr::select(-sum)")
operation <- makeOperation(variable, command)
evaluateOperation(operation)

command <- paste0("rowSums(", variable, ")")
operation <- makeOperation(paste0(variable, "$sum"), command)
evaluateOperation(operation)

# colnames(cigar.GB.odd) <- paste0("GB.mm10.", colnames(cigar.GB.odd))
# colnames(cigar.KA.129.odd) <- paste0("KA.129.", colnames(cigar.KA.129.odd))
# colnames(cigar.KA.CAST.odd) <- paste0("KA.CAST.", colnames(cigar.KA.CAST.odd))

# #  Convert columns from integer to factor modes
# cigar.GB.odd[sapply(cigar.GB.odd, is.integer)] <-
#     lapply(
#         cigar.GB.odd[sapply(cigar.GB.odd, is.integer)],
#         as.factor
#     )
# 
# cigar.KA.129.odd[sapply(cigar.KA.129.odd, is.integer)] <-
#     lapply(
#         cigar.KA.129.odd[sapply(cigar.KA.129.odd, is.integer)],
#         as.factor
#     )

# -------------------------------------
gather.GB <- cigar.GB.odd %>%
    dplyr::select(-sum) %>%
    gather("category", "counts", 1:9)

list.GB <- vector(mode = "list", length = nrow(gather.GB))
for(i in 1:nrow(gather.GB)) {
    list.GB[[i]] <- if(is.na(gather.GB$counts[i])) {
        paste(NA_character_)
    } else if(gather.GB$counts[i] > 0) {
        paste(replicate(gather.GB$counts[i], gather.GB$category[i]))
    } else {
        paste("§")
    }
}

vector.GB <- unlist(list.GB)
remove <- "§"
vector.GB <- vector.GB[!vector.GB %in% remove] %>% tibble::as_tibble_col()
# vector.GB <- vector.GB %>% tibble::as_tibble_col()
vector.GB$value <- vector.GB$value %>%
    forcats::as_factor() %>%
    plyr::revalue(c("M" = "N"))

rm(remove, list.GB)

# -------------------------------------
gather.KA.129 <- cigar.KA.129.odd %>%
    dplyr::select(-sum) %>%
    gather("category", "counts", 1:9)

list.KA.129 <- vector(mode = "list", length = nrow(gather.KA.129))
for(i in 1:nrow(gather.KA.129)) {
    list.KA.129[[i]] <- if(is.na(gather.KA.129$counts[i])) {
        paste(NA_character_)
    } else if(gather.KA.129$counts[i] > 0) {
        paste(replicate(gather.KA.129$counts[i], gather.KA.129$category[i]))
    } else {
        paste("§")
    }
}

vector.KA.129 <- unlist(list.KA.129)
remove <- "§"
vector.KA.129 <- vector.KA.129[!vector.KA.129 %in% remove] %>% tibble::as_tibble_col()
# vector.KA.129 <- vector.KA.129 %>% tibble::as_tibble_col()
vector.KA.129$value <- vector.KA.129$value %>%
    forcats::as_factor() %>%
    plyr::revalue(c("M" = "N"))

rm(remove, list.KA.129)

# -------------------------------------
gather.KA.CAST <- cigar.KA.CAST.odd %>%
    dplyr::select(-sum) %>%
    gather("category", "counts", 1:9)

list.KA.CAST <- vector(mode = "list", length = nrow(gather.KA.CAST))
for(i in 1:nrow(gather.KA.CAST)) {
    list.KA.CAST[[i]] <- if(is.na(gather.KA.CAST$counts[i])) {
        paste(NA_character_)
    } else if(gather.KA.CAST$counts[i] > 0) {
        paste(replicate(gather.KA.CAST$counts[i], gather.KA.CAST$category[i]))
    } else {
        paste("§")
    }
}

vector.KA.CAST <- unlist(list.KA.CAST)
remove <- "§"
vector.KA.CAST <- vector.KA.CAST[!vector.KA.CAST %in% remove] %>% tibble::as_tibble_col()
# vector.KA.CAST <- vector.KA.CAST %>% tibble::as_tibble_col()
vector.KA.CAST$value <- vector.KA.CAST$value %>%
    forcats::as_factor() %>%
    plyr::revalue(c("M" = "N"))

rm(remove, list.KA.CAST)


#  Generate confusion matrix --------------------------------------------------

#  'data' is for predicted data, i.e., the data from my pipeline; 'reference'
#+ is for true results, i.e., the data from GB's pipeline
example <- caret::confusionMatrix(
    data = vector.GB$value,
    reference = vector.KA.129$value,
    dnn = c("GB", "KA")
)
example.table <- example$table %>%
    as.data.frame() %>%
    tibble::as_tibble()
example.table <- example.table %>%
    mutate(KA = factor(KA, rev(unique(KA))))

ggplot(example.table, aes(x = GB, y = KA, fill = Freq)) +
    geom_tile(color = "white") +
    theme_bw() +
    coord_equal() +
    scale_fill_distiller(palette = "Spectral", direction = -1) +
    guides(fill = FALSE) +
    labs(title = "Confusion matrix for N-masked\nmm10 (GB) vs. 129 (KA)") +
    geom_text(aes(label = Freq), color = "black")

rm(example, example.table)


#  Subset GB and KA.129 tibbles to non-NA rows in KA.CAST ---------------------
index <- which(!is.na(cigar.KA.CAST.odd$M), arr.ind = TRUE)
subset.GB.odd <- cigar.GB.odd[index, ]
subset.KA.129.odd <- cigar.KA.129.odd[index, ]
subset.KA.CAST.odd <- cigar.KA.CAST.odd[index, ]

# -------------------------------------
gather.GB <- subset.GB.odd %>%
    dplyr::select(-sum) %>%
    gather("category", "counts", 1:9)

list.GB <- vector(mode = "list", length = nrow(gather.GB))
for(i in 1:nrow(gather.GB)) {
    list.GB[[i]] <- if(is.na(gather.GB$counts[i])) {
        paste(NA_character_)
    } else if(gather.GB$counts[i] > 0) {
        paste(replicate(gather.GB$counts[i], gather.GB$category[i]))
    } else {
        paste("§")
    }
}

vector.GB <- unlist(list.GB)
remove <- "§"
vector.GB <- vector.GB[!vector.GB %in% remove] %>% tibble::as_tibble_col()
# vector.GB <- vector.GB %>% tibble::as_tibble_col()
vector.GB$value <- vector.GB$value %>%
    forcats::as_factor() %>%
    plyr::revalue(c("M" = "N"))

rm(remove, list.GB)

# -------------------------------------
gather.KA.129 <- subset.KA.129.odd %>%
    dplyr::select(-sum) %>%
    gather("category", "counts", 1:9)

list.KA.129 <- vector(mode = "list", length = nrow(gather.KA.129))
for(i in 1:nrow(gather.KA.129)) {
    list.KA.129[[i]] <- if(is.na(gather.KA.129$counts[i])) {
        paste(NA_character_)
    } else if(gather.KA.129$counts[i] > 0) {
        paste(replicate(gather.KA.129$counts[i], gather.KA.129$category[i]))
    } else {
        paste("§")
    }
}

vector.KA.129 <- unlist(list.KA.129)
remove <- "§"
vector.KA.129 <- vector.KA.129[!vector.KA.129 %in% remove] %>% tibble::as_tibble_col()
# vector.KA.129 <- vector.KA.129 %>% tibble::as_tibble_col()
vector.KA.129$value <- vector.KA.129$value %>%
    forcats::as_factor() %>%
    plyr::revalue(c("M" = "N"))

rm(remove, list.KA.129)

# -------------------------------------
gather.KA.CAST <- subset.KA.CAST.odd %>%
    dplyr::select(-sum) %>%
    gather("category", "counts", 1:9)

list.KA.CAST <- vector(mode = "list", length = nrow(gather.KA.CAST))
for(i in 1:nrow(gather.KA.CAST)) {
    list.KA.CAST[[i]] <- if(is.na(gather.KA.CAST$counts[i])) {
        paste(NA_character_)
    } else if(gather.KA.CAST$counts[i] > 0) {
        paste(replicate(gather.KA.CAST$counts[i], gather.KA.CAST$category[i]))
    } else {
        paste("§")
    }
}

vector.KA.CAST <- unlist(list.KA.CAST)
remove <- "§"
vector.KA.CAST <- vector.KA.CAST[!vector.KA.CAST %in% remove] %>% tibble::as_tibble_col()
# vector.KA.CAST <- vector.KA.CAST %>% tibble::as_tibble_col()
vector.KA.CAST$value <- vector.KA.CAST$value %>%
    forcats::as_factor() %>%
    plyr::revalue(c("M" = "N"))

rm(remove, list.KA.CAST)


#  Generate confusion matrix --------------------------------------------------

#  Add N level to factor vector.KA.129$value
levels(vector.KA.129$value) <- c(levels(vector.KA.129$value), "N")

#  'data' is for predicted data, i.e., the data from my pipeline; 'reference'
#+ is for true results, i.e., the data from GB's pipeline
example <- caret::confusionMatrix(
    data = vector.GB$value,
    reference = vector.KA.129$value,
    dnn = c("GB", "KA")
)
example.table <- example$table %>%
    as.data.frame() %>%
    tibble::as_tibble()
example.table <- example.table %>%
    mutate(KA = factor(KA, rev(unique(KA))))

ggplot(example.table, aes(x = GB, y = KA, fill = Freq)) +
    geom_tile(color = "white") +
    theme_bw() +
    coord_equal() +
    scale_fill_distiller(palette = "Spectral", direction = -1) +
    guides(fill = FALSE) +
    labs(title = "For reads not NA in CAST (KA),\nconfusion matrix for N-masked\nmm10 (GB) vs. 129 (KA)") +
    geom_text(aes(label = Freq), color = "black")

rm(example, example.table)

#  'data' is for predicted data, i.e., the data from my pipeline; 'reference'
#+ is for true results, i.e., the data from GB's pipeline
example <- caret::confusionMatrix(
    data = vector.GB$value,
    reference = vector.KA.CAST$value,
    dnn = c("GB", "KA")
)
example.table <- example$table %>%
    as.data.frame() %>%
    tibble::as_tibble()
example.table <- example.table %>%
    mutate(KA = factor(KA, rev(unique(KA))))

ggplot(example.table, aes(x = GB, y = KA, fill = Freq)) +
    geom_tile(color = "white") +
    theme_bw() +
    coord_equal() +
    scale_fill_distiller(palette = "Spectral", direction = -1) +
    guides(fill = FALSE) +
    labs(title = "For reads not NA in CAST (KA),\nconfusion matrix for N-masked\nmm10 (GB) vs. CAST (KA)") +
    geom_text(aes(label = Freq), color = "black")

rm(example, example.table)

#  'data' is for predicted data, i.e., the data from my pipeline; 'reference'
#+ is for true results, i.e., the data from GB's pipeline
example <- caret::confusionMatrix(
    data = vector.KA.CAST$value,
    reference = vector.KA.129$value,
    dnn = c("KA.CAST", "KA.129")
)
example.table <- example$table %>%
    as.data.frame() %>%
    tibble::as_tibble()

ggplot(example.table, aes(x = KA.CAST, y = KA.129, fill = Freq)) +
    geom_tile(color = "white") +
    theme_bw() +
    coord_equal() +
    scale_fill_distiller(palette = "Spectral", direction = -1) +
    guides(fill = FALSE) +
    labs(title = "For reads not NA in CAST (KA),\nconfusion matrix for N-masked\nCAST (KA) vs. 129 (KA)") +
    geom_text(aes(label = Freq), color = "black")

rm(example, example.table)


#  Create 'both' columns (from 'odd' by 'even') in tibble 'bam' ---------------
bam$cigar.GB.both <- paste0(
    bam$cigar.GB.odd, " × ", bam$cigar.GB.even
)
bam$read_sequence.GB.both <- paste0(
    bam$read_sequence.GB.odd, " × ", bam$read_sequence.GB.even
)
bam$matches.GB.both <- paste0(
    bam$matches.GB.odd, " × ", bam$matches.GB.even
)
bam$reference_sequence.GB.both <- paste0(
    bam$reference_sequence.GB.odd, " × ", bam$reference_sequence.GB.even
)

bam$cigar.KA.129.both <- paste0(
    bam$cigar.KA.129.odd, " × ", bam$cigar.KA.129.even
)
bam$read_sequence.KA.129.both <- paste0(
    bam$read_sequence.KA.129.odd, " × ", bam$read_sequence.KA.129.even
)
bam$matches.KA.129.both <- paste0(
    bam$matches.KA.129.odd, " × ", bam$matches.KA.129.even
)
bam$reference_sequence.KA.129.both <- paste0(
    bam$reference_sequence.KA.129.odd, " × ", bam$reference_sequence.KA.129.even
)

bam$cigar.KA.CAST.both <- paste0(
    bam$cigar.KA.CAST.odd, " × ", bam$cigar.KA.CAST.even
)
bam$read_sequence.KA.CAST.both <- paste0(
    bam$read_sequence.KA.CAST.odd, " × ", bam$read_sequence.KA.CAST.even
)
bam$matches.KA.CAST.both <- paste0(
    bam$matches.KA.CAST.odd, " × ", bam$matches.KA.CAST.even
)
bam$reference_sequence.KA.CAST.both <- paste0(
    bam$reference_sequence.KA.CAST.odd, " × ", bam$reference_sequence.KA.CAST.even
)

# category <- category %>% dplyr::relocate(
#     c(
#         M.KA.both, M.GB.both, I.KA.both, I.GB.both, D.KA.both, D.GB.both,
#         N.KA.both, N.GB.both, S.KA.both, S.GB.both, H.KA.both, H.GB.both,
#         P.KA.both, P.GB.both, `=.KA.both`, `=.GB.both`, X.KA.both, X.GB.both
#     ),
#     .after = reference_sequence.GB.both
# )


#  #TODO Name the section -----------------------------------------------------
bam$no_differences.KA.129 <- ifelse(
    is.na(bam$matches.KA.129.odd),
    NA_character_,
    stringi::stri_count(bam$matches.KA.129.both, fixed = "*")
)
bam$no_differences.KA.CAST <- ifelse(
    is.na(bam$matches.KA.CAST.odd),
    NA_character_,
    stringi::stri_count(bam$matches.KA.CAST.both, fixed = "*")
)
bam$no_differences.GB <- ifelse(
    is.na(bam$matches.GB.odd),
    NA_character_,
    stringi::stri_count(bam$matches.GB.both, fixed = "*")
)


#  Looking at all entries in the cell -----------------------------------------
# colnames(bam)
# 
# bam$rname.129.odd %>% is.na() %>% table(useNA = "ifany")
# bam$pos.129.odd %>% is.na() %>% table(useNA = "ifany")
# bam$pos_end.129.odd %>% is.na() %>% table(useNA = "ifany")
# bam$rname.129.even %>% is.na() %>% table(useNA = "ifany")
# bam$pos.129.even %>% is.na() %>% table(useNA = "ifany")
# bam$pos_end.129.even %>% is.na() %>% table(useNA = "ifany")
# bam$cigar.KA.129.both %>% is.na() %>% table(useNA = "ifany")
# bam$read_sequence.KA.129.both %>% is.na() %>% table(useNA = "ifany")

out <- paste0(
    "\n",
    "129 (KA)", "\n",
    "read pair:", " ",
    bam$rname.129.odd, ":",
    bam$pos.129.odd, "-",
    bam$pos_end.129.odd, " × ",
    bam$rname.129.even, ":",
    bam$pos.129.even, "-",
    bam$pos_end.129.even, "\n",
    "    CIGAR:", " ",
    bam$cigar.KA.129.both, "\n",
    "     read:", " ",
    bam$read_sequence.KA.129.both, "\n",
    "  matches:", " ",
    stringr::str_replace(bam$matches.KA.129.both, "×", " "), "\n",
    "reference:", " ",
    bam$reference_sequence.KA.129.both, "\n",
    "          ", " ",
    bam$no_differences.KA.129,
    " general difference(s)", "\n\n",
    "CAST (KA)", "\n",
    "read pair:", " ",
    bam$rname.CAST.odd, ":",
    bam$pos.CAST.odd, "-",
    bam$pos_end.CAST.odd, " × ",
    bam$rname.CAST.even, ":",
    bam$pos.CAST.even, "-",
    bam$pos_end.CAST.even, "\n",
    "    CIGAR:", " ",
    bam$cigar.KA.CAST.both, "\n",
    "     read:", " ",
    bam$read_sequence.KA.CAST.both, "\n",
    "  matches:", " ",
    stringr::str_replace(bam$matches.KA.CAST.both, "×", " "), "\n",
    "reference:", " ",
    bam$reference_sequence.KA.CAST.both, "\n",
    "          ", " ",
    bam$no_differences.KA.CAST,
    " general difference(s)", "\n\n",
    "N-masked mm10 (GB)", "\n",
    "read pair:", " ",
    bam$rname.GB.odd, ":",
    bam$pos.GB.odd, "-",
    bam$pos_end.GB.odd, " × ",
    bam$rname.GB.even, ":",
    bam$pos.GB.even, "-",
    bam$pos_end.GB.even, "\n",
    "    CIGAR:", " ",
    bam$cigar.GB.both, "\n",
    "     read:", " ",
    bam$read_sequence.GB.both, "\n",
    "  matches:", " ",
    stringr::str_replace(bam$matches.GB.both, "×", " "), "\n",
    "reference:", " ",
    bam$reference_sequence.GB.both, "\n",
    "          ", " ",
    bam$no_differences.GB,
    " general difference(s)", "\n\n\n"
)
out %>% cat()

write.table(
    out,
    "2022-0214.KA.129-by-GB.Ambiguous.sample-all.txt",
    row.names = FALSE,
    col.names = FALSE
)

# category.init <- category


#  Looking at all entries with M > 0 ------------------------------------------
category <- category.init[category.init$M.GB.both > 0 | category.init$M.KA.both > 0, ]
# category <- category.init[category.init$M.KA.both > 0, ]
out <- paste0(
    "\n",
    "KA", "\n",
    "read pair:", " ",
    category$rname.129.odd, ":",
    category$pos.129.odd, "-",
    category$pos_end.129.odd, " × ",
    category$rname.129.even, ":",
    category$pos.129.even, "-",
    category$pos_end.129.even, "\n",
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

write.table(
    out,
    "2022-0214.KA.129-by-GB.Ambiguous.sample-all.M-gt-0.txt",
    row.names = FALSE,
    col.names = FALSE
)
# rm(out)
