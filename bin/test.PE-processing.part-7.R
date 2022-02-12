#!/usr/bin/env Rscript

library(BSgenome)
library(caret)
library(ggplot2)
library(Mmusculus129S1SvImJ)
library(MmusculusCASTEiJ)
library(Mmusculus129Inserted)
library(MmusculusCASTInserted)
library(MmusculusCAST129Inserted)
library(Mmusculus129Nmasked)
library(MmusculusCASTNmasked)
library(MmusculusCAST129Nmasked)
library(magrittr)
library(pheatmap)
library(Rsamtools)
library(scales)
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


#  Set up functions -----------------------------------------------------------
`%notin%` <- Negate(`%in%`)


compGenome <- function(rname, pos, pos_end, seq, row, genome) {
    chr <- rname[row] %>% as.character()
    start <- pos[row]
    end <- pos_end[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- seq[row]
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    comparison <- paste0(
        # paste0(chr, ":", start, "-", end, "\n"),
        p$string_1, "\n",
        p$result, "\n",
        p$string_2, "\n",
        p$n_diff, " difference(s)."
    )
    return(comparison)
}


compGenomeDifferences <- function(rname, pos, pos_end, seq, row, genome) {
    chr <- rname[row] %>% as.character()
    start <- pos[row]
    end <- pos_end[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- seq[row]
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    difference <- paste0(
        # paste0(chr, ":", start, "-", end, "\n"),
        # p$string_1, "\n",
        # p$result, "\n",
        # p$string_2, "\n",
        # p$n_diff, " difference(s)."
        p$n_diff
    ) %>%
        as.numeric()
    return(difference)
}


compGenomeTibble <- function(rname, pos, pos_end, seq, row, genome) {
    chr <- rname[row] %>% as.character()
    start <- pos[row]
    end <- pos_end[row]
    
    gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = start, end = end)
    )
    string_1 <<- seq[row]
    string_2 <<- genome %>%
        BSgenome::getSeq(., names = gr) %>%
        as.character()
    
    p <- reticulate::py_run_string(python_code) %>% reticulate::py_to_r()
    tibble <- tibble::tibble(
        row = row,
        chr = chr,
        start = start,
        end = end,
        seq_query = p$string_1,
        seq_genome = p$string_2,
        seq_match = p$result,
        no_difference = p$n_diff
    )
    return(tibble)
}


convertPercent <- function(x) {
    if(is.numeric(x)) {
        ifelse(is.na(x), x, paste0(round(x * 100L, 2), "%"))
    } else x
}


evaluateOperation <- function(operation = operation) {
    return(eval(parse(text = operation), envir = .GlobalEnv))
}


joinAndSelectWorst <- function(
    row, genome,
    rname.odd, pos.odd, pos_end.odd, seq.odd,
    rname.even, pos.even, pos_end.even, seq.even
) {
    # #  Assignments for tests
    # row <-  row
    # genome <- Mmusculus129Nmasked
    # rname.odd <- bam$lO_rname.129S1.odd
    # pos.odd <- bam$lO_pos.129S1.odd
    # pos_end.odd <- as.numeric(bam$lO_pos.129S1.odd + 49)
    # # pos_end.odd <- bam$lO_pos_end.129S1.odd
    # seq.odd <- bam$seq.129S1.odd
    # rname.even <- bam$lO_rname.129S1.even
    # pos.even <- bam$lO_pos.129S1.even
    # pos_end.even <- as.numeric(bam$lO_pos.129S1.even + 49)
    # # pos_end.even <- bam$lO_pos_end.129S1.even
    # seq.even <- bam$seq.129S1.even
    
    #  Generate tibble for 'odd' values
    tbl.odd <- sapply(
        row, compGenomeTibble,
        rname = rname.odd,
        pos = pos.odd,
        pos_end = pos_end.odd,
        seq = seq.odd,
        genome = genome
    ) %>%
        t() %>% 
        dplyr::as_tibble()

    #  Generate tibble for 'even' values
    tbl.even <- sapply(
        row, compGenomeTibble,
        rname = rname.even,
        pos = pos.even,
        pos_end = pos_end.even,
        seq = seq.even,
        genome = genome
    ) %>%
        t() %>% 
        dplyr::as_tibble()

    #  Join the 'odd' and 'even' tibbles by the row column
    tbl.join <- dplyr::full_join(tbl.odd, tbl.even, by = "row") %>%
        tibble::as_tibble()
    rm(tbl.even, tbl.odd)

    colnames(tbl.join) <- colnames(tbl.join) %>%
        stringr::str_replace_all("\\.x", "\\.odd") %>%
        stringr::str_replace_all("\\.y", "\\.even")

    #  For the columns, convert 'list' modes to appropriate modes, e.g.,
    #+ 'integer', 'numeric', 'character'
    tbl.join$row <- tbl.join$row %>% as.integer()
    tbl.join$chr.odd <- tbl.join$chr.odd %>% as.character()
    tbl.join$chr.even <- tbl.join$chr.even %>% as.character()
    tbl.join$start.odd <- tbl.join$start.odd %>% as.numeric()
    tbl.join$start.even <- tbl.join$start.even %>% as.numeric()
    tbl.join$end.odd <- tbl.join$end.odd %>% as.numeric()
    tbl.join$end.even <- tbl.join$end.even %>% as.numeric()
    tbl.join$seq_query.odd <- tbl.join$seq_query.odd %>% as.character()
    tbl.join$seq_query.even <- tbl.join$seq_query.even %>% as.character()
    tbl.join$seq_genome.odd <- tbl.join$seq_genome.odd %>% as.character()
    tbl.join$seq_genome.even <- tbl.join$seq_genome.even %>% as.character()
    tbl.join$seq_match.odd <- tbl.join$seq_match.odd %>% as.character()
    tbl.join$seq_match.even <- tbl.join$seq_match.even %>% as.character()
    tbl.join$no_difference.odd <- tbl.join$no_difference.odd %>% as.integer()
    tbl.join$no_difference.even <- tbl.join$no_difference.even %>% as.integer()

    #  For readability, better organization, change order of columns
    tbl.join <- tbl.join %>%
        dplyr::relocate(
            c(no_difference.odd, no_difference.even),
            .after = row
        ) %>%
        dplyr::relocate(
            c(
                seq_match.odd, seq_genome.odd, seq_query.odd,
                seq_match.even, seq_genome.even, seq_query.even
            ),
            .after = no_difference.even
        ) %>%
        dplyr::relocate(
            c(chr.odd, start.odd, end.odd, chr.even, start.even, end.even),
            .after = seq_query.even
        )

    #  To a column called 'worst', determine which of the two mates (the odd
    #+ or even) has the higher number of differences when compared to the
    #+ reference genome (for more information, see compGenomeTibble())
    tbl.join$worst <- ifelse(
        tbl.join$no_difference.odd > tbl.join$no_difference.even,
        "odd",
        ifelse(
            tbl.join$no_difference.odd < tbl.join$no_difference.even,
            "even",
            ifelse(
                tbl.join$no_difference.odd == tbl.join$no_difference.even,
                "same",
                NA_character_
            )
        )
    )
    tbl.join <- tbl.join %>%
        dplyr::relocate(worst, .after = no_difference.even)

    #  To a column called 'worst.select', assign "odd" to "odd" in 'worst',
    #+ "even" to "even" in 'worst', and randomly select "even" or "odd" for
    #+ "same" in 'worst'
    tbl.join$worst.select <- vector(
        mode = "character", length = nrow(tbl.join)
    )
    for(i in 1:nrow(tbl.join)) {
        if(tbl.join$worst[i] == "odd") {
            tbl.join$worst.select[i] <- "odd"
        } else if(tbl.join$worst[i] == "even") {
            tbl.join$worst.select[i] <- "even"
        } else if(tbl.join$worst[i] == "same") {
            tbl.join$worst.select[i] <- sample(c("odd", "even"), size = 1)
        }
        
    }
    tbl.join <- tbl.join %>%
        dplyr::relocate(worst.select, .after = worst)
    rm(i)
    
    #  Per row, add 'chr', 'start', 'end', 'seq_match', 'seq_genome',
    #+ 'seq_query', and 'no_difference' columns for the value in 'worst.select'
    tbl.join$chr.select <- vector(
        mode = "character", length = nrow(tbl.join)
    )
    tbl.join$start.select <- vector(
        mode = "numeric", length = nrow(tbl.join)
    )
    tbl.join$end.select <- vector(
        mode = "numeric", length = nrow(tbl.join)
    )
    tbl.join$seq_match.select <- vector(
        mode = "character", length = nrow(tbl.join)
    )
    tbl.join$seq_genome.select <- vector(
        mode = "character", length = nrow(tbl.join)
    )
    tbl.join$seq_query.select <- vector(
        mode = "character", length = nrow(tbl.join)
    )
    tbl.join$no_difference.select <- vector(
        mode = "integer", length = nrow(tbl.join)
    )
    for(i in 1:nrow(tbl.join)) {
        if(tbl.join$worst.select[i] == "odd") {
            tbl.join$chr.select[i] <- tbl.join$chr.odd[i]
            tbl.join$start.select[i] <- tbl.join$start.odd[i]
            tbl.join$end.select[i] <- tbl.join$end.odd[i]
            tbl.join$seq_match.select[i] <- tbl.join$seq_match.odd[i]
            tbl.join$seq_genome.select[i] <- tbl.join$seq_genome.odd[i]
            tbl.join$seq_query.select[i] <- tbl.join$seq_genome.odd[i]
            tbl.join$no_difference.select[i] <- tbl.join$no_difference.odd[i]
        } else if(tbl.join$worst.select[i] == "even") {
            tbl.join$chr.select[i] <- tbl.join$chr.even[i]
            tbl.join$start.select[i] <- tbl.join$start.even[i]
            tbl.join$end.select[i] <- tbl.join$end.even[i]
            tbl.join$seq_match.select[i] <- tbl.join$seq_match.even[i]
            tbl.join$seq_genome.select[i] <- tbl.join$seq_genome.even[i]
            tbl.join$seq_query.select[i] <- tbl.join$seq_genome.even[i]
            tbl.join$no_difference.select[i] <- tbl.join$no_difference.even[i]
        }
    }
    tbl.join <- tbl.join %>%
        dplyr::relocate(
            c(
                worst.select, no_difference.select,
                chr.select, start.select, end.select,
                seq_match.select, seq_genome.select, seq_query.select, 
            ),
            .after = row
        )
    rm(i)
    
    #  Subsection: #TODO Written explanation
    tbl.join$contains_N.select <- vector(
        mode = "logical", length = nrow(tbl.join)
    )
    tbl.join$no_N.select <- vector(
        mode = "integer", length = nrow(tbl.join)
    )
    tbl.join$no_diff_not_N.select <- vector(
        mode = "integer", length = nrow(tbl.join)
    )
    tbl.join <- tbl.join %>%
        dplyr::relocate(
            c(contains_N.select, no_N.select, no_diff_not_N.select),
            .after = no_difference.select
        )
    for(i in 1:nrow(tbl.join)) {
        tbl.join$contains_N.select[i] <- stringr::str_detect(
            tbl.join$seq_genome.select[i], "N"
        ) %>%
            unlist(., recursive = TRUE, use.names = TRUE)
        
        if(tbl.join$contains_N.select[i] == TRUE) {
            tbl.join$no_N.select[i] <-
                stringr::str_count(tbl.join$seq_genome.select[i], "N") %>%
                    as.integer()
        } else if(tbl.join$contains_N.select[i] == FALSE) {
            tbl.join$no_N.select[i] <- 0L
        }
        
        tbl.join$no_diff_not_N.select[i] <-
            tbl.join$no_difference.select[i] - tbl.join$no_N.select[i]
    }
    tbl.join <- tbl.join %>%
        dplyr::relocate(
            c(contains_N.select, no_N.select, no_diff_not_N.select),
            .after = no_difference.select
        )
    rm(i)
    
    
    #  New section, written explanation ---------------------------------------
    tbl.join$chr.both <- vector(
        mode = "character", length = nrow(tbl.join)
    )
    tbl.join$start.odd.both <- vector(
        mode = "character", length = nrow(tbl.join)
    )
    tbl.join$start.even.both <- vector(
        mode = "character", length = nrow(tbl.join)
    )
    tbl.join$end.odd.both <- vector(
        mode = "character", length = nrow(tbl.join)
    )
    tbl.join$end.even.both <- vector(
        mode = "character", length = nrow(tbl.join)
    )
    tbl.join$seq_match.both <- vector(
        mode = "character", length = nrow(tbl.join)
    )
    tbl.join$seq_genome.both <- vector(
        mode = "character", length = nrow(tbl.join)
    )
    tbl.join$seq_query.both <- vector(
        mode = "character", length = nrow(tbl.join)
    )
    tbl.join$no_difference.sum.both <- vector(
        mode = "integer", length = nrow(tbl.join)
    )
    tbl.join$no_difference.both <- vector(
        mode = "character", length = nrow(tbl.join)
    )
    for(i in 1:nrow(tbl.join)) {
        tbl.join$chr.both[i] <- tbl.join$chr.odd[i]
        tbl.join$start.odd.both[i] <- tbl.join$start.odd[i]
        tbl.join$start.even.both[i] <- tbl.join$start.even[i]
        tbl.join$end.odd.both[i] <- tbl.join$end.odd[i]
        tbl.join$end.even.both[i] <- tbl.join$end.even[i]
        tbl.join$seq_match.both[i] <- paste0(
            tbl.join$seq_match.odd[i], " × ", tbl.join$seq_match.even[i]
        )
        tbl.join$seq_genome.both[i] <- paste0(
            tbl.join$seq_genome.odd[i], " × ", tbl.join$seq_genome.even[i]
        )
        tbl.join$seq_query.both[i] <- paste0(
            tbl.join$seq_query.odd[i], " × ", tbl.join$seq_query.even[i]
        )
        tbl.join$no_difference.sum.both[i] <-
            tbl.join$no_difference.odd[i] + tbl.join$no_difference.even[i]
        tbl.join$no_difference.both[i] <- paste0(
            paste0(tbl.join$no_difference.odd[i]),
            " + ",
            paste0(tbl.join$no_difference.even[i])
        )
    }
    tbl.join <- tbl.join %>%
        dplyr::relocate(
            c(
                no_difference.sum.both, no_difference.both,
                seq_match.both, seq_genome.both, seq_query.both,
                chr.both,
                start.odd.both, start.even.both,
                end.odd.both, end.even.both
            ),
            .after = row
        )
    rm(i)
    
    #  Subsection: #TODO Written explanation
    tbl.join$contains_N.both <- vector(
        mode = "logical", length = nrow(tbl.join)
    )
    tbl.join$no_N.both <- vector(
        mode = "integer", length = nrow(tbl.join)
    )
    tbl.join$no_diff_not_N.both <- vector(
        mode = "integer", length = nrow(tbl.join)
    )
    for(i in 1:nrow(tbl.join)) {
        tbl.join$contains_N.both[i] <- stringr::str_detect(
            tbl.join$seq_genome.both[i], "N"
        ) %>%
            unlist(., recursive = TRUE, use.names = TRUE)
        
        if(tbl.join$contains_N.both[i] == TRUE) {
            tbl.join$no_N.both[i] <-
                stringr::str_count(tbl.join$seq_genome.both[i], "N") %>%
                    as.integer()
        } else if(tbl.join$contains_N.both[i] == FALSE) {
            tbl.join$no_N.both[i] <- 0L
        }
        
        tbl.join$no_diff_not_N.both[i] <-
            tbl.join$no_difference.sum.both[i] - tbl.join$no_N.both[i]
    }
    tbl.join <- tbl.join %>%
        dplyr::relocate(
            c(contains_N.both, no_N.both, no_diff_not_N.both),
            .after = no_difference.both
        )
    rm(i)
    
    return(tbl.join)
}


loadRDS <- readRDS


loadTibbleFromRDS <- function(variable, file) {
    command <- paste0(variable, " <- loadRDS(file = \"", file, "\")") %>% as.list()
    return(eval(parse(text = command), envir = globalenv()))
}


makeOperation <- function(variable, command) {
    operation <- paste(variable, command)
    return(operation)
}


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


#  Load in the KA and GB datasets ---------------------------------------------

#  Having run parts 5 and 6 first, load *.Rdata files for them
load("test.PE-processing.part-5.Rdata")  # KA
load("test.PE-processing.part-6.Rdata")  # GB

# chromosome <- "chr1"
chromosome <- "chrX"

#  Clean up unneeded variables, etc. from the environment
operation <- paste0(
    "rm(",
        variable, ", ",
        variable_b, ", ",
        variable_m, ", ",
        variable_uniline,
    ")"
)
evaluateOperation(operation)

rm(suffix, variable, variable_b, variable_m, variable_uniline)


#  Script name ----------------------------------------------------------------
script <- "test.PE-processing.part-7.R"


#  Rename variables to make things clearer ------------------------------------
bi.AS.pmin <- b.full
uni.AS.pmin <- AS.pmin

bi.joint.GB <- joint.GB
uni.joint.GB <- uniline.joint.GB

rm(AS.pmin, b.full, joint.GB, uniline.joint.GB)


#  Create factor levels for NA entries ----------------------------------------
uni.AS.pmin$assignment <- uni.AS.pmin$assignment %>% as.character()
uni.AS.pmin$assignment[is.na(uni.AS.pmin$assignment)] <- "NA"
uni.AS.pmin$assignment <- uni.AS.pmin$assignment %>% forcats::as_factor()
uni.AS.pmin$assignment %>% levels()

uni.AS.pmin$assignment <- uni.AS.pmin$assignment %>%
    plyr::revalue(c("Neutral" = "Ambiguous"))


#  Deconcatenate the uni.AS.pmin tibble ----------------------------------------
df <- uni.AS.pmin %>%
    dplyr::select(
        -trinary, -trinary.r,
        -assignment_trinary, -assignment_trinary.r
    )

df <- df %>% dplyr::mutate(
    AS.mm10.odd = AS.mm10,
    AS.129S1.odd = AS.129S1,
    AS.CAST.odd = AS.CAST,
    assignment.odd = assignment,
    AS.mm10.even = AS.mm10,
    AS.129S1.even = AS.129S1,
    AS.CAST.even = AS.CAST,
    assignment.even = assignment
) %>%
    dplyr::select(-AS.mm10, -AS.129S1, -AS.CAST, -assignment) %>%
    dplyr::relocate(
        c(
            AS.mm10.odd, AS.129S1.odd, AS.CAST.odd, assignment.odd,
            AS.mm10.even, AS.129S1.even, AS.CAST.even, assignment.even
        ),
        .before = coordinate.odd
    )

#  Split the uniline tibbles by odd or even status
odd <- df[stringr::str_subset(colnames(df), "\\.odd")]
even <- df[stringr::str_subset(colnames(df), "\\.even")]

#  Strip suffixes from column names
colnames(odd) <- str_replace_all(colnames(odd), "\\.odd", "")
colnames(even) <- str_replace_all(colnames(even), "\\.even", "")

#  Odd tibble is tibble #1, even tibble is tibble #2
odd$tibble <- "1"
even$tibble <- "2"

#  Again, interleave the odd/even rows, then arrange them by group ID and
#+ tibble number; save results to initial "bi.AS.pmin" variable
command <- paste0(
    "<- ", "odd", " %>%",
        "dplyr::mutate(groupid = row_number())", " %>% ",
        "dplyr::bind_rows(even %>% mutate(groupid = row_number()))", " %>% ",
        "dplyr::arrange(groupid, tibble)", " %>% ",
        "dplyr::select(-tibble)"
)
operation <- makeOperation(paste0("bi.AS.pmin"), command)
evaluateOperation(operation)

uni.AS.pmin <- df

# #  Again, check to make sure that there are no more than two entries per
# #+ group ID
# n_occur <- data.frame(table(eval(parse(text = paste0("AS.pmin")))$groupid))
# print(n_occur[n_occur$Freq > 2, ])

#  Clean up
rm(df, n_occur, even, odd)


###############################################################################
names(uni.joint.GB)[names(uni.joint.GB) == "coordinate.GB.even"] <-
    "coordinate.even"
names(uni.joint.GB)[names(uni.joint.GB) == "coordinate.GB.odd"] <-
    "coordinate.odd"

joint <- dplyr::full_join(
    uni.AS.pmin, uni.joint.GB, by = "coordinate.odd"
)
(joint$coordinate.odd %notin% uni.AS.pmin$coordinate.odd) %>%
    table(useNA = "ifany")
(joint$coordinate.odd %notin% uni.joint.GB$coordinate.odd) %>%
    table(useNA = "ifany")
(joint$coordinate.even.x == joint$coordinate.even.y) %>%
    table(useNA = "ifany")

names(joint)[names(joint) == "coordinate.even.x"] <- "coordinate.KA.even"
names(joint)[names(joint) == "coordinate.even.y"] <- "coordinate.GB.even"
names(joint)[names(joint) == "coordinate.odd"] <- "coordinate.KA.odd"

names(joint)[names(joint) == "assignment.odd"] <- "assignment.KA.odd"
names(joint)[names(joint) == "assignment.even"] <- "assignment.KA.even"

joint$coordinate.GB.odd <- joint$coordinate.KA.odd
joint <- joint %>%
    dplyr::relocate(coordinate.GB.odd, .after = qmpos.GB.odd) %>%
    dplyr::select(-discrepancy.GB.odd, -discrepancy.GB.even)

joint$assignment.KA.odd <- paste0("KA.", joint$assignment.KA.odd)
joint$assignment.GB.odd[is.na(joint$assignment.GB.odd)] <- "GB.NA"
joint$assignment.GB.even[is.na(joint$assignment.GB.even)] <- "NA"

joint$assignment.KA.even <- joint$assignment.KA.even %>% as.character()

order <- c("129", "CAST", "Ambiguous", "NA")
joint$assignment.KA.even <- joint$assignment.KA.even %>%
    stringr::str_remove("-EiJ") %>%
    stringr::str_remove("S1-SvImJ") %>%
    factor(., levels = order)
joint$assignment.KA.odd <- joint$assignment.KA.odd %>%
    stringr::str_remove("-EiJ") %>%
    stringr::str_remove("S1-SvImJ") %>%
    factor(., levels = paste0("KA.", order))

order <- c("129", "CAST", "Ambiguous", "Contra", "NA")
joint$assignment.GB.even <- joint$assignment.GB.even %>%
    stringr::str_remove("GB.") %>%
    stringr::str_remove("S1") %>%
    factor(., levels = order)
joint$assignment.GB.odd <- joint$assignment.GB.odd %>%
    stringr::str_remove("S1") %>%
    factor(., levels = paste0("GB.", order))

# #  Check
# joint$assignment.KA.odd %>% head()
# joint$assignment.KA.odd %>% tail()
# joint$assignment.KA.odd %>% levels()
# 
# joint$assignment.KA.even %>% head()
# joint$assignment.KA.even %>% tail()
# joint$assignment.KA.even %>% levels()
# 
# joint$assignment.GB.odd %>% head()
# joint$assignment.GB.odd %>% tail()
# joint$assignment.GB.odd %>% levels()
# 
# joint$assignment.GB.even %>% head()
# joint$assignment.GB.even %>% tail()
# joint$assignment.GB.even %>% levels()

order <- c(
    "KA.129 × GB.Ambiguous",
    "KA.129 × GB.129",
    "KA.129 × GB.CAST",
    "KA.129 × GB.Contra",
    "KA.129 × GB.NA",
    "KA.CAST × GB.129",
    "KA.CAST × GB.CAST",
    "KA.CAST × GB.Ambiguous",
    "KA.CAST × GB.NA",
    "KA.Ambiguous × GB.129",
    "KA.Ambiguous × GB.CAST",
    "KA.Ambiguous × GB.Ambiguous",
    "KA.Ambiguous × GB.Contra",
    "KA.Ambiguous × GB.NA",
    "KA.NA × GB.129",
    "KA.NA × GB.CAST",
    "KA.NA × GB.Ambiguous",
    "KA.NA × GB.Contra"
)
joint$KA_by_GB <- paste0(
    joint$assignment.KA.odd, " × ", joint$assignment.GB.odd
) %>%
    factor(levels = order)

order <- c(
    "GB.129 × KA.129",
    "GB.129 × KA.CAST",
    "GB.129 × KA.Ambiguous",
    "GB.129 × KA.NA",
    "GB.CAST × KA.129",
    "GB.CAST × KA.CAST",
    "GB.CAST × KA.Ambiguous",
    "GB.CAST × KA.NA",
    "GB.Ambiguous × KA.129",
    "GB.Ambiguous × KA.CAST",
    "GB.Ambiguous × KA.Ambiguous",
    "GB.Ambiguous × KA.NA",
    "GB.Contra × KA.129",
    "GB.Contra × KA.Ambiguous",
    "GB.Contra × KA.NA",
    "GB.NA × KA.129",
    "GB.NA × KA.CAST",
    "GB.NA × KA.Ambiguous"
)
joint$GB_by_KA <- paste0(
    joint$assignment.GB.odd, " × ", joint$assignment.KA.odd
) %>%
    factor(levels = order)

# #  Check
# joint$KA_by_GB %>% head()
# joint$GB_by_KA %>% head()
# joint$KA_by_GB %>% levels()
# joint$GB_by_KA %>% levels()

#  For increased readbility, reorder columns
joint <- joint %>%
    dplyr::relocate(c(KA_by_GB, GB_by_KA), .before = AS.mm10.odd) %>%
    dplyr::relocate(
        c(
            rname.mm10.odd, pos.mm10.odd, pos_end.mm10.odd,
            rname.mm10.even, pos.mm10.even, pos_end.mm10.even,
            rname.129S1.odd, pos.129S1.odd, pos_end.129S1.odd,
            rname.129S1.even, pos.129S1.even, pos_end.129S1.even,
            rname.CAST.odd, pos.CAST.odd, pos_end.CAST.odd,
            rname.CAST.even, pos.CAST.even, pos_end.CAST.even,
            lO_rname.129S1.odd, lO_pos.129S1.odd, lO_pos_end.129S1.odd,
            lO_rname.129S1.even, lO_pos.129S1.even, lO_pos_end.129S1.even,
            lO_rname.CAST.odd, lO_pos.CAST.odd, lO_pos_end.CAST.odd,
            lO_rname.CAST.even, lO_pos.CAST.even, lO_pos_end.CAST.even,
            mrnm.mm10.odd, mpos.mm10.odd, mpos_end.mm10.odd,
            mrnm.mm10.even, mpos.mm10.even, mpos_end.mm10.even,
            mrnm.129S1.odd, mpos.129S1.odd, mpos_end.129S1.odd,
            mrnm.129S1.even, mpos.129S1.even, mpos_end.129S1.even,
            mrnm.CAST.odd, mpos.CAST.odd, mpos_end.CAST.odd,
            mrnm.CAST.even, mpos.CAST.even, mpos_end.CAST.even,
            lO_mrnm.129S1.odd, lO_mpos.129S1.odd, lO_mpos_end.129S1.odd,
            lO_mrnm.129S1.even, lO_mpos.129S1.even, lO_mpos_end.129S1.even,
            lO_mrnm.CAST.odd, lO_mpos.CAST.odd, lO_mpos_end.CAST.odd,
            lO_mrnm.CAST.even, lO_mpos.CAST.even, lO_mpos_end.CAST.even
        ),
        .after = seq.CAST.even
    ) %>%
    dplyr::relocate(
        c(
            assignment.GB.odd, assignment.GB.even,
            coordinate.GB.odd, coordinate.GB.even,
            flag.GB.odd, flag.GB.even,
            seq.GB.odd, seq.GB.even,
            rname.GB.odd, pos.GB.odd, pos_end.GB.odd,
            rname.GB.even, pos.GB.even, pos_end.GB.even,
            mrnm.GB.odd, mpos.GB.odd, mpos_end.GB.odd,
            mrnm.GB.even, mpos.GB.even, mpos_end.GB.even,
            qname.GB.odd, qpos.GB.odd, qmpos.GB.odd,
            qname.GB.even, qpos.GB.even, qmpos.GB.even,
            mapq.GB.odd, mapq.GB.even,
            cigar.GB.odd, cigar.GB.even,
            isize.GB.odd, isize_abs.GB.odd,
            isize.GB.even, isize_abs.GB.even,
            qual.GB.odd, qual.GB.even,
            groupid.GB.odd, groupid.GB.even
        ),
        .after = lO_mpos_end.CAST.even
    )
# colnames(joint)


#  Generate confusion matrix --------------------------------------------------

#  'data' is for predicted data, i.e., the data from my pipeline; 'reference'
#+ is for true results, i.e., the data from GB's pipeline
example <- confusionMatrix(
    data = joint$assignment.KA.even,
    reference = joint$assignment.GB.even,
    dnn = c("KA", "GB")
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
    labs(title = "Value distribution") +
    geom_text(aes(label = Freq), color = "black")

rm(example, example.table)


###############################################################################
#  Read in sam2pairwise files in order to reconstruct CIGAR sequences ---------

# Switch to location of sam2pairwise files for GB .bam files
setwd(paste0(path.1, "/", path.2))

#TODO 1/4 Work with processed.mate-paired.mm10.chrX.sam2pairwise.munged.txt,
#TODO 2/4 processed.mate-paired.CAST.chrX.sam2pairwise.munged.txt, and
#TODO 3/4 processed.mate-paired.129S1.chrX.sam2pairwise.munged.txt, not the
#TODO 4/4 "\\.sam2pairwise.munged.all.txt$" files as seen in the code below

#  sam2pairwise files for .bam files output by GB's allele assignment script
file <- list.files(pattern = paste0("\\.sam2pairwise.munged.all.txt$"))
variable <- file %>%
    strsplit(., "\\.") %>%
    lapply(., `[[`, 4) %>% unlist() %>%
    paste0("sam2pairwise.GB.", .)
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

command <- paste0("<- readr::read_tsv(", variable, ")")
operation <- makeOperation(variable, command)
evaluateOperation(operation)

# #  The .bam files themselves, loading in AS and MD fields
# file <- list.files(pattern = paste0("\\.bam$"))
# variable <- file %>%
#     strsplit(., "\\.") %>%
#     lapply(., `[[`, 4) %>% unlist() %>%
#     paste0("bam.GB.", .)
# mapply(
#     assign, variable, file, MoreArgs = list(envir = parent.frame())
# )
# 
# map_params <- Rsamtools::ScanBamParam(tag = c("AS", "MD"))
# command <- paste0(
#     "<- ", variable, " %>% ",
#         "Rsamtools::BamFile(., index = ", index, ", asMates = TRUE)", " %>% ",
#         "Rsamtools::scanBam(param = map_params)", " %>% ",
#         "as.data.frame()", " %>% ",
#         "tibble::as_tibble()"
# )
# operation <- makeOperation(variable, command)
# evaluateOperation(operation)
# rm(map_params)

#  Save the sam2pairwise variables for GB .bam files
variable.tmp <- variable


#  Switch to location of sam2pairwise files for KA .bam files
setwd(paste0(path.1, "/", path.3))

#  sam2pairwise files for .bam files output by KA's allele assignment script
file <- list.files(pattern = paste0("\\.rmdup.extendedCIGAR.sam2pairwise.munged.txt$"))
variable <- file %>%
    strsplit(., "-") %>%
    lapply(., `[[`, 1) %>% unlist() %>%
    paste0("sam2pairwise.KA.", .)
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

command <- paste0("<- readr::read_tsv(", variable, ")")
operation <- makeOperation(variable, command)
evaluateOperation(operation)

variable.tmp.2 <- variable

# #  Load in associated .bam information
# # chromosome <- "chr1"
# chromosome <- "chrX"
# 
# #  .bam files
# file <- list.files(pattern = paste0(
#     "\\", chromosome, ".rmdup.extendedCIGAR.bam$"
# ))
# file <- file %>% stringr::str_subset("mm10\\.", negate = TRUE)
# variable <- file %>%
#     strsplit(., "-") %>%
#     lapply(., `[[`, 1) %>% unlist() %>%
#     paste0("dedup.", .)
# mapply(
#     assign, variable, file, MoreArgs = list(envir = parent.frame())
# )
# 
# #  .bai files
# file <- list.files(pattern = paste0(
#     "\\", chromosome, ".rmdup.extendedCIGAR.bam.bai$"
# ))
# file <- file %>% stringr::str_subset("mm10\\.", negate = TRUE)
# index <- file %>%
#     strsplit(., "-") %>%
#     lapply(., `[[`, 1) %>% unlist() %>%
#     paste0("index.", .)
# mapply(
#     assign, index, file, MoreArgs = list(envir = parent.frame())
# )
# 
# #  Load in standard .bam fields
# command <- paste0(
#     "<- ", variable, " %>% ",
#         "Rsamtools::BamFile(., index = ", index, ", asMates = TRUE)", " %>% ",
#         "Rsamtools::scanBam()", " %>% ",
#         "as.data.frame()", " %>% ",
#         "tibble::as_tibble()"
# )
# operation <- makeOperation(paste0(variable, ".full"), command)
# evaluateOperation(operation)
# 
# #  Load in .bam AS and MD fields
# map_params <- Rsamtools::ScanBamParam(tag = c("AS", "MD"))
# command <- paste0(
#     "<- ", variable, " %>% ",
#         "Rsamtools::BamFile(., index = ", index, ", asMates = TRUE)", " %>% ",
#         "Rsamtools::scanBam(param = map_params)", " %>% ",
#         "as.data.frame()", " %>% ",
#         "tibble::as_tibble()"
# )
# operation <- makeOperation(variable, command)
# evaluateOperation(operation)
# rm(map_params)
# 
# #  Join the standard, AS, and MD fields
# command <- paste0(
#     "<- dplyr::bind_cols(", variable, ".full, ", variable, ")"
# )
# operation <- makeOperation(variable, command)
# evaluateOperation(operation)
# 
# #  Remove variables no longer needed
# operation <- paste0("rm(", variable, ".full, ", index, ")")
# evaluateOperation(operation)
# 
# rm(file, index)
# 
# sam2pairwise.KA.129S1 <- dplyr::left_join(
#     sam2pairwise.KA.129S1,
#     dedup.129S1,
#     by
# )




#  Create a variable for the sam2pairwise variables for GB and KA .bam files
variable <- c(variable.tmp, variable.tmp.2)
rm(variable.tmp, variable.tmp.2)


#  From sam2pairwise variables, filter out rows with mapq less than 30 --------
command <- paste0("<- ", variable, "[", variable, "$mapq >= 30, ]")
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Prior to sorting, deduplicate and mate-label the sam2pairwise tibbles ------

#  Set up $criteria, a variable needed for sorting: qname, flag, pos, mpos
command <- paste0(
    "<- paste0(",
        variable, "$qname, ", "'_', ",
        variable, "$flag, ", "'_', ",
        variable, "$pos, ", "'_', ",
        variable, "$mpos",
    ")"
)
operation <- makeOperation(paste0(variable, "$criteria"), command)
evaluateOperation(operation)

#  Set up $coordinate, a variable needed for sorting too: qname, rname, pos
command <- paste0(
    "<- paste0(",
            variable, "$qname, ", "'_', ",
            variable, "$rname, ", "'_', ",
            variable, "$pos",
    ")"
)
operation <- makeOperation(paste0(variable, "$coordinate"), command)
evaluateOperation(operation)

#  Set up $qpos, a variable needed for sorting
command <- paste0(
    "<- paste0(",
        variable, "$qname, ", "'_', ",
        variable, "$pos",
    ")"
)
operation <- makeOperation(paste0(variable, "$qpos"), command)
evaluateOperation(operation)

#  Set up $qmpos, a variable needed for sorting
command <- paste0(
    "<- paste0(",
        variable, "$qname, ", "'_', ",
        variable, "$mpos",
    ")"
)
operation <- makeOperation(paste0(variable, "$qmpos"), command)
evaluateOperation(operation)

#  Set up $isize_abs, i.e., absolute insert-size values, for sorting
command <- paste0(
    "<- ", variable, "$isize %>% abs() %>% as.numeric()"
)
operation <- makeOperation(paste0(variable, "$isize_abs"), command)
evaluateOperation(operation)

#  To survey duplicates, split out entries in which $qpos is the same as
#+ $qmpos, i.e., entries in which each member of the pair maps to the same
#+ location
command <- paste0(
    "<- ", variable, "[", variable, "$qpos == ", variable, "$qmpos, ]"
)
operation <- makeOperation(paste0("same.", variable), command)
evaluateOperation(operation)

#  Regarding the main tibbles, remove entries in which $qpos is the same as
#+ $qmpos, i.e., entries in which each member of the pair maps to the same
#+ location
command <- paste0(
    "<- ", variable, "[!(", variable, "$qpos == ", variable, "$qmpos), ]"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  Remove variables entries in which $qpos is the same as $qmpos, i.e., entries
#+ in which each member of the pair maps to the same location
operation <- paste0("rm(", "same.", variable, ")")
evaluateOperation(operation)

# #  Order tibble by $criteria and $mapq
# command <- paste0(
#     "<- ", variable, "[",
#         "order(",
#             variable, "$criteria, ",
#             variable, "$mapq",
#         "), ",
#     "]"
# )
# operation <- makeOperation(variable, command)
# evaluateOperation(operation)

#QUESTION? Does the below leave one of the entries post deduplication? Yes.
#+
#+ Deduplicate the main tibbles based on $criteria: There should be no more
#+ than one entry with a given $criteria value (combination of $qname, $flag,
#+ $pos, and $mpos)
command <- paste0(
    "<- ", variable, "[!duplicated(", variable, "$criteria), ]"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


# #  Is $qpos in $qmpos? Vice versa? Create tables and tibbles of answers -------
# createTibbleFromTibbles <- function(vector_string_tibbles) {
#     bind_rows(
#         eval(parse(text = paste0(vector_string_tibbles)[1])),
#         eval(parse(text = paste0(vector_string_tibbles)[2])),
#         eval(parse(text = paste0(vector_string_tibbles)[3])),
#         eval(parse(text = paste0(vector_string_tibbles)[4])),
#         eval(parse(text = paste0(vector_string_tibbles)[5])),
#         eval(parse(text = paste0(vector_string_tibbles)[6])),
#         eval(parse(text = paste0(vector_string_tibbles)[7]))
#     )
# }
# 
# command <- paste0(
#     "<- testPosInMpos(",
#         "pos = ", variable, "$qpos, ",
#         "mpos = ", variable, "$qmpos",
#     ") %>% ",
#         "table()"
# )
# operation <- makeOperation(paste0("PIM.", variable), command)
# evaluateOperation(operation)
# 
# command <- paste0(
#     "<- testMposInPos(",
#         "pos = ", variable, "$qpos, ",
#         "mpos = ", variable, "$qmpos",
#     ") %>% ",
#         "table()"
# )
# operation <- makeOperation(paste0("MIP.", variable), command)
# evaluateOperation(operation)
# 
# #  Convert the tables to tibbles: PIM
# command <- paste0(
#     "<- PIM.", variable, " %>% ",
#         "t() %>% ",
#         "cbind() %>% ",
#         "as_tibble()"
# )
# operation <- makeOperation(paste0("PIM.", variable), command)
# evaluateOperation(operation)
# 
# #  Create a "tibble of tibbles:" PIM
# command <- paste0(
#     "<- ", "tibble::tibble(", "PIM.", variable, ")"
# )
# operation <- makeOperation(paste0("PIM.", variable), command)
# evaluateOperation(operation)
# 
# vector_string_tibbles <- paste0("PIM.", variable)
# z.check.1.PIM <- createTibbleFromTibbles(vector_string_tibbles)
# z.check.1.PIM$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
# z.check.1.PIM <- z.check.1.PIM %>% dplyr::arrange(rownames)
# 
# #  Convert the tables to tibbles: MIP
# command <- paste0(
#     "<- MIP.", variable, " %>% ",
#         "t() %>% ",
#         "cbind() %>% ",
#         "as_tibble()"
# )
# operation <- makeOperation(paste0("MIP.", variable), command)
# evaluateOperation(operation)
# 
# #  Create a "tibble of tibbles:" MIP
# command <- paste0(
#     "<- ", "tibble::tibble(", "MIP.", variable, ")"
# )
# operation <- makeOperation(paste0("MIP.", variable), command)
# evaluateOperation(operation)
# 
# vector_string_tibbles <- paste0("MIP.", variable)
# z.check.1.MIP <- createTibbleFromTibbles(vector_string_tibbles)
# z.check.1.MIP$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
# z.check.1.MIP <- z.check.1.MIP %>% dplyr::arrange(rownames)
# 
# 
# #  Is PIM.*$`FALSE` == MIP.*$`FALSE`? -----------------------------------------
# command <- paste0(
#     "<- PIM.", variable, "$`FALSE` == MIP.", variable, "$`FALSE`"
# )
# operation <- makeOperation(paste0("PeqM.", variable), command)
# evaluateOperation(operation)
# 
# #  Create a "tibble of tibbles"
# command <- paste0(
#     "<- ", "tibble::tibble(", "PeqM.", variable, ")"
# )
# operation <- makeOperation(paste0("PeqM.", variable), command)
# evaluateOperation(operation)
# 
# operation <- paste0("colnames(", "PeqM.", variable, ") ", "<- ", "\"n\"")
# evaluateOperation(operation)
# 
# vector_string_tibbles <- paste0("PeqM.", variable)
# z.check.2.PeqM <- createTibbleFromTibbles(vector_string_tibbles)
# z.check.2.PeqM$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
# z.check.2.PeqM <- z.check.2.PeqM %>% dplyr::arrange(rownames)
# 
# 
# #  Is PIM.*$`FALSE` > MIP.*$`FALSE`? ------------------------------------------
# command <- paste0(
#     "<- PIM.", variable, "$`FALSE` > MIP.", variable, "$`FALSE`"
# )
# operation <- makeOperation(paste0("PgtM.", variable), command)
# evaluateOperation(operation)
# 
# #  Create a tibble
# command <- paste0(
#     "<- ", "tibble::tibble(", "PgtM.", variable, ")"
# )
# operation <- makeOperation(paste0("PgtM.", variable), command)
# evaluateOperation(operation)
# 
# operation <- paste0("colnames(", "PgtM.", variable, ") ", "<- ", "\"n\"")
# evaluateOperation(operation)
# 
# vector_string_tibbles <- paste0("PgtM.", variable)
# z.check.2.PgtM <- createTibbleFromTibbles(vector_string_tibbles)
# z.check.2.PgtM$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
# z.check.2.PgtM <- z.check.2.PgtM %>% dplyr::arrange(rownames)
# 
# 
# #  Is PIM.*$`FALSE` < MIP.*$`FALSE`? ------------------------------------------
# command <- paste0(
#     "<- PIM.", variable, "$`FALSE` < MIP.", variable, "$`FALSE`"
# )
# operation <- makeOperation(paste0("PltM.", variable), command)
# evaluateOperation(operation)
# 
# #  Create a tibble
# command <- paste0(
#     "<- ", "tibble::tibble(", "PltM.", variable, ")"
# )
# operation <- makeOperation(paste0("PltM.", variable), command)
# evaluateOperation(operation)
# 
# operation <- paste0("colnames(", "PltM.", variable, ") ", "<- ", "\"n\"")
# evaluateOperation(operation)
# 
# vector_string_tibbles <- paste0("PltM.", variable)
# z.check.2.PltM <- createTibbleFromTibbles(vector_string_tibbles)
# z.check.2.PltM$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
# z.check.2.PltM <- z.check.2.PltM %>% dplyr::arrange(rownames)
# 
# #  Remove unneeded variables
# command <- paste0(
#     "rm(",
#         "PIM.", variable, ", ",
#         "MIP.", variable, ", ",
#         "PeqM.", variable, ", ",
#         "PltM.", variable, ", ",
#         "PgtM.", variable,
#     ")"
# )
# evaluateOperation(command)
# rm(vector_string_tibbles)
# 
# 
# #  Make logical vectors for $qpos in $qmpos and vice versa --------------------
# 
# #  PIM
# command <- paste0(
#     "<- testPosInMpos(pos = ", variable, "$qpos, mpos = ", variable, "$qmpos)"
# )
# operation <- makeOperation(paste0("vPIM.", variable), command)
# evaluateOperation(operation)
# 
# #  MIP
# command <- paste0(
#     "<- testMposInPos(pos = ", variable, "$qpos, mpos = ", variable, "$qmpos)"
# )
# operation <- makeOperation(paste0("vMIP.", variable), command)
# evaluateOperation(operation)
# 
# #NOTE
# #  Based on PeqM.*, PgtM.*, and PltM.* results, use the following vectors:
# #+ All vPIM.* are fine, apparently
# 
# #  Remove unneeded variables
# command <- paste0(
#     "rm(", "vMIP.", variable, ")"
# )
# evaluateOperation(command)
# 
# 
# #  Subset tibbles based on logical vectors for $qpos in $qmpos ----------------
# 
# #NOTE vPIM: logical vectors for $qpos %in% $qmpos
# command <- paste0("<- ", variable, "[", "vPIM.", variable, ", ]")
# operation <- makeOperation(paste0("sub.", variable), command)
# evaluateOperation(operation)
# 
# #  Remove unneeded vPIM vectors
# command <- paste0("rm(", "vPIM.", variable, ")")
# evaluateOperation(command)
# 
# 
# #  For the subsetted tibbles, save qpos %in$ qmpos, vice versa in vectors -----
# variable <- paste0("sub.", variable)
# command <- paste0(  # Table qpos %in$ qmpos
#     "<- testPosInMpos(",
#         "pos = ", variable, "$qpos, ",
#         "mpos = ", variable, "$qmpos",
#     ") %>% ",
#         "table()"
# )
# operation <- makeOperation(paste0("PIM.", variable), command)
# evaluateOperation(operation)
# 
# command <- paste0(  # Table qmpos %in$ qpos
#     "<- testMposInPos(",
#         "pos = ", variable, "$qpos, ",
#         "mpos = ", variable, "$qmpos",
#     ") %>% ",
#         "table()"
# )
# operation <- makeOperation(paste0("MIP.", variable), command)
# evaluateOperation(operation)
# 
# #  Convert the tables to tibbles
# command <- paste0(
#     "<- PIM.", variable, " %>% ",
#         "t() %>% ",
#         "cbind() %>% ",
#         "as_tibble()"
# )
# operation <- makeOperation(paste0("PIM.", variable), command)
# evaluateOperation(operation)
# 
# command <- paste0(
#     "<- MIP.", variable, " %>% ",
#         "t() %>% ",
#         "cbind() %>% ",
#         "as_tibble()"
# )
# operation <- makeOperation(paste0("MIP.", variable), command)
# evaluateOperation(operation)
# 
# #  Collect pertinent readouts for the subsetted tibbles in one tibble; i.e.,
# #+ create a "tibble of tibbles"
# 
# #  PIM
# vector_string_tibbles <- paste0("PIM.", variable)
# z.check.sub.0.PIM <- createTibbleFromTibbles(vector_string_tibbles)
# z.check.sub.0.PIM$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
# z.check.sub.0.PIM <- z.check.sub.0.PIM %>% dplyr::arrange(rownames)
# 
# #  Remove unneeded variables
# command <- paste0("rm(", "PIM.", variable, ")")
# evaluateOperation(command)
# rm(vector_string_tibbles)
# 
# #  MIP
# vector_string_tibbles <- paste0("MIP.", variable)
# z.check.sub.0.MIP <- createTibbleFromTibbles(vector_string_tibbles)
# z.check.sub.0.MIP$rownames <- stringr::str_remove(vector_string_tibbles, "tbl.")
# z.check.sub.0.MIP <- z.check.sub.0.MIP %>% dplyr::arrange(rownames)
# 
# #  Remove unneeded variables
# command <- paste0("rm(", "MIP.", variable, ")")
# evaluateOperation(command)
# rm(vector_string_tibbles)
# 
# 
# #  Order the tibbles based on $isize_abs, $qpos, and $qmpos -------------------
# command <- paste0(
#     "<- ", variable, "[",
#         "order(",
#             variable, "$isize_abs, ",
#             variable, "$qpos, ",
#             variable, "$qmpos",
#         "), ",
#     "]"
# )
# operation <- makeOperation(variable, command)
# evaluateOperation(operation)


# -----------------------------------------------------------------------------
test <- joint


# -----------------------------------------------------------------------------
sam2pairwise.GB <- dplyr::bind_rows(
    sam2pairwise.GB.alt,
    sam2pairwise.GB.ambig,
    sam2pairwise.GB.contra,
    sam2pairwise.GB.ref
)


# -----------------------------------------------------------------------------
getSam2pairwiseGBodd <- function(tibble) {
    tbl <- tibble %>%
    dplyr::select(
        coordinate, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        coordinate.GB.odd = coordinate,
        cigar.GB.odd = cigar,
        read_sequence.GB.odd = read_sequence,
        matches.GB.odd = matches,
        reference_sequence.GB.odd = reference_sequence
    )
    return(tbl)
}


getSam2pairwiseGBeven <- function(tibble) {
    tbl <- tibble %>%
    dplyr::select(
        coordinate, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        coordinate.GB.even = coordinate,
        cigar.GB.even = cigar,
        read_sequence.GB.even = read_sequence,
        matches.GB.even = matches,
        reference_sequence.GB.even = reference_sequence
    )
    return(tbl)
}

test <- dplyr::left_join(
    test,
    getSam2pairwiseGBodd(sam2pairwise.GB),
    by = "coordinate.GB.odd"
) %>% dplyr::rename(
        cigar.GB.odd.2 = cigar.GB.odd.y
    ) %>% dplyr::rename(
        cigar.GB.odd.1 = cigar.GB.odd.x
    )

test <- dplyr::left_join(
    test,
    getSam2pairwiseGBeven(sam2pairwise.GB),
    by = "coordinate.GB.even"
) %>% dplyr::rename(
        cigar.GB.even.2 = cigar.GB.even.y
    ) %>% dplyr::rename(
        cigar.GB.even.1 = cigar.GB.even.x
    )


#  Add coordinate information for 129, CAST, and mm10 -------------------------
test$coordinate.KA.129S1.odd <- paste0(
    stringr::str_split(test$coordinate.KA.odd, "_") %>%
        lapply(., `[[`, 1) %>%
        unlist(), "_",
    test$rname.129S1.odd, "_",
    test$pos.129S1.odd
)
test$coordinate.KA.129S1.even <- paste0(
    stringr::str_split(test$coordinate.KA.even, "_") %>%
        lapply(., `[[`, 1) %>%
        unlist(), "_",
    test$rname.129S1.even, "_",
    test$pos.129S1.even
)
test$coordinate.KA.CAST.odd <- paste0(
    stringr::str_split(test$coordinate.KA.odd, "_") %>%
        lapply(., `[[`, 1) %>%
        unlist(), "_",
    test$rname.CAST.odd, "_",
    test$pos.CAST.odd
)
test$coordinate.KA.CAST.even <- paste0(
    stringr::str_split(test$coordinate.KA.even, "_") %>%
        lapply(., `[[`, 1) %>%
        unlist(), "_",
    test$rname.CAST.even, "_",
    test$pos.CAST.even
)
test$coordinate.KA.mm10.odd <- paste0(
    stringr::str_split(test$coordinate.KA.odd, "_") %>%
        lapply(., `[[`, 1) %>%
        unlist(), "_",
    test$rname.mm10.odd, "_",
    test$pos.mm10.odd
)
test$coordinate.KA.mm10.even <- paste0(
    stringr::str_split(test$coordinate.KA.even, "_") %>%
        lapply(., `[[`, 1) %>%
        unlist(), "_",
    test$rname.mm10.even, "_",
    test$pos.mm10.even
)

test <- test %>%
    dplyr::relocate(
        c(
            coordinate.KA.129S1.odd,
            coordinate.KA.129S1.even,
            coordinate.KA.CAST.odd,
            coordinate.KA.CAST.even,
            coordinate.KA.mm10.odd,
            coordinate.KA.mm10.even
        ),
        .after = coordinate.KA.even
    )


###############################################################################
#  Join all of the sam2pairwise.KA tibbles ------------------------------------
tmp.odd.129S1 <- sam2pairwise.KA.129S1 %>%
    dplyr::select(
        qname, coordinate, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        coordinate.KA.129S1.odd = coordinate,
        cigar.KA.129S1.odd = cigar,
        read_sequence.KA.129S1.odd = read_sequence,
        matches.KA.129S1.odd = matches,
        reference_sequence.KA.129S1.odd = reference_sequence
    ) %>%
    dplyr::select(-qname)

tmp.odd.CAST <- sam2pairwise.KA.CAST %>%
    dplyr::select(
        qname, coordinate, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        coordinate.KA.CAST.odd = coordinate,
        cigar.KA.CAST.odd = cigar,
        read_sequence.KA.CAST.odd = read_sequence,
        matches.KA.CAST.odd = matches,
        reference_sequence.KA.CAST.odd = reference_sequence
    ) %>%
    dplyr::select(-qname)

tmp.odd.mm10 <- sam2pairwise.KA.mm10 %>%
    dplyr::select(
        qname, coordinate, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        coordinate.KA.mm10.odd = coordinate,
        cigar.KA.mm10.odd = cigar,
        read_sequence.KA.mm10.odd = read_sequence,
        matches.KA.mm10.odd = matches,
        reference_sequence.KA.mm10.odd = reference_sequence
    ) %>%
    dplyr::select(-qname)

tmp.even.129S1 <- tmp.odd.129S1
tmp.even.CAST <- tmp.odd.CAST
tmp.even.mm10 <- tmp.odd.mm10

colnames(tmp.even.129S1) <- gsub("odd", "even", colnames(tmp.even.129S1))
colnames(tmp.even.CAST) <- gsub("odd", "even", colnames(tmp.even.CAST))
colnames(tmp.even.mm10) <- gsub("odd", "even", colnames(tmp.even.mm10))


# -----------------------------------------------------------------------------
test.2 <- test
test.2 <- dplyr::left_join(
    test.2,
    tmp.odd.129S1,
    by = "coordinate.KA.129S1.odd"
) %>% 
    dplyr::left_join(
        .,
        tmp.even.129S1,
        by = "coordinate.KA.129S1.even"
    ) %>% 
    dplyr::left_join(
        .,
        tmp.odd.CAST,
        by = "coordinate.KA.CAST.odd"
    ) %>%
    dplyr::left_join(
        .,
        tmp.even.CAST,
        by = "coordinate.KA.CAST.even"
    ) %>%
    dplyr::left_join(
        .,
        tmp.odd.mm10,
        by = "coordinate.KA.mm10.odd"
    ) %>%
    dplyr::left_join(
        .,
        tmp.even.mm10,
        by = "coordinate.KA.mm10.even"
    )

(test.2$read_sequence.GB.odd == test.2$read_sequence.KA.mm10.odd) %>% table(., useNA = "ifany")
colnames(test.2)
test.odd.FALSE <- test.2[
    !(test.2$read_sequence.GB.odd == test.2$read_sequence.KA.mm10.odd) &
    !is.na(test.2$read_sequence.GB.odd == test.2$read_sequence.KA.mm10.odd), 
]
test.odd.FALSE$GB_by_KA %>% table(., useNA = "ifany")
test.odd.FALSE.trim <- test.odd.FALSE[, c(129:160)]
colnames(test.odd.FALSE.trim)
test.odd.FALSE.trim.trim <- test.odd.FALSE.trim[, c(1:8, 25:32)]
colnames(test.odd.FALSE.trim.trim)
test.odd.FALSE.trim.trim.trim <- test.odd.FALSE.trim.trim[, c(1, 9, 2, 10)]


# -----------------------------------------------------------------------------
tbl <- sam2pairwise.KA.CAST %>%
    dplyr::select(
        coordinate, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        coordinate.KA.CAST.odd = coordinate,
        cigar.KA.CAST.odd = cigar,
        read_sequence.KA.CAST.odd = read_sequence,
        matches.KA.CAST.odd = matches,
        reference_sequence.KA.CAST.odd = reference_sequence
    )
test.CAST <- test.CAST %>% dplyr::left_join(
    dplyr::distinct(
        tbl,
        coordinate.KA.CAST.odd,
        .keep_all = TRUE
    )
)

tbl <- sam2pairwise.KA.CAST %>%
    dplyr::select(
        coordinate, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        coordinate.KA.CAST.even = coordinate,
        cigar.KA.CAST.even = cigar,
        read_sequence.KA.CAST.even = read_sequence,
        matches.KA.CAST.even = matches,
        reference_sequence.KA.CAST.even = reference_sequence
    )
test.CAST <- test.CAST %>% dplyr::left_join(
    dplyr::distinct(
        tbl,
        coordinate.KA.CAST.even,
        .keep_all = TRUE
    )
)


# -----------------------------------------------------------------------------
tbl <- sam2pairwise.KA.mm10 %>%
    dplyr::select(
        coordinate, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        coordinate.KA.mm10.odd = coordinate,
        cigar.KA.mm10.odd = cigar,
        read_sequence.KA.mm10.odd = read_sequence,
        matches.KA.mm10.odd = matches,
        reference_sequence.KA.mm10.odd = reference_sequence
    )
test.Ambiguous <- test.Ambiguous %>% dplyr::left_join(
    dplyr::distinct(
        tbl,
        coordinate.KA.mm10.odd
    )
)

tbl <- sam2pairwise.KA.mm10 %>%
    dplyr::select(
        coordinate, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        coordinate.KA.mm10.even = coordinate,
        cigar.KA.mm10.even = cigar,
        read_sequence.KA.mm10.even = read_sequence,
        matches.KA.mm10.even = matches,
        reference_sequence.KA.mm10.even = reference_sequence
    )
test.Ambiguous <- test.Ambiguous %>% dplyr::left_join(
    dplyr::distinct(
        tbl,
        coordinate.KA.mm10.even
    )
)


# -----------------------------------------------------------------------------
tbl <- sam2pairwise.KA.mm10 %>%
    dplyr::select(
        coordinate, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        coordinate.KA.mm10.odd = coordinate,
        cigar.KA.mm10.odd = cigar,
        read_sequence.KA.mm10.odd = read_sequence,
        matches.KA.mm10.odd = matches,
        reference_sequence.KA.mm10.odd = reference_sequence
    )
test.NA <- test.NA %>% dplyr::left_join(
    dplyr::distinct(
        tbl,
        coordinate.KA.mm10.odd
    )
)

tbl <- sam2pairwise.KA.mm10 %>%
    dplyr::select(
        coordinate, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        coordinate.KA.mm10.even = coordinate,
        cigar.KA.mm10.even = cigar,
        read_sequence.KA.mm10.even = read_sequence,
        matches.KA.mm10.even = matches,
        reference_sequence.KA.mm10.even = reference_sequence
    )
test.NA <- test.NA %>% dplyr::left_join(
    dplyr::distinct(
        tbl,
        coordinate.KA.mm10.even
    )
)

rm(tbl)


# -----------------------------------------------------------------------------
test <- dplyr::bind_rows(
    test.129,
    test.CAST,
    test.Ambiguous,
    test.NA
)
colnames(test)

#  Remove unneeded variables
rm(
    bi.AS.pmin, bi.joint.GB, joint,
    test.129, test.Ambiguous, test.CAST, test.NA,
    sam2pairwise.GB, uni.AS.pmin, uni.joint.GB
)

operation <- paste0("rm(", variable, ")")
evaluateOperation(operation)


#  Divide main variable into list variable, listed by KA_by_GB categories -----
categories <- test %>%
    dplyr::group_by(KA_by_GB) %>%
    dplyr::group_split()


#  Save the image -------------------------------------------------------------
save.image(stringr::str_remove(script, ".R") %>% paste0(., ".Rdata"))


#  Script is completed; clean up environment ----------------------------------
rm(list = ls())

#  Go to 'test.PE-processing.part-8.R'
