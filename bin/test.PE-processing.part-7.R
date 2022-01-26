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
    
    #  Subsection: Written explanation
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
    
    #  Subsection: Written explanation
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
    operation <- paste(variable, command) #%>% paste(., collapse = "; ")
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


#  Script name ----------------------------------------------------------------
script <- "test.PE-processing.part-7.R"


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
#+ tibble number; save results to initial "tbl.*" variables
command <- paste0(
    "<- odd %>%",
        "dplyr::mutate(groupid = row_number()) %>% ",
        "dplyr::bind_rows(even %>% mutate(groupid = row_number())) %>% ",
        "dplyr::arrange(groupid, tibble) %>% ",
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


#  Create version of 'joint' with tibble for each "KA.* × GB.*" entry ---------
joint.assign <- joint %>%
    dplyr::group_by(KA_by_GB) %>%
    dplyr::group_split()

#  Clean out unneeded variables
# rm(bi.AS.pmin, bi.joint.GB, uni.AS.pmin, uni.joint.GB)


#  Cell-wise analyses of SNPs, etc. for [[5]] KA.129 × GB.NA ------------------
bam <- joint.assign[[5]] %>% dplyr::select(sort(dplyr::current_vars()))

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

#  Generate tibble of comparisons, metrics, etc. for [[5]] KA.129 × GB.NA
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

