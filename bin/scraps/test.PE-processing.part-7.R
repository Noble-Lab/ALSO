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
    # rname.odd <- bam$lO_rname.129.odd
    # pos.odd <- bam$lO_pos.129.odd
    # pos_end.odd <- as.numeric(bam$lO_pos.129.odd + 49)
    # # pos_end.odd <- bam$lO_pos_end.129.odd
    # seq.odd <- bam$seq.129.odd
    # rname.even <- bam$lO_rname.129.even
    # pos.even <- bam$lO_pos.129.even
    # pos_end.even <- as.numeric(bam$lO_pos.129.even + 49)
    # # pos_end.even <- bam$lO_pos_end.129.even
    # seq.even <- bam$seq.129.even
    
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
directory_bam <- "results/kga0/2022-0214_segregated_reads.thresh_SNP_1.thresh_Q_30"
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
bi.AS.pmin <- a.full
uni.AS.pmin <- AS.pmin

bi.joint.GB <- joint.GB
uni.joint.GB <- uniline.joint.GB

rm(AS.pmin, a.full, joint.GB, uniline.joint.GB)

# ###############################################################################
# # CHECK START
# ###############################################################################
# colnames(uni.AS.pmin)
# 
# t <- uni.AS.pmin %>% select(c(pos.129.odd, pos.129.even))
# u <- uni.AS.pmin$pos.129.odd[is.na(uni.AS.pmin$pos.129.even)] %>% as_tibble_col()
# v <- uni.AS.pmin$pos.129.even[is.na(uni.AS.pmin$pos.129.odd)] %>% as_tibble_col()
# 
# w <- uni.AS.pmin %>% select(c(pos.CAST.odd, pos.CAST.even))
# x <- uni.AS.pmin$pos.CAST.odd[is.na(uni.AS.pmin$pos.CAST.even)] %>% as_tibble_col()
# y <- uni.AS.pmin$pos.CAST.even[is.na(uni.AS.pmin$pos.CAST.odd)] %>% as_tibble_col()
# 
# is.na(t$pos.129.odd) %>% table(useNA = "ifany")  # 111430 + 15313
# is.na(t$pos.129.even) %>% table(useNA = "ifany")  # 111430 + 15313
# is.na(w$pos.CAST.odd) %>% table(useNA = "ifany")  # 108755 + 17988
# is.na(w$pos.CAST.even) %>% table(useNA = "ifany")  # 108755 + 17988
# 
# t$pos.129.odd %in% t$pos.129.even %>% table(useNA = "ifany")  # 110707 + 16036
# w$pos.CAST.odd %in% w$pos.CAST.even %>% table(useNA = "ifany")  # 108050 + 18693
# 
# rm(t, u, v, w, x, y)
# ###############################################################################
# # CHECK END
# ###############################################################################


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
# colnames(df)

df <- df %>% dplyr::mutate(
    AS.mm10.odd = AS.mm10,
    AS.129.odd = AS.129,
    AS.CAST.odd = AS.CAST,
    assignment.odd = assignment,
    AS.mm10.even = AS.mm10,
    AS.129.even = AS.129,
    AS.CAST.even = AS.CAST,
    assignment.even = assignment
) %>%
    dplyr::select(-AS.mm10, -AS.129, -AS.CAST, -assignment) %>%
    dplyr::relocate(
        c(
            AS.mm10.odd, AS.129.odd, AS.CAST.odd, assignment.odd,
            AS.mm10.even, AS.129.even, AS.CAST.even, assignment.even
        ),
        .before = lO_criteria.odd
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
    "odd", " %>%",
        "dplyr::mutate(groupid = row_number())", " %>% ",
        "dplyr::bind_rows(even %>% mutate(groupid = row_number()))", " %>% ",
        "dplyr::arrange(groupid, tibble)", " %>% ",
        "dplyr::select(-tibble)"
)
operation <- makeOperation("bi.AS.pmin", command)
evaluateOperation(operation)

uni.AS.pmin <- df

# #  Again, check to make sure that there are no more than two entries per
# #+ group ID
# n_occur <- data.frame(table(eval(parse(text = paste0("AS.pmin")))$groupid))
# print(n_occur[n_occur$Freq > 2, ])

#  Clean up
rm(df, n_occur, even, odd)


###############################################################################
colnames(uni.joint.GB)[colnames(uni.joint.GB) == "lO_criteria.GB.odd"] <-
    "lO_criteria.odd"
colnames(uni.joint.GB)[colnames(uni.joint.GB) == "lO_criteria.GB.even"] <-
    "lO_criteria.even"

joint <- dplyr::full_join(
    uni.AS.pmin, uni.joint.GB, by = "lO_criteria.odd"
)
(joint$lO_criteria.odd %notin% uni.AS.pmin$lO_criteria.odd) %>%
    table(useNA = "ifany")  # FALSE: 126743
(joint$lO_criteria.odd %notin% uni.joint.GB$lO_criteria.odd) %>%
    table(useNA = "ifany")  # FALSE: 117341 TRUE: 9402

colnames(joint)[colnames(joint) == "lO_criteria.even.x"] <-
    "lO_criteria.even"

joint <- joint %>%
    dplyr::select(
        -lO_criteria.even.y,
        -coordinate.GB.odd, -coordinate.GB.even,
        -qpos.GB.odd, -qpos.GB.even,
        -qmpos.GB.odd, -qmpos.GB.even,
        -discrepancy.GB.odd, -discrepancy.GB.even,
        -groupid.GB.odd, -groupid.GB.even
    )

colnames(joint)[colnames(joint) == "assignment.odd"] <- "assignment.KA.odd"
colnames(joint)[colnames(joint) == "assignment.even"] <- "assignment.KA.even"

joint$assignment.KA.odd <- paste0("KA.", joint$assignment.KA.odd)
joint$assignment.GB.odd[is.na(joint$assignment.GB.odd)] <- "GB.NA"
joint$assignment.GB.even[is.na(joint$assignment.GB.even)] <- "NA"

joint$assignment.KA.even <- joint$assignment.KA.even %>% as.character()

order <- c("129", "CAST", "Ambiguous", "NA")
joint$assignment.KA.odd <- joint$assignment.KA.odd %>%
    stringr::str_remove("-EiJ") %>%
    stringr::str_remove("S1-SvImJ") %>%
    factor(., levels = paste0("KA.", order))
joint$assignment.KA.even <- joint$assignment.KA.even %>%
    stringr::str_remove("-EiJ") %>%
    stringr::str_remove("S1-SvImJ") %>%
    factor(., levels = order)


order <- c("129", "CAST", "Ambiguous", "Contra", "NA")
joint$assignment.GB.odd <- joint$assignment.GB.odd %>%
    stringr::str_remove("S1") %>%
    factor(., levels = paste0("GB.", order))
joint$assignment.GB.even <- joint$assignment.GB.even %>%
    stringr::str_remove("GB.") %>%
    stringr::str_remove("S1") %>%
    factor(., levels = order)

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
            rname.129.odd, pos.129.odd, pos_end.129.odd,
            rname.129.even, pos.129.even, pos_end.129.even,
            rname.CAST.odd, pos.CAST.odd, pos_end.CAST.odd,
            rname.CAST.even, pos.CAST.even, pos_end.CAST.even,
            lO_rname.129.odd, lO_pos.129.odd, lO_pos_end.129.odd,
            lO_rname.129.even, lO_pos.129.even, lO_pos_end.129.even,
            lO_rname.CAST.odd, lO_pos.CAST.odd, lO_pos_end.CAST.odd,
            lO_rname.CAST.even, lO_pos.CAST.even, lO_pos_end.CAST.even,
            mrnm.mm10.odd, mpos.mm10.odd, mpos_end.mm10.odd,
            mrnm.mm10.even, mpos.mm10.even, mpos_end.mm10.even,
            mrnm.129.odd, mpos.129.odd, mpos_end.129.odd,
            mrnm.129.even, mpos.129.even, mpos_end.129.even,
            mrnm.CAST.odd, mpos.CAST.odd, mpos_end.CAST.odd,
            mrnm.CAST.even, mpos.CAST.even, mpos_end.CAST.even,
            lO_mrnm.129.odd, lO_mpos.129.odd, lO_mpos_end.129.odd,
            lO_mrnm.129.even, lO_mpos.129.even, lO_mpos_end.129.even,
            lO_mrnm.CAST.odd, lO_mpos.CAST.odd, lO_mpos_end.CAST.odd,
            lO_mrnm.CAST.even, lO_mpos.CAST.even, lO_mpos_end.CAST.even
        ),
        .after = seq.CAST.even
    ) %>%
    dplyr::relocate(
        c(
            assignment.GB.odd, assignment.GB.even,
            flag.GB.odd, flag.GB.even,
            seq.GB.odd, seq.GB.even,
            rname.GB.odd, pos.GB.odd, pos_end.GB.odd,
            rname.GB.even, pos.GB.even, pos_end.GB.even,
            mrnm.GB.odd, mpos.GB.odd, mpos_end.GB.odd,
            mrnm.GB.even, mpos.GB.even, mpos_end.GB.even,
            qname.GB.odd, qname.GB.even,
            mapq.GB.odd, mapq.GB.even,
            cigar.GB.odd, cigar.GB.even,
            isize.GB.odd, isize_abs.GB.odd,
            isize.GB.even, isize_abs.GB.even,
            qual.GB.odd, qual.GB.even
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
getwd()

#  Load in sam2pairwise files for .bam files output by GB's allele assignment
#+ script
file <- Sys.glob(
    paste0("*.", chromosome, "*.extendedCIGAR.sam2pairwise.munged.txt")
)
variable <- file %>%
    strsplit(., "\\.") %>%
    lapply(., `[[`, 5) %>% unlist() %>%
    paste0("sam2pairwise.GB.", chromosome, ".", .)
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

command <- paste0("readr::read_tsv(", variable, ")")
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  Save the sam2pairwise variables for GB .bam files
variable.tmp <- variable

#  Switch to location of sam2pairwise files for KA .bam files
setwd(paste0(path.1, "/", path.3))

#  sam2pairwise files for .bam files output by KA's allele assignment script
file <- Sys.glob(
    paste0(
        "processed.*", chromosome, "*.extendedCIGAR.sam2pairwise.munged.txt"
    )
)
variable <- file %>%
    strsplit(., "\\.") %>%
    lapply(., `[[`, 3) %>% unlist() %>%
    paste0("sam2pairwise.KA.", chromosome, ".", .) %>%
    stringr::str_remove("S1")
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

command <- paste0("readr::read_tsv(", variable, ")")
operation <- makeOperation(variable, command)
evaluateOperation(operation)

variable.tmp.2 <- variable

#  Create a variable for the sam2pairwise variables for GB and KA .bam files
variable <- c(variable.tmp, variable.tmp.2)
rm(variable.tmp, variable.tmp.2)


#  Create $lO_criteria columns for subsequent join operations -----------------
command <- paste0(
    variable, " %>% ",
        "tidyr::unite(",
            "lO_criteria, ",
            "c(\"qname\", \"flag\", \"pos\", \"mpos\"), ",
            "sep = \"_\", ",
            "remove = FALSE",
        ")"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Combine sam2pairwise.GB tibbles into one tibble ----------------------------
command <- paste0(
    "dplyr::bind_rows(",
        stringr::str_subset(variable, "KA", negate = TRUE)[1], ", ",
        stringr::str_subset(variable, "KA", negate = TRUE)[2], ", ",
        stringr::str_subset(variable, "KA", negate = TRUE)[3], ", ",
        stringr::str_subset(variable, "KA", negate = TRUE)[4],
    ")"
)
operation <- makeOperation(paste0("sam2pairwise.GB.", chromosome), command)
evaluateOperation(operation)

#  Clean out no-longer-needed variables
operation <- paste0(
    "rm(",
        stringr::str_subset(variable, "KA", negate = TRUE)[1], ", ",
        stringr::str_subset(variable, "KA", negate = TRUE)[2], ", ",
        stringr::str_subset(variable, "KA", negate = TRUE)[3], ", ",
        stringr::str_subset(variable, "KA", negate = TRUE)[4],
    ")"
)
evaluateOperation(operation)


#  Get sam2pairwise.GB tibble mate pairs into one row each --------------------
getSam2pairwiseGBodd <- function(tibble) {
    tbl <- tibble %>%
    dplyr::select(
        lO_criteria, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        lO_criteria.odd = lO_criteria,
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
        lO_criteria, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        lO_criteria.even = lO_criteria,
        cigar.GB.even = cigar,
        read_sequence.GB.even = read_sequence,
        matches.GB.even = matches,
        reference_sequence.GB.even = reference_sequence
    )
    return(tbl)
}


joint <- dplyr::left_join(
    joint,
    getSam2pairwiseGBodd(eval(parse(
        text = paste0("sam2pairwise.GB.", chromosome)
    ))),
    by = "lO_criteria.odd"
) %>% dplyr::rename(
        cigar.GB.odd.2 = cigar.GB.odd.y
    ) %>% dplyr::rename(
        cigar.GB.odd.1 = cigar.GB.odd.x
    )

joint <- dplyr::left_join(
    joint,
    getSam2pairwiseGBeven(eval(parse(
        text = paste0("sam2pairwise.GB.", chromosome)
    ))),
    by = "lO_criteria.even"
) %>% dplyr::rename(
        cigar.GB.even.2 = cigar.GB.even.y
    ) %>% dplyr::rename(
        cigar.GB.even.1 = cigar.GB.even.x
    )


# #  Add coordinate information for 129, CAST, and mm10 -------------------------
# colnames(joint) <- colnames(joint) %>% stringr::str_remove("S1")
# 
# joint$coordinate.KA.129.odd <- paste0(
#     stringr::str_split(joint$coordinate.KA.odd, "_") %>%
#         lapply(., `[[`, 1) %>%
#         unlist(), "_",
#     joint$rname.129.odd, "_",
#     joint$pos.129.odd
# )
# joint$coordinate.KA.129.even <- paste0(
#     stringr::str_split(joint$coordinate.KA.even, "_") %>%
#         lapply(., `[[`, 1) %>%
#         unlist(), "_",
#     joint$rname.129.even, "_",
#     joint$pos.129.even
# )
# joint$coordinate.KA.CAST.odd <- paste0(
#     stringr::str_split(joint$coordinate.KA.odd, "_") %>%
#         lapply(., `[[`, 1) %>%
#         unlist(), "_",
#     joint$rname.CAST.odd, "_",
#     joint$pos.CAST.odd
# )
# joint$coordinate.KA.CAST.even <- paste0(
#     stringr::str_split(joint$coordinate.KA.even, "_") %>%
#         lapply(., `[[`, 1) %>%
#         unlist(), "_",
#     joint$rname.CAST.even, "_",
#     joint$pos.CAST.even
# )
# joint$coordinate.KA.mm10.odd <- paste0(
#     stringr::str_split(joint$coordinate.KA.odd, "_") %>%
#         lapply(., `[[`, 1) %>%
#         unlist(), "_",
#     joint$rname.mm10.odd, "_",
#     joint$pos.mm10.odd
# )
# joint$coordinate.KA.mm10.even <- paste0(
#     stringr::str_split(joint$coordinate.KA.even, "_") %>%
#         lapply(., `[[`, 1) %>%
#         unlist(), "_",
#     joint$rname.mm10.even, "_",
#     joint$pos.mm10.even
# )
# # test <- joint[, c(131:136)]
# # rm(test)
# 
# joint <- joint %>%
#     dplyr::relocate(
#         c(
#             coordinate.KA.129.odd,
#             coordinate.KA.129.even,
#             coordinate.KA.CAST.odd,
#             coordinate.KA.CAST.even,
#             coordinate.KA.mm10.odd,
#             coordinate.KA.mm10.even
#         ),
#         .after = coordinate.KA.even
#     )


#  Join all of the sam2pairwise.KA tibbles ------------------------------------
suffix <- c("129", "CAST", "mm10")
tmp.odd.129 <- eval(parse(
    text = paste0("sam2pairwise.KA.", chromosome, ".", suffix[1])
)) %>%
    dplyr::select(
        qname, lO_criteria, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        lO_criteria.odd = lO_criteria,
        cigar.KA.129.odd = cigar,
        read_sequence.KA.129.odd = read_sequence,
        matches.KA.129.odd = matches,
        reference_sequence.KA.129.odd = reference_sequence
    ) %>%
    dplyr::select(-qname)

tmp.odd.CAST <- eval(parse(
    text = paste0("sam2pairwise.KA.", chromosome, ".", suffix[2])
)) %>%
    dplyr::select(
        qname, lO_criteria, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        lO_criteria.odd = lO_criteria,
        cigar.KA.CAST.odd = cigar,
        read_sequence.KA.CAST.odd = read_sequence,
        matches.KA.CAST.odd = matches,
        reference_sequence.KA.CAST.odd = reference_sequence
    ) %>%
    dplyr::select(-qname)

tmp.odd.mm10 <- eval(parse(
    text = paste0("sam2pairwise.KA.", chromosome, ".", suffix[3])
)) %>%
    dplyr::select(
        qname, lO_criteria, cigar, read_sequence, matches, reference_sequence
    ) %>%
    dplyr::rename(
        lO_criteria.odd = lO_criteria,
        cigar.KA.mm10.odd = cigar,
        read_sequence.KA.mm10.odd = read_sequence,
        matches.KA.mm10.odd = matches,
        reference_sequence.KA.mm10.odd = reference_sequence
    ) %>%
    dplyr::select(-qname)

tmp.even.129 <- tmp.odd.129
tmp.even.CAST <- tmp.odd.CAST
tmp.even.mm10 <- tmp.odd.mm10

colnames(tmp.even.129) <- gsub("odd", "even", colnames(tmp.even.129))
colnames(tmp.even.CAST) <- gsub("odd", "even", colnames(tmp.even.CAST))
colnames(tmp.even.mm10) <- gsub("odd", "even", colnames(tmp.even.mm10))


# -----------------------------------------------------------------------------
joint <- joint %>%
    dplyr::left_join(
        .,
        tmp.odd.129,
        by = "lO_criteria.odd"
    ) %>% 
    dplyr::left_join(
        .,
        tmp.even.129,
        by = "lO_criteria.even"
    ) %>% 
    dplyr::left_join(
        .,
        tmp.odd.CAST,
        by = "lO_criteria.odd"
    ) %>%
    dplyr::left_join(
        .,
        tmp.even.CAST,
        by = "lO_criteria.even"
    ) %>%
    dplyr::left_join(
        .,
        tmp.odd.mm10,
        by = "lO_criteria.odd"
    ) %>%
    dplyr::left_join(
        .,
        tmp.even.mm10,
        by = "lO_criteria.even"
    )
rm(list = ls(pattern = "tmp."))


# # -----------------------------------------------------------------------------
# tbl <- sam2pairwise.KA.CAST %>%
#     dplyr::select(
#         coordinate, cigar, read_sequence, matches, reference_sequence
#     ) %>%
#     dplyr::rename(
#         coordinate.KA.CAST.odd = coordinate,
#         cigar.KA.CAST.odd = cigar,
#         read_sequence.KA.CAST.odd = read_sequence,
#         matches.KA.CAST.odd = matches,
#         reference_sequence.KA.CAST.odd = reference_sequence
#     )
# test.CAST <- test.CAST %>% dplyr::left_join(
#     dplyr::distinct(
#         tbl,
#         coordinate.KA.CAST.odd,
#         .keep_all = TRUE
#     )
# )
# 
# tbl <- sam2pairwise.KA.CAST %>%
#     dplyr::select(
#         coordinate, cigar, read_sequence, matches, reference_sequence
#     ) %>%
#     dplyr::rename(
#         coordinate.KA.CAST.even = coordinate,
#         cigar.KA.CAST.even = cigar,
#         read_sequence.KA.CAST.even = read_sequence,
#         matches.KA.CAST.even = matches,
#         reference_sequence.KA.CAST.even = reference_sequence
#     )
# test.CAST <- test.CAST %>% dplyr::left_join(
#     dplyr::distinct(
#         tbl,
#         coordinate.KA.CAST.even,
#         .keep_all = TRUE
#     )
# )
# 
# 
# # -----------------------------------------------------------------------------
# tbl <- sam2pairwise.KA.mm10 %>%
#     dplyr::select(
#         coordinate, cigar, read_sequence, matches, reference_sequence
#     ) %>%
#     dplyr::rename(
#         coordinate.KA.mm10.odd = coordinate,
#         cigar.KA.mm10.odd = cigar,
#         read_sequence.KA.mm10.odd = read_sequence,
#         matches.KA.mm10.odd = matches,
#         reference_sequence.KA.mm10.odd = reference_sequence
#     )
# test.Ambiguous <- test.Ambiguous %>% dplyr::left_join(
#     dplyr::distinct(
#         tbl,
#         coordinate.KA.mm10.odd
#     )
# )
# 
# tbl <- sam2pairwise.KA.mm10 %>%
#     dplyr::select(
#         coordinate, cigar, read_sequence, matches, reference_sequence
#     ) %>%
#     dplyr::rename(
#         coordinate.KA.mm10.even = coordinate,
#         cigar.KA.mm10.even = cigar,
#         read_sequence.KA.mm10.even = read_sequence,
#         matches.KA.mm10.even = matches,
#         reference_sequence.KA.mm10.even = reference_sequence
#     )
# test.Ambiguous <- test.Ambiguous %>% dplyr::left_join(
#     dplyr::distinct(
#         tbl,
#         coordinate.KA.mm10.even
#     )
# )
# 
# 
# # -----------------------------------------------------------------------------
# tbl <- sam2pairwise.KA.mm10 %>%
#     dplyr::select(
#         coordinate, cigar, read_sequence, matches, reference_sequence
#     ) %>%
#     dplyr::rename(
#         coordinate.KA.mm10.odd = coordinate,
#         cigar.KA.mm10.odd = cigar,
#         read_sequence.KA.mm10.odd = read_sequence,
#         matches.KA.mm10.odd = matches,
#         reference_sequence.KA.mm10.odd = reference_sequence
#     )
# test.NA <- test.NA %>% dplyr::left_join(
#     dplyr::distinct(
#         tbl,
#         coordinate.KA.mm10.odd
#     )
# )
# 
# tbl <- sam2pairwise.KA.mm10 %>%
#     dplyr::select(
#         coordinate, cigar, read_sequence, matches, reference_sequence
#     ) %>%
#     dplyr::rename(
#         coordinate.KA.mm10.even = coordinate,
#         cigar.KA.mm10.even = cigar,
#         read_sequence.KA.mm10.even = read_sequence,
#         matches.KA.mm10.even = matches,
#         reference_sequence.KA.mm10.even = reference_sequence
#     )
# test.NA <- test.NA %>% dplyr::left_join(
#     dplyr::distinct(
#         tbl,
#         coordinate.KA.mm10.even
#     )
# )
# 
# rm(tbl)
# 
# 
# # -----------------------------------------------------------------------------
# test <- dplyr::bind_rows(
#     test.129,
#     test.CAST,
#     test.Ambiguous,
#     test.NA
# )
# colnames(test)

#  Remove unneeded variables
rm(list = ls(pattern = "bi."))
rm(list = ls(pattern = "sam2pairwise."))
rm(list = ls(pattern = "uni."))
rm(order, dedup.129S1)


#  Divide main variable into list variable, listed by KA_by_GB categories -----
categories <- joint %>%
    dplyr::group_by(KA_by_GB) %>%
    dplyr::group_split()


#  Save the image -------------------------------------------------------------
save.image(stringr::str_remove(script, ".R") %>% paste0(., ".Rdata"))


#  Script is completed; clean up environment ----------------------------------
rm(list = ls())

#  Go to 'test.PE-processing.part-8.R'
