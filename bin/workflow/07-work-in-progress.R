
library("tidyverse")

# dir <- "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
dir <- "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/files_bam"
setwd(dir)

arguments <- list()
arguments$sample_1 <- "mm10"
arguments$sample_2 <- "CAST"
arguments$threshold <- 0

file.mm10 <- paste0(dir, "/Disteche_sample_1.dedup.mm10.corrected.AS.txt.gz")
file.CAST <- paste0(dir, "/Disteche_sample_1.dedup.CAST.corrected.AS.txt.gz")


#  Functions ------------------------------------------------------------------
readLine <- function(x, y, z) {
    #TODO Description of function
    #
    # :param x: txt/txt.gz file, including path (chr)
    # :param y: line number preceding line of interest (int)
    # :param z: identifier for sample (chr)
    # :return: #TODO
    a <- read.delim(
        x,
        header = FALSE,
        skip = y,
        nrows = 1
    ) %>%
        tibble::as_tibble() %>%
        dplyr::rename("qname" = "V1") %>%
        dplyr::rename(!!paste0("AS.", z) := "V2")
    return(a)
}


initializeTibble2 <- function(x) {
    #TODO Description of function
    #
    # :param x: identifier for sample (chr)
    # :return: #TODO
    y <- dplyr::tibble(
        qname = character(),
        !!paste0("AS.", x) := numeric()
    )
    return(y)
}


initializeTibble5 <- function(x, y) {
    #TODO Description of function
    #
    # :param x: identifier for sample #1 (chr)
    # :param x: identifier for sample #2 (chr)
    # :return: #TODO
    z <- dplyr::tibble(
        qname = character(),
        !!paste0("AS.", x) := numeric(),
        !!paste0("AS.", y) := numeric(),
        difference = numeric(),
        assignment = character()
    )
    return(z)
}


makeAssignmentViaAS <- function(x, y) {
    #TODO Description of function
    #
    # :param x: tibble comprised of vairables 'qname' (chr) and 'AS' (int)
    # :param y: threshold for assignments (int)
    x %>% 
        dplyr::mutate(difference = .[[2]] - .[[3]]) %>% 
        dplyr::mutate(
            assignment = case_when(
                difference >= (-1 * y) & difference <= y ~ "ambiguous",
                difference > y ~ gsub("AS\\.", "", colnames(x[2])),
                difference < (-1 * y) ~ gsub("AS\\.", "", colnames(x[3]))
            )
        )
}


#  Initialize tibbles for storing values
sample_1 <- sample_2 <- ambiguous <- initializeTibble5(
    arguments$sample_1,
    arguments$sample_2
)

sample_1_alone <- initializeTibble2(arguments$sample_1)
sample_2_alone <- initializeTibble2(arguments$sample_2)

#  Read in line of interest
i <- 2  #TODO Enclose the below in for loop
mm10 <- readLine(file.mm10, i, arguments$sample_1)
CAST <- readLine(file.CAST, i, arguments$sample_2)

#  If QNAMEs are the same, perform logic for AS comparisons; if not, then store
#+ values for additional logic
if(mm10[1] == CAST[1]) {
    AS <- dplyr::bind_cols(mm10, CAST[2])
    
    int <- arguments$threshold
    AS <- AS %>% 
        dplyr::mutate(difference = .[[2]] - .[[3]]) %>% 
        dplyr::mutate(
            assignment = case_when(
                difference >= (-1 * int) & difference <= int ~ "ambiguous",
                difference > int ~ gsub("AS\\.", "", colnames(AS[2])),
                difference < (-1 * int) ~ gsub("AS\\.", "", colnames(AS[3]))
            )
        )
    
    if(AS[5] == arguments$sample_1) {
        sample_1 <- dplyr::bind_rows(sample_1, AS)
    }
    
    if(AS[5] == arguments$sample_2) {
        sample_2 <- dplyr::bind_rows(sample_2, AS)
    }
    
    if(AS[5] == "ambiguous") {
        ambiguous <- dplyr::bind_rows(ambiguous, AS)
    }
} else {
    sample_1_alone <- dplyr::bind_rows(sample_1_alone, mm10)
    sample_2_alone <- dplyr::bind_rows(sample_2_alone, CAST)
    
    if(sample_1_alone$qname %in% sample_2_alone$qname) {
        cat("True")
        salvage_1 <- sample_1_alone[sample_1_alone$qname %in% sample_2_alone$qname, ]
        salvage_2 <- sample_2_alone[sample_2_alone$qname %in% sample_1_alone$qname, ]
    } else {

    }
}

# sample_1_alone <- dplyr::bind_rows(sample_1_alone, CAST)

line.mm10.u <- 
line.CAST.u <- read.delim(CAST, header = FALSE, skip = 0, nrows = 1) %>% tibble::as_tibble()


#  Notes on the logic
# do lexicographic sort of ♂.AS.txt.gz, ♀.AS.txt.gz
# 
# 	♀	♂
# a	1	1
# b	2	2
# c	3	4
# d	4	5
# e	5	6
# f	7	7
# g	9	10
# h	10	11
# ...
# 
# line a: match, ∴ do comparison and store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# line b: match, ∴ do comparison and store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# line c: mismatch, ∴ store 3 in tmp.♀, store 4 in tmp.♂; compare tmp.♀ with tmp.♂, find no matches; move on to next line
# line d: mismatch, ∴ store 4 in tmp.♀, store 5 in tmp.♂; compare tmp.♀ with tmp.♂, find that 4♀ and 4♂ match, store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# line e: mismatch, ∴ store 5 in tmp.♀, store 6 in tmp.♂; compare tmp.♀ with tmp.♂, find that 5♀ and 5♂ match, store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# line f: match, ∴ do comparison and store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# line g:	mismatch, ∴ store 9 in tmp.♀, store 10 in tmp.♂; compare tmp.♀ with tmp.♂, find no matches; move on to next line
# line h: mismatch, ∴ store 10 in tmp.♀, store 11 in tmp.♂; compare tmp.♀ with tmp.♂, find that 10♀ and 10♂ match, store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# 
# #  pseudocode
# read in line a of ♂.AS.txt.gz
# read in line a of ♀.AS.txt.gz
# 
# test if QNAME matches
# 	if so, do
# 		AS comparison
# 			if ♀.AS > ♂.AS, append QNAME to final.♀.txt.gz
# 			if ♀.AS < ♂.AS, append QNAME to final.♂.txt.gz
# 			if ♀.AS = ♂.AS, append QNAME to final.⚥.txt.gz
# 
# 	if not, do
# 		append QNAME.♀ to object tmp.♀ (in memory)
# 		append QNAME.♂ to object tmp.♂ (in memory)
# 
# 		test if QNAME in tmp.♀ matches QNAME in tmp.♂
# 			if so, do
# 				AS comparison
# 					if ♀.AS > ♂.AS, append QNAME to final.♀.txt.gz
# 					if ♀.AS < ♂.AS, append QNAME to final.♂.txt.gz
# 					if ♀.AS = ♂.AS, append QNAME to final.⚥.txt.gz
# 
# 			if not, do
# 				nothing
# 
# read in line b of ♂.AS.txt.gz
# read in line b of ♀.AS.txt.gz
# 
# test if QNAME matches
# ...
# 
# How to handle if lines remain in one file but all lines have been read in the other file?
# How to handle if, after comparisons between tmp.♀ and tmp.♂, there are no more matches?
