
library("tidyverse")


#  Functions ------------------------------------------------------------------
countRecords <- function(x) {
    # Count number of records in txt.gz file
    # 
    # :param x: txt.gz file, including path (chr)
    # :return y: number of records in txt.gz file (int)
    y <- system(paste0("zcat ", x, " | wc -l"), intern = TRUE) %>%
        as.integer()
    return(y)
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


makeAssignment <- function(x, y) {
    #TODO Description of function
    #
    # :param x: tibble comprised of the following variables: qname (chr), AS
    #           for sample 1 (int), AS for sample 2 (int)
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


writeList <- function(x, y, z, a) {
    #TODO Description of function
    #
    # :param x: tibble
    # :param y: outdirectory, including path (chr)
    # :param z: outfile prefix (chr)
    # :param a: sample string (chr)
    readr::write_tsv(
        dplyr::bind_rows(x),
        paste0(y, "/", z, ".assign-", a, ".txt.gz"),
        append = TRUE
    )
}


#  Location, variables --------------------------------------------------------
# dir <- "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
dir <- "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/files_bam"
setwd(dir)

arguments <- list()
arguments$sample_1 <- paste0(dir, "/Disteche_sample_1.dedup.mm10.corrected.AS.2600.txt.gz")
arguments$sample_2 <- paste0(dir, "/Disteche_sample_1.dedup.CAST.corrected.AS.2500.txt.gz")
arguments$string_1 <- "mm10"
arguments$string_2 <- "CAST"
arguments$threshold <- 0
arguments$outdir <- dir
arguments$outprefix <- unlist(stringr::str_split(basename(arguments$sample_1), "\\."))[1]

count.sample_1 <- countRecords(arguments$sample_1)
count.sample_2 <- countRecords(arguments$sample_2)

int <- arguments$threshold


#  Initialize tibbles for storing values
sample_1 <- sample_2 <- ambiguous <- AS <- salvage_AS <- initializeTibble5(
    arguments$string_1,
    arguments$string_2
)

mm10 <- alone_1 <- salvage_1 <- initializeTibble2(arguments$string_1)
CAST <- alone_2 <- salvage_2 <- initializeTibble2(arguments$string_2)


#  Process lines of interest --------------------------------------------------
n <- pmax(count.sample_1, count.sample_2)

# test.mm10 <- readLine(arguments$sample_1, n - 1, arguments$string_1)
# n <- 1000L
bar <- utils::txtProgressBar(min = 0, max = n, initial = 0, style = 3)
for(i in 1:n) {
    #  Read in line of txt.gz file only if that line number is present in the
    #+ txt.gz file; i.e., don't read in a line number greater than the total
    #+ number of lines in the file
    if(i < count.sample_1) {
        mm10 <- readLine(arguments$sample_1, i, arguments$string_1)
    }
    if(i < count.sample_2) {
        CAST <- readLine(arguments$sample_2, i, arguments$string_2)
    }
    
    #  Test: Do QNAMEs match on a per-line basis? If they match, then perform
    #+ logic for AS comparisons (see function makeAssignment()); if not, then
    #+ store values in "alone_" variables (further description below)
    if(mm10[1] == CAST[1]) {
        AS <- dplyr::bind_cols(mm10, CAST[2]) %>% makeAssignment(., int)
        
        if(AS[5] == arguments$string_1) {
            sample_1 <- dplyr::bind_rows(sample_1, AS)
            writeList(
                AS,
                arguments$outdir,
                arguments$outprefix,
                arguments$string_1
            )
        } else if(AS[5] == arguments$string_2) {
            sample_2 <- dplyr::bind_rows(sample_2, AS)
            writeList(
                AS,
                arguments$outdir,
                arguments$outprefix,
                arguments$string_2
            )
        } else if(AS[5] == "ambiguous") {
            ambiguous <- dplyr::bind_rows(ambiguous, AS)
            writeList(
                AS,
                arguments$outdir,
                arguments$outprefix,
                "ambiguous"
            )
        }
    } else {
        if(i < count.sample_1) alone_1 <- dplyr::bind_rows(alone_1, mm10)
        if(i < count.sample_2) alone_2 <- dplyr::bind_rows(alone_2, CAST)
        
        #  Attempt to "salvage" unmatched entries
        if(sum(alone_1$qname %in% alone_2$qname) == 1) {
            salvage_1 <- alone_1[alone_1$qname %in% alone_2$qname, ]
            salvage_2 <- alone_2[alone_2$qname %in% alone_1$qname, ]
            salvage_AS <- dplyr::bind_cols(salvage_1, salvage_2[2]) %>%
                makeAssignment(., int)
            
            if(salvage_AS[5] == arguments$string_1) {
                sample_1 <- dplyr::bind_rows(sample_1, salvage_AS)
                writeList(
                    salvage_AS,
                    arguments$outdir,
                    arguments$outprefix,
                    arguments$string_1
                )
            } else if(salvage_AS[5] == arguments$string_2) {
                sample_2 <- dplyr::bind_rows(sample_2, salvage_AS)
                writeList(
                    salvage_AS,
                    arguments$outdir,
                    arguments$outprefix,
                    arguments$string_2
                )
            } else if(salvage_AS[5] == "ambiguous") {
                ambiguous <- dplyr::bind_rows(ambiguous, salvage_AS)
                writeList(
                    salvage_AS,
                    arguments$outdir,
                    arguments$outprefix,
                    "ambiguous"
                )
            }
            
            #  To prevent duplicate searches, and since they're no longer
            #+ "alone", remove QNAME from alone_1 and alone_2: The
            #+ dplyr::filter approach is finding all rows not equal to the one
            #+ to exclude, then rewriting those to the tibble; I think this
            #+ will slow down the script as alone_1 and alone_2 increase in
            #+ size
            alone_1 <- alone_1 %>%
                dplyr::filter(qname != as.character(salvage_1[1]))
            alone_2 <- alone_2 %>%
                dplyr::filter(qname != as.character(salvage_2[1]))
        }
        
        #  Only write out the lists of "alone" QNAMEs (i.e., in one sample but
        #+ not other) after having read in and run logic on all lines; #TODO
        #+ May need to consider a strategy in which we write x number of
        #+ entries and remove them from memory after reaching a certain number
        #+ of lines in the object or after reaching a certain size in RAM;
        #+ e.g., if the "alone" object reaches 1GB in size, then write the
        #+ first x number of lines to disc and remove them from memory (I think
        #+ this approach should be find because of the lexical-sorted nature
        #+ of the input: earlier entries in one sample are very likely to be
        #+ with earlier entries of other sample and very unlikely to be with
        #+ later entries in that sample)
        if(i == n) {
            writeList(
                alone_1,
                arguments$outdir,
                arguments$outprefix,
                arguments$string_1
            )
            writeList(
                alone_2,
                arguments$outdir,
                arguments$outprefix,
                arguments$string_2
            )
        }
    }
    utils::setTxtProgressBar(bar, i)
}
mm10
pryr::object_size(mm10)
CAST
pryr::object_size(CAST)
AS
pryr::object_size(AS)

alone_1
pryr::object_size(alone_1)
alone_2
pryr::object_size(alone_2)

ambiguous
pryr::object_size(ambiguous)
sample_1
pryr::object_size(sample_1)
sample_2
pryr::object_size(sample_2)

i
pryr::object_size(i)


# ###########################################################################
# #  Just for testing: Add an extra value
# test <- CAST %>% dplyr::rename(AS.mm10 = AS.CAST)
# test[, 2] <- -2
# alone_1 <- dplyr::bind_rows(alone_1, test)
# ###########################################################################

#  Notes on the logic
# do lexicographic sort of ♂.AS.txt.gz, ♀.AS.txt.gz
# 
# 	♀	♂
# 1	1	1
# 2	2	2
# 3	3	4
# 4	4	5
# 5	5	6
# 6	7	7
# 7	9	10
# 8	10	11
# ...
# 
# line 1: match, ∴ do comparison and store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# line 2: match, ∴ do comparison and store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# line 3: mismatch, ∴ store 3 in tmp.♀, store 4 in tmp.♂; compare tmp.♀ with tmp.♂, find no matches; move on to next line
# line 4: mismatch, ∴ store 4 in tmp.♀, store 5 in tmp.♂; compare tmp.♀ with tmp.♂, find that 4♀ and 4♂ match, store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# line 5: mismatch, ∴ store 5 in tmp.♀, store 6 in tmp.♂; compare tmp.♀ with tmp.♂, find that 5♀ and 5♂ match, store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# line 6: match, ∴ do comparison and store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# line 7:	mismatch, ∴ store 9 in tmp.♀, store 10 in tmp.♂; compare tmp.♀ with tmp.♂, find no matches; move on to next line
# line 8: mismatch, ∴ store 10 in tmp.♀, store 11 in tmp.♂; compare tmp.♀ with tmp.♂, find that 10♀ and 10♂ match, store results (final.♀.txt.gz, final.♂.txt.gz, or final.⚥.txt.gz); move on to next line
# 
# 



#  pseudocode takes in two files as input

#  file 1 is, e.g., ♂.AS.txt.gz
#+ (header and first line of ♂.AS.txt.gz)
#+ qname   AS
#+ AACCAATAACAAGTCAACGGAATATAGCGGCAAGTCAGCA:1417636611   -5

#  file 2 is, e.g., ♀.AS.txt.gz
#+ (header and first line of ♀.AS.txt.gz)
#+ qname   AS
#+ AACCAATAACAAGTCAACGGAATATAGCGGCAAGTCAGCA:1417636611   0

#  find the max number of records between ♂.AS.txt.gz and ♀.AS.txt.gz
# n = maximum of ♂.AS.txt.gz and ♀.AS.txt.gz

#  loop over n
# for i in n, do
#     read in line i of ♂.AS.txt.gz
#     read in line i of ♀.AS.txt.gz
#
#     test if ♀.AS QNAME matches ♂.AS QNAME
#     	if so, do
#     		AS comparison
#     			if ♀.AS > ♂.AS, append QNAME to final.♀.txt.gz (write to disc)
#     			if ♀.AS < ♂.AS, append QNAME to final.♂.txt.gz (write to disc)
#     			if ♀.AS = ♂.AS, append QNAME to final.⚥.txt.gz (write to disc)
#
#     	if not, do
#     		append QNAME.♀ to object tmp.♀ (store in memory)
#     		append QNAME.♂ to object tmp.♂ (store in memory)
#
#     		test if QNAME in tmp.♀ matches QNAME in tmp.♂
#     			if so, do
#     				AS comparison
#     					if ♀.AS > ♂.AS,
#                            append QNAME to final.♀.txt.gz (write to disc),
#                            then
#                                remove QNAME from both tmp.♀ and tmp.♂ (remove from memory)
#     					if ♀.AS < ♂.AS,
#                            append QNAME to final.♂.txt.gz (write to disc),
#                            then
#                                remove QNAME from both tmp.♀ and tmp.♂  (remove from memory)
#     					if ♀.AS = ♂.AS,
#                            append QNAME to final.⚥.txt.gz (write to disc),
#                            then
#                                remove QNAME from both tmp.♀ and tmp.♂  (remove from memory)
#
#     			if not, do
#     				nothing

# How to handle if lines remain in one file but all lines have been read in the other file?  #DONE, see code
# How to handle if, after comparisons between tmp.♀ and tmp.♂, there are no more matches?  #DONE, see code
