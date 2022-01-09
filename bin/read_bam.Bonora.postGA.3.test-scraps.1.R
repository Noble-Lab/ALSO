#!/usr/bin/env Rscript

library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)


# Set up work directory (location TB∆) ----------------------------------------
directory_user <- "/Users/kalavattam"
directory_base <- "Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do"
directory_work <- "get_unique_fragments/Bonora"
setwd(
    paste0(directory_user, "/", directory_base, "/", directory_work)
)
rm(directory_base, directory_work)


#  Set up functions -----------------------------------------------------------
`%notin%` <- Negate(`%in%`)


evaluateOperation <- function(operation = operation) {
    return(eval(parse(text = operation), envir = .GlobalEnv))
}


makeOperation <- function(variable, command) {
    operation <- paste(variable, command) #%>% paste(., collapse = "; ")
    return(operation)
}


#  Assign variables necessary for subsequent commands -------------------------
# chromosome <- "chr1"
chromosome <- "chrX"
file <- list.files(pattern = paste0("\\", chromosome, ".bam$")) %>%
    stringr::str_subset("^dedup\\.", negate = TRUE) %>%
    stringr::str_subset("^mm10\\.", negate = TRUE)
variable <- file %>%
    stringr::str_split(., "-") %>%
    lapply(., `[[`, 1) %>%
    unlist() %>%
    paste0("dedup.", .)

#  Check to see what initial .bam files the variables correspond to
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)


#  Load the tibbles containing all liftOver assignments as .rds files ---------
command <- paste0(
    "<- ", "readRDS(",
        "file = \"", variable, ".", chromosome, ".post-liftOvers.rds\"",
    ")"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Work with only the "liftOver successful" observations ----------------------
l.variable <- paste0("l.", str_subset(variable, "mm10", negate = TRUE))
command <- paste0(
    "<- ", str_subset(variable, "mm10", negate = TRUE), " %>% ",
        "group_by(liftOver_reason) %>% ",
        "group_split()"
)
operation <- makeOperation(l.variable, command)
evaluateOperation(operation)


#  Create "minimal tibbles" ("m") of only ... and AS values -------------------
#+ 
#+ Rename the AS columns to identify their alignment of origin
s.variable <- c(
    paste0(l.variable, "[[4]]"),
    str_subset(variable, "mm10", negate = FALSE)
)
m.variable <- paste0("m.", variable)
suffix <- m.variable %>%
    strsplit(., "\\.") %>%
    lapply(., `[[`, 3) %>%
    unlist()
command <- paste0(
    "<- ", s.variable, " %>% ",
        "dplyr::select(coordinate, tag.AS) %>% ",
        "dplyr::rename(AS.", suffix, " = tag.AS)"
)
operation <- makeOperation(m.variable, command)
evaluateOperation(operation)

rm(l.dedup.129S1, l.dedup.CAST, l.variable, m.variable, s.variable)


#  Fully join the "m's" by ... ------------------------------------------------
#+ 
#+ e.g., see the following URL:
#+ https://stat545.com/join-cheatsheet.html#full_joinsuperheroes-publishers
m.full <- full_join(m.dedup.mm10, m.dedup.129S1, by = "coordinate") %>% 
    full_join(., m.dedup.CAST, by = "coordinate")


#  Calculate the differences between 129 and CAST alignment scores ------------
#+ 
#+ Initially assign categories to the difference based on values w/r/t/the
#+ variable int (see below)
int <- 0 %>% as.integer()  #TODO Make the integer an argument
m.full <- m.full %>% 
    mutate(difference = AS.129S1 - AS.CAST) %>% 
    dplyr::mutate(
        assignment.initial = case_when(
            difference >= (-1 * int) & difference <= int ~ "Neutral",
            difference > int ~ "129S1-SvImJ",
            difference < (-1 * int) ~ "CAST-EiJ"
        )
    )


#  Formal assignments of categories: See the block of code --------------------
m.full$assignment <- ifelse(
    is.na(m.full$assignment.initial),
    ifelse(
        !is.na(m.full$AS.129S1),
        "129S1-SvImJ",
        ifelse(
            !is.na(m.full$AS.CAST),
            "CAST-EiJ",
            m.full$assignment.initial
        )
    ),
    m.full$assignment.initial
) %>%
    forcats::as_factor()


#  Remove all but read, AS, and formal assignment columns ---------------------
m.full <- m.full %>%
    select(coordinate, AS.mm10, AS.129S1, AS.CAST, assignment)


#  Assign 0/1 categories for what aligned to what -----------------------------
m.full$tmp.mm10 <- ifelse(is.na(m.full$AS.mm10), "0", "1")
m.full$dedup.129S1 <- ifelse(is.na(m.full$AS.129S1), "0", "1")
m.full$dedup.CAST <- ifelse(is.na(m.full$AS.CAST), "0", "1")

m.full$trinary <- paste0(
    m.full$tmp.mm10,
    m.full$dedup.129S1,
    m.full$dedup.CAST
) %>%
    forcats::as_factor()

m.full$assignment_trinary <- paste(
    m.full$assignment,
    m.full$trinary
) %>%
    as_factor()


#  Remove all but read, AS, and formal assignment, and "trinary" columns ------
m.full <- m.full %>%
    select(
        coordinate,
        AS.mm10,
        AS.129S1,
        AS.CAST,
        assignment,
        trinary,
        assignment_trinary
    )


#  Create individual tibbles for each possible "assignment" value -------------
m.full$assignment <- m.full$assignment %>% as.character()
m.full$assignment[is.na(m.full$assignment)] <- "NA"
variable <- list(
    c("assignment.129S1", "129S1-SvImJ"),
    c("assignment.CAST", "CAST-EiJ"),
    c("assignment.Neutral", "Neutral"),
    c("assignment.NA", "NA")
)
command <- paste0(
    "<- m.full %>% ",
        "dplyr::filter_at(",
            "vars(assignment), ",
            "any_vars(. %in% \"", lapply(variable, `[[`, 2), "\")",
        ")"
)
operation <- makeOperation(lapply(variable, `[[`, 1), command)
evaluateOperation(operation)


#  Make $assignment factors again ---------------------------------------------
command <- paste0(
    "<- ", lapply(variable, `[[`, 1), "$assignment %>% ",
        "as.factor()"
)
operation <- makeOperation(
    paste0(lapply(variable, `[[`, 1), "$assignment"), command
)
evaluateOperation(operation)

assignment.NA$assignment <- NA
assignment.NA$assignment <- NA_character_

rm(variable, m.dedup.129S1, m.dedup.CAST, m.dedup.mm10, m.full)

bak.assignment.129S1 <- assignment.129S1
bak.assignment.CAST <- assignment.CAST
bak.assignment.Neutral <- assignment.Neutral
bak.assignment.NA <- assignment.NA

assignment.129S1 <- bak.assignment.129S1
assignment.CAST <- bak.assignment.CAST
assignment.Neutral <- bak.assignment.Neutral
assignment.NA <- bak.assignment.NA


#  Full joins for 129S1 and CAST ----------------------------------------------

#  Assignment.129S1
assignment.129S1 <- full_join(
    dedup.129S1, assignment.129S1, by = "coordinate"
) %>%
    tidyr::drop_na("AS.129S1") %>% 
    dplyr::select(-dplyr::one_of("AS.129S1")) %>%
    dplyr::rename(
        AS.129S1 = tag.AS,
        XS.129S1 = tag.XS,
        NM.129S1 = tag.NM
    ) %>%
    dplyr::select(-XS.129S1, -NM.129S1) %>%
    dplyr::relocate(AS.mm10, .after = seq) %>%
    dplyr::rename(  # new_name = old_name
        c(
            coordinate.129S1 = old_coordinate,
            qname.129S1 = qname,
            flag.129S1 = flag,
            rname.129S1 = old_rname,
            strand.129S1 = strand,
            pos.129S1 = old_pos,
            pos_end.129S1 = old_pos_end,
            rname.liftOver_129S1_mm10 = rname,
            pos.liftOver_129S1_mm10 = pos,
            pos_end.liftOver_129S1_mm10 = pos_end,
            mapq.129S1 = mapq,
            cigar.129S1 = cigar,
            mrnm.129S1 = mrnm,
            seq.129S1 = seq,
            isize.129S1 = isize,
            coordinate.129S1 = old_coordinate,
            rname.129S1 = old_rname,
            mpos.129S1 = old_mpos,
            mpos_end.129S1 = old_mpos_end
        )
    )

#  Assignment.CAST
assignment.CAST <- full_join(
    dedup.CAST, assignment.CAST, by = "coordinate"
) %>%
    tidyr::drop_na("AS.CAST") %>% 
    dplyr::select(-dplyr::one_of("AS.CAST")) %>%
    dplyr::rename(
        AS.CAST = tag.AS,
        XS.CAST = tag.XS,
        NM.CAST = tag.NM
    ) %>%
    dplyr::select(-XS.CAST, -NM.CAST) %>%
    dplyr::relocate(AS.CAST, .after = AS.129S1) %>%
    dplyr::rename(  # new_name = old_name
        c(
            coordinate.CAST = old_coordinate,
            qname.CAST = qname,
            flag.CAST = flag,
            rname.CAST = old_rname,
            strand.CAST = strand,
            pos.CAST = old_pos,
            pos_end.CAST = old_pos_end,
            rname.liftOver_CAST_mm10 = rname,
            pos.liftOver_CAST_mm10 = pos,
            pos_end.liftOver_CAST_mm10 = pos_end,
            mapq.CAST = mapq,
            cigar.CAST = cigar,
            mrnm.CAST = mrnm,
            seq.CAST = seq,
            isize.CAST = isize,
            coordinate.CAST = old_coordinate,
            rname.CAST = old_rname,
            mpos.CAST = old_mpos,
            mpos_end.CAST = old_mpos_end
        )
    )

#  Assignment.Neutral
assignment.Neutral <- full_join(
    dedup.129S1, assignment.Neutral, by = "coordinate"
) %>%
    tidyr::drop_na("AS.129S1") %>% 
    dplyr::select(-dplyr::one_of("AS.129S1")) %>%
    dplyr::rename(
        AS.129S1 = tag.AS,
        XS.129S1 = tag.XS,
        NM.129S1 = tag.NM
    ) %>%
    dplyr::select(-XS.129S1, -NM.129S1) %>%
    dplyr::relocate(AS.mm10, .before = AS.129S1) %>%
    dplyr::rename(  # new_name = old_name
        c(
            qname.129S1 = qname,
            flag.129S1 = flag,
            strand.129S1 = strand,
            rname.liftOver_129S1_mm10 = rname,
            pos.liftOver_129S1_mm10 = pos,
            pos_end.liftOver_129S1_mm10 = pos_end,
            mapq.129S1 = mapq,
            cigar.129S1 = cigar,
            mrnm.129S1 = mrnm,
            seq.129S1 = seq,
            isize.129S1 = isize,
            coordinate.129S1 = old_coordinate,
            rname.129S1 = old_rname,
            pos.129S1 = old_pos,
            pos_end.129S1 = old_pos_end,
            mpos.129S1 = old_mpos,
            mpos_end.129S1 = old_mpos_end
        )
    ) %>% full_join(
    .,
    dedup.CAST %>% dplyr::select(
        -liftOver_reason,
        -tag.XS,
        -tag.NM
    ),
    by = "coordinate"
    ) %>%
    tidyr::drop_na("AS.CAST") %>% 
    dplyr::select(-dplyr::one_of("AS.CAST")) %>%
    dplyr::rename(
        AS.CAST = tag.AS
    ) %>%
    dplyr::rename(  # new_name = old_name
        c(
            qname.CAST = qname,
            flag.CAST = flag,
            strand.CAST = strand,
            rname.liftOver_CAST_mm10 = rname,
            pos.liftOver_CAST_mm10 = pos,
            pos_end.liftOver_CAST_mm10 = pos_end,
            mapq.CAST = mapq,
            cigar.CAST = cigar,
            mrnm.CAST = mrnm,
            seq.CAST = seq,
            isize.CAST = isize,
            coordinate.CAST = old_coordinate,
            rname.CAST = old_rname,
            pos.CAST = old_pos,
            pos_end.CAST = old_pos_end,
            mpos.CAST = old_mpos,
            mpos_end.CAST = old_mpos_end
        )
    ) %>%
    dplyr::relocate(AS.CAST, .after = AS.129S1) %>%
    dplyr::relocate(coordinate.CAST, .after = coordinate.129S1) %>%
    dplyr::relocate(qname.CAST, .after = qname.129S1) %>%
    dplyr::relocate(flag.CAST, .after = flag.129S1) %>%
    dplyr::relocate(rname.CAST, .after = rname.129S1) %>%
    dplyr::relocate(strand.CAST, .after = strand.129S1) %>%
    dplyr::relocate(pos.CAST, .after = pos_end.129S1) %>%
    dplyr::relocate(pos_end.CAST, .after = pos.CAST) %>%
    dplyr::relocate(rname.liftOver_CAST_mm10, .after = pos_end.liftOver_129S1_mm10) %>%
    dplyr::relocate(pos.liftOver_CAST_mm10, .after = rname.liftOver_CAST_mm10) %>%
    dplyr::relocate(pos_end.liftOver_CAST_mm10, .after = pos.liftOver_CAST_mm10) %>%
    dplyr::relocate(mapq.CAST, .after = mapq.129S1) %>%
    dplyr::relocate(cigar.CAST, .after = cigar.129S1) %>%
    dplyr::relocate(mrnm.CAST, .after = mrnm.129S1) %>%
    dplyr::relocate(seq.CAST, .after = seq.129S1) %>%
    dplyr::relocate(isize.CAST, .after = isize.129S1) %>%
    dplyr::relocate(coordinate.CAST, .after = coordinate.129S1) %>%
    dplyr::relocate(rname.CAST, .after = rname.129S1) %>%
    dplyr::relocate(pos.CAST, .after = pos.129S1) %>%
    dplyr::relocate(pos_end.CAST, .after = pos_end.129S1) %>%
    dplyr::relocate(mpos.CAST, .after = mpos.129S1) %>%
    dplyr::relocate(mpos_end.CAST, .after = mpos_end.129S1)

#  Assignment.NA
assignment.NA <- full_join(dedup.mm10, assignment.NA, by = "coordinate") %>%
    tidyr::drop_na("AS.mm10") %>% 
    dplyr::select(-dplyr::one_of("AS.mm10")) %>%
    dplyr::rename(
        AS.mm10 = tag.AS,
        XS.mm10 = tag.XS,
        NM.mm10 = tag.NM
    ) %>%
    dplyr::select(-old_coordinate, -XS.mm10, -NM.mm10) %>%
    dplyr::rename(  # new_name = old_name
        c(
            qname.mm10 = qname,
            flag.mm10 = flag,
            rname.mm10 = rname,
            strand.mm10 = strand,
            pos.mm10 = pos,
            pos_end.mm10 = pos_end,
            mapq.mm10 = mapq,
            cigar.mm10 = cigar,
            mrnm.mm10 = mrnm,
            seq.mm10 = seq,
            isize.mm10 = isize,
            mpos.mm10 = mpos,
            mpos_end.mm10 = mpos_end
        )
    )


#  Save the "assignment tibbles" as RDS
variable <- paste0("assignment.", c("129S1", "CAST", "NA", "Neutral"))
operation <- paste0(
    "saveRDS(",
        variable, ", file = ", "\"", variable, ".", chromosome, ".rds\"",
    ")"
)
evaluateOperation(operation)
#TODO 1/2 Example of what these output files look like and/or some systematic
#TODO 2/2 way to tell they're from this script/step in the pipeline

#  Output is, e.g., "assignment.Neutral.chrX.rds"


#  Temporary, to be deleted ---------------------------------------------------
ne <- assignment.Neutral
which(ne$cigar.129S1 != ne$cigar.CAST) %>% length()  # 5
which(ne$mapq.129S1 != ne$mapq.CAST) %>% length()    # 10222
which(ne$seq.129S1 != ne$seq.CAST) %>% length()      # 22
which(ne$flag.129S1 != ne$flag.CAST) %>% length()    # 3

tmp <- ne[ne$seq.129S1 != ne$seq.CAST, ]
tmp <- ne[ne$mapq.129S1 != ne$mapq.CAST, ]  #TODO Possibility: Make consideration of MAPQ part of my assignment steps
tmp <- ne[ne$flag.129S1 != ne$flag.CAST, ]

rm(ne, tmp)

rm(list = ls())
