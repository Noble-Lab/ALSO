#!/usr/bin/env Rscript

# library(BSgenome)
# library(Mmusculus129S1SvImJ)
# library(MmusculusCASTEiJ)
# library(Mmusculus129Inserted)
# library(MmusculusCASTInserted)
# library(MmusculusCAST129Inserted)
# library(Mmusculus129Nmasked)
# library(MmusculusCASTNmasked)
# library(MmusculusCAST129Nmasked)

library(ggplot2)
library(scales)
library(stringr)
library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)


#  Functions ------------------------------------------------------------------
mungeMatesIntoOneRow <- function(variable_in, sort_by) {
    #  Input should be of class "character"
    if (class(variable_in) != "character") {
        stop("Argument 'variable_in' is not class 'character'. Exiting.")
    } else {
        variable_out <- eval(parse(text = variable_in))
    }
    
    if (class(sort_by) != "character") {
        stop("Argument 'sort_by' should be class 'character'. Exiting.")
    }

    odd.seq <- seq(1, nrow(variable_out), 2)
    even.seq <- (seq(1, nrow(variable_out), 2)) + 1
    
    odd <- variable_out[odd.seq, ]
    even <- variable_out[even.seq, ]
    
    odd <- odd %>%
        dplyr::mutate(groupid = row_number()) #%>%
        #dplyr::relocate(groupid, .after = assignment)
    even <- even %>%
        dplyr::mutate(groupid = row_number()) #%>%
        #dplyr::relocate(groupid, .after = assignment)
    
    variable_out <- full_join(odd, even, by = "groupid") %>%
        dplyr::arrange(sort_by) %>%  # Sort by pos.x
        dplyr::rename(groupid.x = groupid)  # new_name = old_name
    variable_out$groupid.y <- variable_out$groupid.x
    variable_out <- variable_out #%>%
        #dplyr::relocate(groupid.y, .after = assignment.y)
    
    colnames(variable_out) <- str_replace_all(
        colnames(variable_out), "\\.x", "\\.odd"
    )
    colnames(variable_out) <- str_replace_all(
        colnames(variable_out), "\\.y", "\\.even"
    )
    
    return(variable_out)
}


mungeMatesIntoTwoRows <- function(variable_in, sort_by) {
    #  All function arguments should be of class "character"
    if (class(variable_in) != "character") {
        stop("Argument 'variable_in' is not class 'character'. Exiting.")
    } else {
        variable_out <- eval(parse(text = variable_in))
    }
    
    if (class(sort_by) != "character") {
        stop("Argument 'sort_by' should be class 'character'. Exiting.")
    }
    
    #  Sort the uniline tibbles by...
    variable_out <- variable_out %>% dplyr::arrange(sort_by)

    #  Split the uniline tibbles by odd or even status
    odd <- variable_out[
        stringr::str_subset(colnames(variable_out), "\\.odd")
    ]
    even <- variable_out[
        stringr::str_subset(colnames(variable_out), "\\.even")
    ]

    #  Strip suffixes from column names
    colnames(odd) <- str_replace_all(colnames(odd), "\\.odd", "")
    colnames(even) <- str_replace_all(colnames(even), "\\.even", "")

    #  Odd tibble is tibble #1, even tibble is tibble #2
    odd$tibble <- "1"
    even$tibble <- "2"

    #  Again, interleave the odd/even rows, then arrange them by group ID and
    #+ tibble number; save results to initial "tbl.*" variables
    variable_out <- odd %>%
            dplyr::mutate(groupid = row_number()) %>%
            dplyr::bind_rows(even %>% mutate(groupid = row_number())) %>%
            dplyr::arrange(groupid, tibble) %>%
            # dplyr::relocate(
            #     groupid, .after = tidyselect::all_of(position_after)
            # ) %>%
            dplyr::select(-tibble)

    return(variable_out)
}


#  Set up work directory (locations TB∆) --------------------------------------
directory_user <- "/Users/kalavattam"
directory_base <- "Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do"
directory_work <- "get_unique_fragments/Bonora"

setwd(
    paste0(directory_user, "/", directory_base, "/", directory_work)
)
rm(directory_user, directory_base, directory_work)


#  Having run part 3 first, load part-3 .Rdata into environment ---------------
load("test.PE-processing.part-4.Rdata")
#  Doing so loads in appropriate functions, variables, etc.


#  Script name ----------------------------------------------------------------
script <- "test.PE-processing.part-5.R"


#  Create "minimal tibbles" ("m") of only ... and AS values -------------------
#+ 
#+ Rename the AS columns to identify their alignment of origin
variable <- paste0("uniline.tbl.", c("129S1", "CAST", "mm10"))
variable_m <- paste0("m.", stringr::str_remove(variable, "uniline.tbl."))
suffix <- variable_m %>%
    strsplit(., "\\.") %>%
    lapply(., `[[`, 2) %>%
    unlist()

command <- paste0(
    "<- ", variable, " %>% ",
        "dplyr::rename(AS.odd = tag.AS.odd)", " %>% ",
        "dplyr::rename(AS.even = tag.AS.even)", " %>% ",
        "dplyr::rename(MD.odd = tag.MD.odd)", " %>% ",
        "dplyr::rename(MD.even = tag.MD.even)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

#  129, CAST --------------------------
command <- paste0(
    "<- ",
        variable[stringr::str_detect(variable, "mm10", negate = TRUE)], " %>% ",
        "dplyr::select(",
            "coordinate.odd, groupid.odd, AS.odd, ",
            "coordinate.even, groupid.even, AS.even, ",
            "flag.odd, seq.odd, flag.even, seq.even, ",
            "rname.odd, pos.odd, pos_end.odd, ",
            "mrnm.odd, mpos.odd, mpos_end.odd, ",
            "rname.even, pos.even, pos_end.even, ",
            "mrnm.even, mpos.even, mpos_end.even, ",
            "lO_rname.odd, lO_pos.odd, lO_pos_end.odd, ",
            "lO_mrnm.odd, lO_mpos.odd, lO_mpos_end.odd, ",
            "lO_rname.even, lO_pos.even, lO_pos_end.even, ",
            "lO_mrnm.even, lO_mpos.even, lO_mpos_end.even",
        ")", " %>% ",
        "dplyr::rename(",
            "AS.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = AS.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "AS.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = AS.even",
        ")", " %>% ",
        "dplyr::rename(",
            "groupid.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = groupid.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "groupid.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = groupid.even",
        ")", " %>% ",
        "dplyr::rename(",
            "flag.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = flag.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "flag.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = flag.even",
        ")", " %>% ",
        "dplyr::rename(",
            "seq.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = seq.odd",
        ")", " %>% ",
        "dplyr::rename(", 
            "seq.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = seq.even",
        ")", " %>% ",
        "dplyr::rename(",
            "rname.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = rname.odd",
        ")", " %>% ",
        "dplyr::rename(", 
            "rname.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = rname.even",
        ")", " %>% ",
        "dplyr::rename(",
            "pos.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = pos.odd",
        ")", " %>% ",
        "dplyr::rename(", 
            "pos.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = pos.even",
        ")", " %>% ",
        "dplyr::rename(",
            "pos_end.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = pos_end.odd",
        ")", " %>% ",
        "dplyr::rename(", 
            "pos_end.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = pos_end.even",
        ")", " %>% ",
        "dplyr::rename(",
            "mrnm.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = mrnm.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "mrnm.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = mrnm.even",
        ")", " %>% ",
        "dplyr::rename(",
            "mpos.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = mpos.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "mpos.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = mpos.even",
        ")", " %>% ",
        "dplyr::rename(",
            "mpos_end.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = mpos_end.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "mpos_end.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = mpos_end.even",
        ")", " %>% ",
        "dplyr::rename(",
            "lO_rname.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = lO_rname.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "lO_rname.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = lO_rname.even",
        ")", " %>% ",
        "dplyr::rename(",
            "lO_pos.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = lO_pos.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "lO_pos.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = lO_pos.even",
        ")", " %>% ",
        "dplyr::rename(",
            "lO_pos_end.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = lO_pos_end.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "lO_pos_end.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = lO_pos_end.even",
        ")", " %>% ",
        "dplyr::rename(",
            "lO_mrnm.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = lO_mrnm.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "lO_mrnm.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = lO_mrnm.even",
        ")", " %>% ",
        "dplyr::rename(",
            "lO_mpos.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = lO_mpos.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "lO_mpos.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = lO_mpos.even",
        ")", " %>% ",
        "dplyr::rename(",
            "lO_mpos_end.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = lO_mpos_end.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "lO_mpos_end.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = TRUE)],
            " = lO_mpos_end.even",
        ")"
)
operation <- makeOperation(
    variable_m[stringr::str_detect(variable_m, "mm10", negate = TRUE)],
    command
)
evaluateOperation(operation)

#  mm10 -------------------------------
command <- paste0(
    "<- ", variable[stringr::str_detect(variable, "mm10", negate = FALSE)], " %>% ",
        "dplyr::select(",
            "coordinate.odd, groupid.odd, AS.odd, ",
            "coordinate.even, groupid.even, AS.even, ",
            "flag.odd, seq.odd, flag.even, seq.even, ",
            "rname.odd, pos.odd, pos_end.odd, ",
            "mrnm.odd, mpos.odd, mpos_end.odd, ",
            "rname.even, pos.even, pos_end.even, ",
            "mrnm.even, mpos.even, mpos_end.even",
        ")", " %>% ",
        "dplyr::rename(",
            "AS.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = AS.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "AS.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = AS.even",
        ")", " %>% ",
        "dplyr::rename(",
            "groupid.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = groupid.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "groupid.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = groupid.even",
        ")", " %>% ",
        "dplyr::rename(",
            "flag.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = flag.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "flag.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = flag.even",
        ")", " %>% ",
        "dplyr::rename(",
            "seq.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = seq.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "seq.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = seq.even",
        ")", " %>% ",
        "dplyr::rename(",
            "rname.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = rname.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "rname.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = rname.even",
        ")", " %>% ",
        "dplyr::rename(",
            "pos.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = pos.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "pos.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = pos.even",
        ")", " %>% ",
        "dplyr::rename(",
            "pos_end.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = pos_end.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "pos_end.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = pos_end.even",
        ")", " %>% ",
        "dplyr::rename(",
            "mrnm.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = mrnm.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "mrnm.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = mrnm.even",
        ")", " %>% ",
        "dplyr::rename(",
            "mpos.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = mpos.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "mpos.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = mpos.even",
        ")", " %>% ",
        "dplyr::rename(",
            "mpos_end.odd.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = mpos_end.odd",
        ")", " %>% ",
        "dplyr::rename(",
            "mpos_end.even.",
            suffix[stringr::str_detect(suffix, "mm10", negate = FALSE)],
            " = mpos_end.even",
        ")"
)
operation <- makeOperation(
    variable_m[stringr::str_detect(variable_m, "mm10", negate = FALSE)],
    command
)
evaluateOperation(operation)

#  Sort and deconcatenate uniline tibbles for 129, CAST
variable_b <- paste0("b.", suffix)
command <- paste0(" <- \"\"")
operation <- makeOperation(variable_b, command)
evaluateOperation(operation)

command <- paste0(
    "<- ", "mungeMatesIntoTwoRows(\"",
        variable_m, "\", ",
        "\"coordinate.odd\"",
    ")"
)
operation <- makeOperation(variable_b, command)
evaluateOperation(operation)


#  Fully join the "u" tibbles by ... ------------------------------------------
b.full <- full_join(b.mm10, b.129S1, by = "coordinate") %>%
    full_join(., b.CAST, by = "coordinate") %>%
    dplyr::rename("groupid.2.mm10" = "groupid.x") %>%
    dplyr::rename("groupid.2.129S1" = "groupid.y") %>%
    dplyr::rename("groupid.2.CAST" = "groupid") %>%
    dplyr::relocate(
        paste0(
            "groupid.",
            c("mm10", "2.mm10", "129S1", "2.129S1", "CAST", "2.CAST")
        ),
        .after = "lO_mpos_end.CAST"
    )


#  b.full tibble: Concatenate even row to preceding odd row -------------------
variable_uniline <- "uniline.b.full"
mapply(
    assign, variable_uniline, "", MoreArgs = list(envir = parent.frame())
)

uniline.b.full <- mungeMatesIntoOneRow("b.full", "coordinate")


#  For readability, relocate columns in the tibble
uniline.b.full <- uniline.b.full %>%
    dplyr::relocate(
        c(
            AS.mm10.odd, AS.mm10.even,
            AS.129S1.odd, AS.129S1.even,
            AS.CAST.odd, AS.CAST.even
        ),
        .before = coordinate.odd) %>%
    dplyr::relocate(coordinate.even, .after = coordinate.odd) %>%
    dplyr::relocate(
        c(
            flag.mm10.odd, seq.mm10.odd,
            flag.mm10.even, seq.mm10.even,
            flag.129S1.odd, seq.129S1.odd,
            flag.129S1.even, seq.129S1.even,
            flag.CAST.odd, seq.CAST.odd,
            flag.CAST.even, seq.CAST.even
        ),
        .after = coordinate.even
    ) %>%
    dplyr::relocate(c(groupid.even, groupid.odd), .after = seq.CAST.even) %>%
    dplyr::rename(ID.even = groupid.even ) %>%
    dplyr::rename(ID.odd = groupid.odd) %>%
    dplyr::select(-contains("groupid.")) %>%
    dplyr::rename(groupid.even = ID.even) %>%
    dplyr::rename(groupid.odd = ID.odd)


#  Take the AS min for each fragment ------------------------------------------
#+ 
#+ Taking the AS min (pmin()) gives us the worst AS score per read pair; our
#+ logic is that the read pair should be represented by the worst of its two
#+ AS's

#  mm10
command <- paste0(
    "<- ", "pmin(",
        variable_uniline, "$AS.mm10.odd, ", variable_uniline, "$AS.mm10.even",
    ")"
)
operation <- makeOperation(paste0(variable_uniline, "$AS.mm10.pmin"), command)
evaluateOperation(operation)

#  129S1
command <- paste0(
    "<- ", "pmin(",
        variable_uniline, "$AS.129S1.odd, ", variable_uniline, "$AS.129S1.even",
    ")"
)
operation <- makeOperation(paste0(variable_uniline, "$AS.129S1.pmin"), command)
evaluateOperation(operation)

#  CAST
command <- paste0(
    "<- ", "pmin(",
        variable_uniline, "$AS.CAST.odd, ", variable_uniline, "$AS.CAST.even",
    ")"
)
operation <- makeOperation(paste0(variable_uniline, "$AS.CAST.pmin"), command)
evaluateOperation(operation)

#  Check
(uniline.b.full$seq.mm10.odd == uniline.b.full$seq.129S1.odd) %>%
    table(useNA = "ifany")
(uniline.b.full$seq.mm10.odd == uniline.b.full$seq.CAST.odd) %>%
    table(useNA = "ifany")
(uniline.b.full$seq.129S1.odd == uniline.b.full$seq.CAST.odd) %>%
    table(useNA = "ifany")
(uniline.b.full$seq.mm10.even == uniline.b.full$seq.129S1.even) %>%
    table(useNA = "ifany")
(uniline.b.full$seq.mm10.even == uniline.b.full$seq.CAST.even) %>%
    table(useNA = "ifany")
(uniline.b.full$seq.129S1.even == uniline.b.full$seq.CAST.even) %>%
    table(useNA = "ifany")


#  Create minimal tibbles for AS.*.pmin values --------------------------------
AS.pmin <- uniline.b.full %>%
    dplyr::select(
        coordinate.odd, coordinate.even,
        AS.mm10.pmin, AS.129S1.pmin, AS.CAST.pmin,
        flag.mm10.odd, flag.mm10.even,
        seq.mm10.odd, seq.mm10.even,
        flag.129S1.odd, flag.129S1.even,
        seq.129S1.odd, seq.129S1.even,
        flag.CAST.odd, flag.CAST.even,
        seq.CAST.odd, seq.CAST.even,
        rname.mm10.odd, pos.mm10.odd, pos_end.mm10.odd,
        rname.mm10.even, pos.mm10.even, pos_end.mm10.even,
        rname.129S1.odd, pos.129S1.odd, pos_end.129S1.odd,
        rname.129S1.even, pos.129S1.even, pos_end.129S1.even,
        lO_rname.129S1.odd, lO_pos.129S1.odd, lO_pos_end.129S1.odd,
        lO_rname.129S1.even, lO_pos.129S1.even, lO_pos_end.129S1.even,
        rname.CAST.odd, pos.CAST.odd, pos_end.CAST.odd,
        rname.CAST.even, pos.CAST.even, pos_end.CAST.even,
        lO_rname.CAST.odd, lO_pos.CAST.odd, lO_pos_end.CAST.odd,
        lO_rname.CAST.even, lO_pos.CAST.even, lO_pos_end.CAST.even,
        mrnm.mm10.odd, mpos.mm10.odd, mpos_end.mm10.odd,
        mrnm.mm10.even, mpos.mm10.even, mpos_end.mm10.even,
        mrnm.129S1.odd, mpos.129S1.odd, mpos_end.129S1.odd,
        mrnm.129S1.even, mpos.129S1.even, mpos_end.129S1.even,
        lO_mrnm.129S1.odd, lO_mpos.129S1.odd, lO_mpos_end.129S1.odd,
        lO_mrnm.129S1.even, lO_mpos.129S1.even, lO_mpos_end.129S1.even,
        mrnm.CAST.odd, mpos.CAST.odd, mpos_end.CAST.odd,
        mrnm.CAST.even, mpos.CAST.even, mpos_end.CAST.even,
        lO_mrnm.CAST.odd, lO_mpos.CAST.odd, lO_mpos_end.CAST.odd,
        lO_mrnm.CAST.even, lO_mpos.CAST.even, lO_mpos_end.CAST.even
    )
colnames(AS.pmin) <- gsub(".pmin", "", colnames(AS.pmin))


#  Calculate the differences between 129 and CAST alignment scores ------------
#+ 
#+ Initially assign categories to the difference based on values w/r/t/the
#+ variable int (see below)
int <- 0L
AS.pmin <- AS.pmin %>% 
    dplyr::mutate(difference = AS.129S1 - AS.CAST) %>% 
    dplyr::mutate(
        assignment.initial = case_when(
            difference >= (-1 * int) & difference <= int ~ "Neutral",
            difference > int ~ "129S1-SvImJ",
            difference < (-1 * int) ~ "CAST-EiJ"
        )
    )
rm(int)


#  Formal assignments of categories: See the block of code --------------------
AS.pmin$assignment <- ifelse(
    is.na(AS.pmin$assignment.initial),
    ifelse(
        !is.na(AS.pmin$AS.129S1),
        "129S1-SvImJ",
        ifelse(
            !is.na(AS.pmin$AS.CAST),
            "CAST-EiJ",
            AS.pmin$assignment.initial
        )
    ),
    AS.pmin$assignment.initial
) %>%
    forcats::as_factor()


#  Remove all but read, AS, and formal assignment columns ---------------------
AS.pmin <- AS.pmin %>%
    dplyr::select(
        coordinate.odd, coordinate.even,
        AS.mm10, AS.129S1, AS.CAST, assignment,
        flag.mm10.odd, flag.mm10.even,
        seq.mm10.odd, seq.mm10.even,
        flag.129S1.odd, flag.129S1.even,
        seq.129S1.odd, seq.129S1.even,
        flag.CAST.odd, flag.CAST.even,
        seq.CAST.odd, seq.CAST.even,
        rname.mm10.odd, pos.mm10.odd, pos_end.mm10.odd,
        rname.mm10.even, pos.mm10.even, pos_end.mm10.even,
        rname.129S1.odd, pos.129S1.odd, pos_end.129S1.odd,
        rname.129S1.even, pos.129S1.even, pos_end.129S1.even,
        lO_rname.129S1.odd, lO_pos.129S1.odd, lO_pos_end.129S1.odd,
        lO_rname.129S1.even, lO_pos.129S1.even, lO_pos_end.129S1.even,
        rname.CAST.odd, pos.CAST.odd, pos_end.CAST.odd,
        rname.CAST.even, pos.CAST.even, pos_end.CAST.even,
        lO_rname.CAST.odd, lO_pos.CAST.odd, lO_pos_end.CAST.odd,
        lO_rname.CAST.even, lO_pos.CAST.even, lO_pos_end.CAST.even,
        mrnm.mm10.odd, mpos.mm10.odd, mpos_end.mm10.odd,
        mrnm.mm10.even, mpos.mm10.even, mpos_end.mm10.even,
        mrnm.129S1.odd, mpos.129S1.odd, mpos_end.129S1.odd,
        mrnm.129S1.even, mpos.129S1.even, mpos_end.129S1.even,
        lO_mrnm.129S1.odd, lO_mpos.129S1.odd, lO_mpos_end.129S1.odd,
        lO_mrnm.129S1.even, lO_mpos.129S1.even, lO_mpos_end.129S1.even,
        mrnm.CAST.odd, mpos.CAST.odd, mpos_end.CAST.odd,
        mrnm.CAST.even, mpos.CAST.even, mpos_end.CAST.even,
        lO_mrnm.CAST.odd, lO_mpos.CAST.odd, lO_mpos_end.CAST.odd,
        lO_mrnm.CAST.even, lO_mpos.CAST.even, lO_mpos_end.CAST.even
    )


#  Assign 0/1 categories for what aligned to what -----------------------------
AS.pmin$ref.mm10 <- ifelse(is.na(AS.pmin$AS.mm10), "0", "1")
AS.pmin$ref.129S1 <- ifelse(is.na(AS.pmin$AS.129S1), "0", "1")
AS.pmin$ref.CAST <- ifelse(is.na(AS.pmin$AS.CAST), "0", "1")

AS.pmin$trinary <- paste0(
    AS.pmin$ref.mm10,
    AS.pmin$ref.129S1,
    AS.pmin$ref.CAST
) %>%
    forcats::as_factor()

AS.pmin$assignment_trinary <- paste(
    AS.pmin$assignment,
    AS.pmin$trinary
) %>%
    forcats::as_factor()

AS.pmin$trinary.r <- AS.pmin$trinary
AS.pmin$assignment_trinary.r <- AS.pmin$assignment_trinary

AS.pmin$trinary.r <- AS.pmin$trinary.r %>% plyr::revalue(  # old_name = new_name
    .,
    c(
        "100" = "N-masked mm10 reference",
        "111" = "All three references",
        "110" = "N-masked mm10 and 129 references",
        "101" = "N-masked mm10 and CAST references",
        "011" = "129 and CAST references",
        "010" = "129 reference",
        "001" = "CAST reference"
    )
)
AS.pmin$assignment_trinary.r <- AS.pmin$assignment_trinary.r %>% plyr::revalue(
    .,
    c(
        "NA 100" = "assignment: NA; alignment: N-masked mm10 reference only",
        "Neutral 111" = "assignment: Neutral; alignment: all three references",
        "129S1-SvImJ 111" = "assignment: 129; alignment: all three references",
        "CAST-EiJ 111" = "assignment: CAST; alignment: all three references",
        "129S1-SvImJ 110" = "assignment: 129; alignment: N-masked mm10 and 129 references",
        "CAST-EiJ 101" = "assignment: CAST; alignment: N-masked mm10 and CAST references",
        "CAST-EiJ 011" = "assignment: CAST; alignment: 129 and CAST references",
        "129S1-SvImJ 010" = "assignment: 129; alignment: 129 reference only",
        "Neutral 011" = "assignment: Neutral; alignment: 129 and CAST references",
        "129S1-SvImJ 011" = "assignment: 129; alignment: 129 and CAST references",
        "CAST-EiJ 001" = "assignment: CAST; alignment: CAST reference only"
    )
)


#  Remove all but read, AS, formal assignment, and "trinary" columns ----------
AS.pmin <- AS.pmin %>%
    select(
        coordinate.odd, coordinate.even,
        AS.mm10, AS.129S1, AS.CAST, assignment,
        trinary, trinary.r,
        assignment_trinary, assignment_trinary.r,
        flag.mm10.odd, flag.mm10.even,
        seq.mm10.odd, seq.mm10.even,
        flag.129S1.odd, flag.129S1.even,
        seq.129S1.odd, seq.129S1.even,
        flag.CAST.odd, flag.CAST.even,
        seq.CAST.odd, seq.CAST.even,
        rname.mm10.odd, pos.mm10.odd, pos_end.mm10.odd,
        rname.mm10.even, pos.mm10.even, pos_end.mm10.even,
        rname.129S1.odd, pos.129S1.odd, pos_end.129S1.odd,
        rname.129S1.even, pos.129S1.even, pos_end.129S1.even,
        lO_rname.129S1.odd, lO_pos.129S1.odd, lO_pos_end.129S1.odd,
        lO_rname.129S1.even, lO_pos.129S1.even, lO_pos_end.129S1.even,
        rname.CAST.odd, pos.CAST.odd, pos_end.CAST.odd,
        rname.CAST.even, pos.CAST.even, pos_end.CAST.even,
        lO_rname.CAST.odd, lO_pos.CAST.odd, lO_pos_end.CAST.odd,
        lO_rname.CAST.even, lO_pos.CAST.even, lO_pos_end.CAST.even,
        mrnm.mm10.odd, mpos.mm10.odd, mpos_end.mm10.odd,
        mrnm.mm10.even, mpos.mm10.even, mpos_end.mm10.even,
        mrnm.129S1.odd, mpos.129S1.odd, mpos_end.129S1.odd,
        mrnm.129S1.even, mpos.129S1.even, mpos_end.129S1.even,
        lO_mrnm.129S1.odd, lO_mpos.129S1.odd, lO_mpos_end.129S1.odd,
        lO_mrnm.129S1.even, lO_mpos.129S1.even, lO_mpos_end.129S1.even,
        mrnm.CAST.odd, mpos.CAST.odd, mpos_end.CAST.odd,
        mrnm.CAST.even, mpos.CAST.even, mpos_end.CAST.even,
        lO_mrnm.CAST.odd, lO_mpos.CAST.odd, lO_mpos_end.CAST.odd,
        lO_mrnm.CAST.even, lO_mpos.CAST.even, lO_mpos_end.CAST.even
    )


#  Save the environment -------------------------------------------------------
save.image(stringr::str_remove(script, ".R") %>% paste0(., ".Rdata"))

rm(list = ls())

#  Before moving on to 'test.PE-processing.part-6.R', process .bam files made in
#+ 'test.PE-processing.part-4.R'

#  For example, for processing, use 'run_GB-assignment-pipeline.2022-0207.sh',
#+ 'generate-sam2pairwise-files_4dn-mouse-cross_GB.sh', and
#+ 'generate-sam2pairwise-files_4dn-mouse-cross_KA.sh'


# #  Create individual tibbles for each possible "assignment" value -------------
# AS.pmin$assignment <- AS.pmin$assignment %>% as.character()
# AS.pmin$assignment[is.na(AS.pmin$assignment)] <- "NA"
# variable <- list(
#     c("AS.pmin.assignment.129S1", "129S1-SvImJ"),
#     c("AS.pmin.assignment.CAST", "CAST-EiJ"),
#     c("AS.pmin.assignment.Neutral", "Neutral"),
#     c("AS.pmin.assignment.NA", "NA")
# )
# command <- paste0(
#     "<- AS.pmin %>% ",
#         "dplyr::filter_at(",
#             "vars(assignment), ",
#             "any_vars(. %in% \"", lapply(variable, `[[`, 2), "\")",
#         ")"
# )
# operation <- makeOperation(lapply(variable, `[[`, 1), command)
# evaluateOperation(operation)
# 
# AS.pmax$assignment <- AS.pmax$assignment %>% as.character()
# AS.pmax$assignment[is.na(AS.pmax$assignment)] <- "NA"
# variable <- list(
#     c("AS.pmax.assignment.129S1", "129S1-SvImJ"),
#     c("AS.pmax.assignment.CAST", "CAST-EiJ"),
#     c("AS.pmax.assignment.Neutral", "Neutral"),
#     c("AS.pmax.assignment.NA", "NA")
# )
# command <- paste0(
#     "<- AS.pmax %>% ",
#         "dplyr::filter_at(",
#             "vars(assignment), ",
#             "any_vars(. %in% \"", lapply(variable, `[[`, 2), "\")",
#         ")"
# )
# operation <- makeOperation(lapply(variable, `[[`, 1), command)
# evaluateOperation(operation)
# 
# 
# #  Make $assignment factors again ---------------------------------------------
# variable <- list(
#     c("AS.pmin.assignment.129S1", "129S1-SvImJ"),
#     c("AS.pmin.assignment.CAST", "CAST-EiJ"),
#     c("AS.pmin.assignment.Neutral", "Neutral"),
#     c("AS.pmin.assignment.NA", "NA")
# )
# command <- paste0(
#     "<- ", lapply(variable, `[[`, 1), "$assignment %>% ",
#         "as.factor()"
# )
# operation <- makeOperation(
#     paste0(lapply(variable, `[[`, 1), "$assignment"), command
# )
# evaluateOperation(operation)
# 
# AS.pmin.assignment.NA$assignment <- NA
# AS.pmin.assignment.NA$assignment <- NA_character_
# 
# variable <- list(
#     c("AS.pmax.assignment.129S1", "129S1-SvImJ"),
#     c("AS.pmax.assignment.CAST", "CAST-EiJ"),
#     c("AS.pmax.assignment.Neutral", "Neutral"),
#     c("AS.pmax.assignment.NA", "NA")
# )
# command <- paste0(
#     "<- ", lapply(variable, `[[`, 1), "$assignment %>% ",
#     "as.factor()"
# )
# operation <- makeOperation(
#     paste0(lapply(variable, `[[`, 1), "$assignment"), command
# )
# evaluateOperation(operation)
# 
# AS.pmax.assignment.NA$assignment <- NA
# AS.pmax.assignment.NA$assignment <- NA_character_
# 
# rm(variable)


# # -----------------------------------------------------------------------------
# ggplot(AS.pmin, aes(x = trinary.r)) +
#     geom_bar(alpha = 0.5) +
#     geom_text(stat = "count", aes(label = scales::comma(..count..)), size = 3, vjust = -0.5) +
#     ggtitle(returnObjectName(AS.pmin)) +
#     ylab("") +
#     xlab("fragment alignment status") +
#     scale_y_continuous(labels = comma) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
# ggplot(AS.pmax, aes(x = trinary.r)) +
#     geom_bar(alpha = 0.5) +
#     geom_text(stat = "count", aes(label = scales::comma(..count..)), size = 3, vjust = -0.5) +
#     ggtitle(returnObjectName(AS.pmax)) +
#     ylab("") +
#     xlab("fragment alignment status") +
#     scale_y_continuous(labels = comma) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
#
# ggplot(AS.pmin, aes(x = assignment)) +
#     geom_bar(alpha = 0.5) +
#     geom_text(stat = "count", aes(label = scales::comma(..count..)), size = 3, vjust = -0.5) +
#     ggtitle(returnObjectName(AS.pmin)) +
#     ylab("") +
#     xlab("fragment reference assignment") +
#     scale_y_continuous(labels = comma) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
# ggplot(AS.pmax, aes(x = assignment)) +
#     geom_bar(alpha = 0.5) +
#     geom_text(stat = "count", aes(label = scales::comma(..count..)), size = 3, vjust = -0.5) +
#     ggtitle(returnObjectName(AS.pmax)) +
#     ylab("") +
#     xlab("fragment reference assignment") +
#     scale_y_continuous(labels = comma) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
#
# ggplot(AS.pmin, aes(x = assignment_trinary.r)) +
#     geom_bar(alpha = 0.5) +
#     geom_text(stat = "count", aes(label = scales::comma(..count.., accuracy = NULL)), size = 2.5, vjust = -0.5) +  #TODO Figure out why decimals are showing up
#     ggtitle(returnObjectName(AS.pmin)) +
#     ylab("") +
#     xlab("fragment reference assignment and fragment alignment status") +
#     scale_y_continuous(labels = comma) +
#     theme(axis.text.x = element_text(angle = -90, hjust = 0))
#
# ggplot(AS.pmax, aes(x = assignment_trinary.r)) +
#     geom_bar(alpha = 0.5) +
#     geom_text(stat = "count", aes(label = scales::comma(..count.., accuracy = NULL)), size = 2.5, vjust = -0.5) +  #TODO Figure out why decimals are showing up
#     ggtitle(returnObjectName(AS.pmax)) +
#     ylab("") +
#     xlab("fragment reference assignment and fragment alignment status") +
#     scale_y_continuous(labels = comma) +
#     theme(axis.text.x = element_text(angle = -90, hjust = 0))
#
#
# # -----------------------------------------------------------------------------
# ggplot(AS.pmin, aes(x = assignment, color = trinary.r, fill = trinary.r)) +
#     geom_bar(alpha = 0.5, width = 0.9) +
#     ggtitle(returnObjectName(AS.pmin)) +
#     ylab("") +
#     xlab("assignment") +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma)
##
# ggplot(AS.pmax, aes(x = assignment, color = trinary.r, fill = trinary.r)) +
#     geom_bar(alpha = 0.5, width = 0.9) +
#     ggtitle(returnObjectName(AS.pmax)) +
#     ylab("") +
#     xlab("assignment") +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma)
#
#
# #  Cool
# ggplot(AS.pmin, aes(x = assignment, color = trinary.r, fill = trinary.r)) +
#     geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
#     ggtitle(returnObjectName(AS.pmin)) +
#     ylab("") +
#     xlab("assignment") +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma)
#
# ggplot(AS.pmax, aes(x = assignment, color = trinary.r, fill = trinary.r)) +
#     geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
#     ggtitle(returnObjectName(AS.pmax)) +
#     ylab("") +
#     xlab("assignment") +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma)
#
#
# ggplot(AS.pmin, aes(x = trinary.r, color = assignment, fill = assignment)) +
#     geom_bar(alpha = 0.5, width = 0.9) +
#     ggtitle(returnObjectName(AS.pmin)) +
#     ylab("") +
#     xlab("alignment") +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
# ggplot(AS.pmax, aes(x = trinary.r, color = assignment, fill = assignment)) +
#     geom_bar(alpha = 0.5, width = 0.9) +
#     ggtitle(returnObjectName(AS.pmax)) +
#     ylab("") +
#     xlab("alignment") +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
#
# #  Cool
# ggplot(AS.pmin, aes(x = trinary.r, color = assignment, fill = assignment)) +
#     geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
#     ggtitle(returnObjectName(AS.pmin)) +
#     ylab("") +
#     xlab("assembly-alignment combination") +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
# ggplot(AS.pmax, aes(x = trinary.r, color = assignment, fill = assignment)) +
#     geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
#     ggtitle(returnObjectName(AS.pmax)) +
#     ylab("") +
#     xlab("assembly-alignment combination") +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))


# # -----------------------------------------------------------------------------
# # -------
# ggplot(AS.pmin, aes(x = assignment, fill = as_factor(AS.CAST))) +
#     geom_bar(alpha = 0.9, width = 0.9) +
#     ggtitle("AS from alignment to CAST") +
#     ggtitle(returnObjectName(AS.pmin)) +
#     ylab("") +
#     xlab("assignment") +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma)
#
# ggplot(AS.pmax, aes(x = assignment, fill = as_factor(AS.CAST))) +
#     geom_bar(alpha = 0.9, width = 0.9) +
#     ggtitle("AS from alignment to CAST") +
#     ggtitle(returnObjectName(AS.pmax)) +
#     ylab("") +
#     xlab("assignment") +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma)
#
#
# # -------
# ggplot(AS.pmin, aes(x = assignment, fill = as_factor(AS.129S1))) +
#     geom_bar(alpha = 0.9, width = 0.9) +
#     ggtitle(paste0(returnObjectName(AS.pmin), ": AS from alignment to 129")) +
#     ylab("") +
#     xlab("assignment") +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma)
#
# ggplot(AS.pmax, aes(x = assignment, fill = as_factor(AS.129S1))) +
#     geom_bar(alpha = 0.9, width = 0.9) +
#     ggtitle(paste0(returnObjectName(AS.pmax), ": AS from alignment to 129")) +
#     ylab("") +
#     xlab("assignment") +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma)
#
#
# # -------
# ggplot(AS.pmin, aes(x = assignment, fill = as_factor(AS.mm10))) +
#     geom_bar(alpha = 0.9, width = 0.9) +
#     ggtitle(paste0(returnObjectName(AS.pmin), ": AS from alignment to mm10-CAST-129-N-masked")) +
#     ylab("") +
#     xlab("assignment") +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma)
#
# ggplot(AS.pmax, aes(x = assignment, fill = as_factor(AS.mm10))) +
#     geom_bar(alpha = 0.9, width = 0.9) +
#     ggtitle(paste0(returnObjectName(AS.pmax), ": AS from alignment to mm10-CAST-129-N-masked")) +
#     ylab("") +
#     xlab("assignment") +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma)
#
#
# # -----------------------------------------------------------------------------
# # ggplot(AS.pmin, aes(x = CAST.AS, color = assignment)) +
# #     geom_histogram(binwidth = 1, fill = "white", alpha = 0.5)
# # ggplot(AS.pmin, aes(x = oneTwentyNine.AS, color = assignment)) +
# #     geom_histogram(binwidth = 1, fill = "white", alpha = 0.5)
# # ggplot(AS.pmin, aes(x = mm10_CAST_129_Nmasked.AS, color = assignment)) +
# #     geom_histogram(binwidth = 1, fill = "white", alpha = 0.5)
#
# ggplot(AS.pmin, aes(x = AS.CAST, color = assignment, fill = assignment)) +
#     geom_histogram(binwidth = 1, alpha = 0.5) +
#     ggtitle("AS from alignment to CAST")
# ggplot(AS.pmin, aes(x = AS.129S1, color = assignment, fill = assignment)) +
#     geom_histogram(binwidth = 1, alpha = 0.5) +
#     ggtitle("AS from alignment to 129")
# ggplot(AS.pmin, aes(x = AS.mm10, color = assignment, fill = assignment)) +
#     geom_histogram(binwidth = 1, alpha = 0.5) +
#     ggtitle("AS from alignment to mm10-CAST-129-N-masked")
#
# ggplot(AS.pmin, aes(x = AS.CAST, color = assignment, fill = assignment)) +
#     geom_histogram(binwidth = 1, alpha = 0.5) +
#     ggtitle("AS from alignment to CAST") +
#     facet_grid(assignment ~ .)
# ggplot(AS.pmin, aes(x = AS.129S1, color = assignment, fill = assignment)) +
#     geom_histogram(binwidth = 1, alpha = 0.5) +
#     ggtitle("AS from alignment to 129") +
#     facet_grid(assignment ~ .)
# ggplot(AS.pmin, aes(x = AS.mm10, color = assignment, fill = assignment)) +
#     geom_histogram(binwidth = 1, alpha = 0.5) +
#     ggtitle("AS from alignment to mm10-CAST-129-N-masked") +
#     facet_grid(assignment ~ .)
#
#
# ggplot(AS.pmax, aes(x = AS.CAST, color = assignment, fill = assignment)) +
#     geom_histogram(binwidth = 1, alpha = 0.5) +
#     ggtitle("AS from alignment to CAST")
# ggplot(AS.pmax, aes(x = AS.129S1, color = assignment, fill = assignment)) +
#     geom_histogram(binwidth = 1, alpha = 0.5) +
#     ggtitle("AS from alignment to 129")
# ggplot(AS.pmax, aes(x = AS.mm10, color = assignment, fill = assignment)) +
#     geom_histogram(binwidth = 1, alpha = 0.5) +
#     ggtitle("AS from alignment to mm10-CAST-129-N-masked")
#
# ggplot(AS.pmax, aes(x = AS.CAST, color = assignment, fill = assignment)) +
#     geom_histogram(binwidth = 1, alpha = 0.5) +
#     ggtitle("AS from alignment to CAST") +
#     facet_grid(assignment ~ .)
# ggplot(AS.pmax, aes(x = AS.129S1, color = assignment, fill = assignment)) +
#     geom_histogram(binwidth = 1, alpha = 0.5) +
#     ggtitle("AS from alignment to 129") +
#     facet_grid(assignment ~ .)
# ggplot(AS.pmax, aes(x = AS.mm10, color = assignment, fill = assignment)) +
#     geom_histogram(binwidth = 1, alpha = 0.5) +
#     ggtitle("AS from alignment to mm10-CAST-129-N-masked") +
#     facet_grid(assignment ~ .)
