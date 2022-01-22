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
library(stringr)
library(tidyverse)
library(scales)

options(pillar.sigfig = 8, scipen = 10000)
set.seed(24)


#  Set up work directory (locations TB∆) --------------------------------------
directory_user <- "/Users/kalavattam"
directory_base <- "Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do"
directory_work <- "get_unique_fragments/Bonora"

setwd(
    paste0(directory_user, "/", directory_base, "/", directory_work)
)
rm(directory_user, directory_base, directory_work)


#  Having run part 3 first, load part-3 .Rdata into environment ---------------
load("test.PE-processing.2021-1221.2022-0108.part-4.Rdata")
#  Doing so loads in appropriate functions, variables, etc.


#  Script name ----------------------------------------------------------------
script <- "test.PE-processing.2021-1221.2022-0121.part-5.R"


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
        "dplyr::select(",
            "coordinate.odd, groupid.odd, AS.odd, ",
            "coordinate.even, groupid.even, AS.even",
        ")", " %>% ",
        "dplyr::rename(AS.odd.", suffix, " = AS.odd)", " %>% ",
        "dplyr::rename(AS.even.", suffix, " = AS.even)", " %>% ",
        "dplyr::rename(groupid.odd.", suffix, " = groupid.odd)", " %>% ",
        "dplyr::rename(groupid.even.", suffix, " = groupid.even)"
)
operation <- makeOperation(variable_m, command)
evaluateOperation(operation)

#  Sort and deconcatenate uniline tibbles for 129, CAST
variable_b <- paste0("b.", suffix)
command <- paste0(" <- \"\"")
operation <- makeOperation(variable_b, command)
evaluateOperation(operation)

for (i in 1:length(variable_m)) {
    df <- eval(parse(text = variable_m[i]))
    
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
    operation <- makeOperation(variable_b[i], command)
    evaluateOperation(operation)
    
    #  Again, check to make sure that there are no more than two entries per
    #+ group ID
    n_occur <- data.frame(table(eval(parse(text = variable_b[i]))$groupid))
    print(n_occur[n_occur$Freq > 2, ])
    
    #  Clean up
    rm(df, n_occur, even, odd)
}
rm(i)


#  Fully join the "u" tibbles by ... ------------------------------------------
b.full <- full_join(b.mm10, b.129S1, by = "coordinate") %>%
    full_join(., b.CAST, by = "coordinate") %>%
    dplyr::relocate(
        paste0("groupid.", c("mm10", "129S1", "CAST")), .after = coordinate
    )


#  b.full tibble: Concatenate even row to preceding odd row -------------------
variable_uniline <- "uniline.b.full"
mapply(
    assign, variable_uniline, "", MoreArgs = list(envir = parent.frame())
)

df <- b.full

odd.seq <- seq(1, nrow(df), 2)
even.seq <- (seq(1, nrow(df), 2)) + 1

#  Split df
odd <- df[odd.seq, ]
even <- df[even.seq, ]

#  Number each tibble
odd$tibble <- "1"
even$tibble <- "2"

#  Put back together, except now with $groupid
df <- odd %>%
    dplyr::mutate(groupid = row_number()) %>%
    dplyr::bind_rows(even %>% dplyr::mutate(groupid = row_number())) %>%
    dplyr::arrange(groupid, tibble) %>%
    dplyr::select(-tibble)

#  Split df again
odd <- df[odd.seq, ]
even <- df[even.seq, ]

df <- full_join(odd, even, by = "groupid") %>%
    dplyr::select(-groupid.x.x, -groupid.y.x, -groupid.x.y, -groupid.y.y)

colnames(df) <- str_replace_all(colnames(df), "\\.x", "\\.odd")
colnames(df) <- str_replace_all(colnames(df), "\\.y", "\\.even")

df$groupid.even <- df$groupid
df <- df %>%
    dplyr::rename(groupid.odd = groupid) %>%
    dplyr::relocate(groupid.odd, .after = coordinate.odd) %>%
    dplyr::relocate(groupid.even, .after = coordinate.even)

command <- paste0("<- ", "df")
operation <- makeOperation(variable_uniline, command)
evaluateOperation(operation)

#  Remove unneeded variables
rm(df, even, even.seq, odd, odd.seq)

#  For readability, relocate columns in the tibble
uniline.b.full <- uniline.b.full %>%
    dplyr::relocate(c(
        AS.mm10.odd, AS.mm10.even,
        AS.129S1.odd, AS.129S1.even,
        AS.CAST.odd, AS.CAST.even
    )) %>%
    dplyr::relocate(
        c(coordinate.odd, coordinate.even), .before = AS.mm10.odd
    ) %>%
    dplyr::relocate(c(groupid.odd, groupid.even), .after = AS.CAST.even) %>%
    dplyr::relocate(
        c(groupid.mm10.odd, groupid.mm10.even), .after = groupid.even
    ) %>%
    dplyr::relocate(
        c(groupid.129S1.odd, groupid.129S1.even), .after = groupid.mm10.even
    ) %>%
    dplyr::relocate(
        c(groupid.CAST.odd, groupid.CAST.even), .after = groupid.129S1.even
    )


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


# #  Take the AS max for each fragment ------------------------------------------
# #+ 
# #+ Taking the AS max (pmax()) gives us the best AS score per read pair
# 
# #  mm10
# command <- paste0(
#     "<- ", "pmax(",
#         variable_uniline, "$AS.mm10.odd, ",
#         variable_uniline, "$AS.mm10.even",
#     ")"
# )
# operation <- makeOperation(paste0(variable_uniline, "$AS.mm10.pmax"), command)
# evaluateOperation(operation)
# 
# #  129S1
# command <- paste0(
#     "<- ", "pmax(",
#         variable_uniline, "$AS.129S1.odd, ",
#         variable_uniline, "$AS.129S1.even",
#     ")"
# )
# operation <- makeOperation(paste0(variable_uniline, "$AS.129S1.pmax"), command)
# evaluateOperation(operation)
# 
# #  CAST
# command <- paste0(
#     "<- ", "pmax(",
#         variable_uniline, "$AS.CAST.odd, ",
#         variable_uniline, "$AS.CAST.even",
#     ")"
# )
# operation <- makeOperation(paste0(variable_uniline, "$AS.CAST.pmax"), command)
# evaluateOperation(operation)


#  Create minimal tibbles for AS.*.pmin and AS.*.pmax values
AS.pmin <- uniline.b.full %>%
    dplyr::select(
        coordinate.odd, coordinate.even, AS.mm10.pmin, AS.129S1.pmin, AS.CAST.pmin
    )
colnames(AS.pmin) <- gsub(".pmin", "", colnames(AS.pmin))

# AS.pmax <- uniline.m.full %>%
#     dplyr::select(
#         coordinate.odd, coordinate.even, AS.mm10.pmax, AS.129S1.pmax, AS.CAST.pmax
#     )
# colnames(AS.pmax) <- gsub(".pmax", "", colnames(AS.pmax))


#  Calculate the differences between 129 and CAST alignment scores ------------
#+ 
#+ Initially assign categories to the difference based on values w/r/t/the
#+ variable int (see below)
int <- 0 %>% as.integer()  #TODO Make the integer an argument
AS.pmin <- AS.pmin %>% 
    dplyr::mutate(difference = AS.129S1 - AS.CAST) %>% 
    dplyr::mutate(
        assignment.initial = case_when(
            difference >= (-1 * int) & difference <= int ~ "Neutral",
            difference > int ~ "129S1-SvImJ",
            difference < (-1 * int) ~ "CAST-EiJ"
        )
    )

# AS.pmax <- AS.pmax %>% 
#     dplyr::mutate(difference = AS.129S1 - AS.CAST) %>% 
#     dplyr::mutate(
#         assignment.initial = case_when(
#             difference >= (-1 * int) & difference <= int ~ "Neutral",
#             difference > int ~ "129S1-SvImJ",
#             difference < (-1 * int) ~ "CAST-EiJ"
#         )
#     )

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

# AS.pmax$assignment <- ifelse(
#     is.na(AS.pmax$assignment.initial),
#     ifelse(
#         !is.na(AS.pmax$AS.129S1),
#         "129S1-SvImJ",
#         ifelse(
#             !is.na(AS.pmax$AS.CAST),
#             "CAST-EiJ",
#             AS.pmax$assignment.initial
#         )
#     ),
#     AS.pmax$assignment.initial
# ) %>%
#     forcats::as_factor()


#  Remove all but read, AS, and formal assignment columns ---------------------
AS.pmin <- AS.pmin %>%
    dplyr::select(
        coordinate.odd, coordinate.even, AS.mm10, AS.129S1, AS.CAST, assignment
    )

# AS.pmax <- AS.pmax %>%
#     dplyr::select(
#         coordinate.odd, coordinate.even, AS.mm10, AS.129S1, AS.CAST, assignment
#     )


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

# AS.pmax$ref.mm10 <- ifelse(is.na(AS.pmax$AS.mm10), "0", "1")
# AS.pmax$ref.129S1 <- ifelse(is.na(AS.pmax$AS.129S1), "0", "1")
# AS.pmax$ref.CAST <- ifelse(is.na(AS.pmax$AS.CAST), "0", "1")
# 
# AS.pmax$trinary <- paste0(
#     AS.pmax$ref.mm10,
#     AS.pmax$ref.129S1,
#     AS.pmax$ref.CAST
# ) %>%
#     forcats::as_factor()
# 
# AS.pmax$assignment_trinary <- paste(
#     AS.pmax$assignment,
#     AS.pmax$trinary
# ) %>%
#     forcats::as_factor()
# 
# AS.pmax$trinary.r <- AS.pmax$trinary
# AS.pmax$assignment_trinary.r <- AS.pmax$assignment_trinary
# 
# AS.pmax$trinary.r <- AS.pmax$trinary.r %>% plyr::revalue(  # old_name = new_name
#     .,
#     c(
#         "100" = "N-masked mm10 reference",
#         "111" = "All three references",
#         "110" = "N-masked mm10 and 129 references",
#         "101" = "N-masked mm10 and CAST references",
#         "011" = "129 and CAST references",
#         "010" = "129 reference",
#         "001" = "CAST reference"
#     )
# )
# AS.pmax$assignment_trinary.r <- AS.pmax$assignment_trinary.r %>% plyr::revalue(
#     .,
#     c(
#         "NA 100" = "assignment: NA; alignment: N-masked mm10 reference only",
#         "Neutral 111" = "assignment: Neutral; alignment: all three references",
#         "129S1-SvImJ 111" = "assignment: 129; alignment: all three references",
#         "CAST-EiJ 111" = "assignment: CAST; alignment: all three references",
#         "129S1-SvImJ 110" = "assignment: 129; alignment: N-masked mm10 and 129 references",
#         "CAST-EiJ 101" = "assignment: CAST; alignment: N-masked mm10 and CAST references",
#         "CAST-EiJ 011" = "assignment: CAST; alignment: 129 and CAST references",
#         "129S1-SvImJ 010" = "assignment: 129; alignment: 129 reference only",
#         "Neutral 011" = "assignment: Neutral; alignment: 129 and CAST references",
#         "129S1-SvImJ 011" = "assignment: 129; alignment: 129 and CAST references",
#         "CAST-EiJ 001" = "assignment: CAST; alignment: CAST reference only"
#     )
# )


#  Remove all but read, AS, formal assignment, and "trinary" columns ----------
AS.pmin <- AS.pmin %>%
    select(
        coordinate.odd,
        coordinate.even,
        AS.mm10,
        AS.129S1,
        AS.CAST,
        assignment,
        trinary,
        trinary.r,
        assignment_trinary,
        assignment_trinary.r
    )


#  Save the environment -------------------------------------------------------
save.image(stringr::str_remove(script, ".R") %>% paste0(., ".Rdata"))


#  More stuff, TBD ------------------------------------------------------------
# AS.pmax <- AS.pmax %>%
#     select(
#         coordinate.odd,
#         coordinate.even,
#         AS.mm10,
#         AS.129S1,
#         AS.CAST,
#         assignment,
#         trinary,
#         trinary.r,
#         assignment_trinary,
#         assignment_trinary.r
#     )


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


# -----------------------------------------------------------------------------
ggplot(AS.pmin, aes(x = trinary.r)) +
    geom_bar(alpha = 0.5) +
    geom_text(stat = "count", aes(label = scales::comma(..count..)), size = 3, vjust = -0.5) +
    ggtitle(returnObjectName(AS.pmin)) +
    ylab("") +
    xlab("fragment alignment status") +
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(AS.pmax, aes(x = trinary.r)) +
    geom_bar(alpha = 0.5) +
    geom_text(stat = "count", aes(label = scales::comma(..count..)), size = 3, vjust = -0.5) +
    ggtitle(returnObjectName(AS.pmax)) +
    ylab("") +
    xlab("fragment alignment status") +
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(AS.pmin, aes(x = assignment)) +
    geom_bar(alpha = 0.5) +
    geom_text(stat = "count", aes(label = scales::comma(..count..)), size = 3, vjust = -0.5) +
    ggtitle(returnObjectName(AS.pmin)) +
    ylab("") +
    xlab("fragment reference assignment") +
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(AS.pmax, aes(x = assignment)) +
    geom_bar(alpha = 0.5) +
    geom_text(stat = "count", aes(label = scales::comma(..count..)), size = 3, vjust = -0.5) +
    ggtitle(returnObjectName(AS.pmax)) +
    ylab("") +
    xlab("fragment reference assignment") +
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(AS.pmin, aes(x = assignment_trinary.r)) +
    geom_bar(alpha = 0.5) +
    geom_text(stat = "count", aes(label = scales::comma(..count.., accuracy = NULL)), size = 2.5, vjust = -0.5) +  #TODO Figure out why decimals are showing up
    ggtitle(returnObjectName(AS.pmin)) +
    ylab("") +
    xlab("fragment reference assignment and fragment alignment status") +
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0))

ggplot(AS.pmax, aes(x = assignment_trinary.r)) +
    geom_bar(alpha = 0.5) +
    geom_text(stat = "count", aes(label = scales::comma(..count.., accuracy = NULL)), size = 2.5, vjust = -0.5) +  #TODO Figure out why decimals are showing up
    ggtitle(returnObjectName(AS.pmax)) +
    ylab("") +
    xlab("fragment reference assignment and fragment alignment status") +
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0))


# -----------------------------------------------------------------------------
ggplot(AS.pmin, aes(x = assignment, color = trinary.r, fill = trinary.r)) +
    geom_bar(alpha = 0.5, width = 0.9) +
    ggtitle(returnObjectName(AS.pmin)) +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

ggplot(AS.pmax, aes(x = assignment, color = trinary.r, fill = trinary.r)) +
    geom_bar(alpha = 0.5, width = 0.9) +
    ggtitle(returnObjectName(AS.pmax)) +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)


#  Cool
ggplot(AS.pmin, aes(x = assignment, color = trinary.r, fill = trinary.r)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ggtitle(returnObjectName(AS.pmin)) +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

ggplot(AS.pmax, aes(x = assignment, color = trinary.r, fill = trinary.r)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ggtitle(returnObjectName(AS.pmax)) +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)


ggplot(AS.pmin, aes(x = trinary.r, color = assignment, fill = assignment)) +
    geom_bar(alpha = 0.5, width = 0.9) +
    ggtitle(returnObjectName(AS.pmin)) +
    ylab("") +
    xlab("alignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(AS.pmax, aes(x = trinary.r, color = assignment, fill = assignment)) +
    geom_bar(alpha = 0.5, width = 0.9) +
    ggtitle(returnObjectName(AS.pmax)) +
    ylab("") +
    xlab("alignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


#  Cool
ggplot(AS.pmin, aes(x = trinary.r, color = assignment, fill = assignment)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ggtitle(returnObjectName(AS.pmin)) +
    ylab("") +
    xlab("assembly-alignment combination") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(AS.pmax, aes(x = trinary.r, color = assignment, fill = assignment)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ggtitle(returnObjectName(AS.pmax)) +
    ylab("") +
    xlab("assembly-alignment combination") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


# -----------------------------------------------------------------------------
# -------
ggplot(AS.pmin, aes(x = assignment, fill = as_factor(AS.CAST))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle("AS from alignment to CAST") +
    ggtitle(returnObjectName(AS.pmin)) +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

ggplot(AS.pmax, aes(x = assignment, fill = as_factor(AS.CAST))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle("AS from alignment to CAST") +
    ggtitle(returnObjectName(AS.pmax)) +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)


# -------
ggplot(AS.pmin, aes(x = assignment, fill = as_factor(AS.129S1))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle(paste0(returnObjectName(AS.pmin), ": AS from alignment to 129")) +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

ggplot(AS.pmax, aes(x = assignment, fill = as_factor(AS.129S1))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle(paste0(returnObjectName(AS.pmax), ": AS from alignment to 129")) +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)


# -------
ggplot(AS.pmin, aes(x = assignment, fill = as_factor(AS.mm10))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle(paste0(returnObjectName(AS.pmin), ": AS from alignment to mm10-CAST-129-N-masked")) +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

ggplot(AS.pmax, aes(x = assignment, fill = as_factor(AS.mm10))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle(paste0(returnObjectName(AS.pmax), ": AS from alignment to mm10-CAST-129-N-masked")) +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)


# -----------------------------------------------------------------------------
# ggplot(AS.pmin, aes(x = CAST.AS, color = assignment)) +
#     geom_histogram(binwidth = 1, fill = "white", alpha = 0.5)
# ggplot(AS.pmin, aes(x = oneTwentyNine.AS, color = assignment)) +
#     geom_histogram(binwidth = 1, fill = "white", alpha = 0.5)
# ggplot(AS.pmin, aes(x = mm10_CAST_129_Nmasked.AS, color = assignment)) +
#     geom_histogram(binwidth = 1, fill = "white", alpha = 0.5)

ggplot(AS.pmin, aes(x = AS.CAST, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to CAST")
ggplot(AS.pmin, aes(x = AS.129S1, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to 129")
ggplot(AS.pmin, aes(x = AS.mm10, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to mm10-CAST-129-N-masked")

ggplot(AS.pmin, aes(x = AS.CAST, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to CAST") +
    facet_grid(assignment ~ .)
ggplot(AS.pmin, aes(x = AS.129S1, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to 129") +
    facet_grid(assignment ~ .)
ggplot(AS.pmin, aes(x = AS.mm10, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to mm10-CAST-129-N-masked") +
    facet_grid(assignment ~ .)


ggplot(AS.pmax, aes(x = AS.CAST, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to CAST")
ggplot(AS.pmax, aes(x = AS.129S1, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to 129")
ggplot(AS.pmax, aes(x = AS.mm10, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to mm10-CAST-129-N-masked")

ggplot(AS.pmax, aes(x = AS.CAST, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to CAST") +
    facet_grid(assignment ~ .)
ggplot(AS.pmax, aes(x = AS.129S1, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to 129") +
    facet_grid(assignment ~ .)
ggplot(AS.pmax, aes(x = AS.mm10, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to mm10-CAST-129-N-masked") +
    facet_grid(assignment ~ .)
