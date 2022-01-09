
library(ggplot2)
library(scales)
library(tidyverse)

options(pillar.sigfig = 8, scipen = 10000)


# -----------------------------------------------------------------------------
#  Set up functions
`%notin%` <- Negate(`%in%`)


convertPercent <- function(x) {
    if(is.numeric(x)) {
        ifelse(is.na(x), x, paste0(round(x * 100L, 2), "%")) 
    } else x 
}


loadRDS <- readRDS


loadTibbleFromRDS <- function(variable, file) {
    command <- paste0(variable, " <- loadRDS(file = \"", file, "\")") %>% as.list()
    return(eval(parse(text = command), envir = globalenv()))
}


makeOperation <- function(variable, command) {
    operation <- paste(variable, command) %>% paste(., collapse = "; ")
    return(operation)
}


mungeTmpGB <- function(tbl) {
    #  Rename and select columns of interest
    tbl %>%
        plyr::rename(
            .,
            c(
                "rname" = "GB_rname",
                "strand" = "GB_strand",
                "pos" = "GB_pos",
                "mrnm" = "GB_mrnm",
                "mpos" = "GB_mpos"
            )
        ) %>% 
        dplyr::select(
            coordinate,
            GB_rname,
            GB_strand,
            GB_pos,
            GB_mrnm,
            GB_mpos
        )
}


#  Set up work directories (locations to be changed) --------------------------
path.1 <- "/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora"
path.2 <- "/segregatedReads.SNPTHRESH1.Q30"


#  Load tibbles (already sorted) from .rds files ------------------------------

#  KA assignments ---------------------
path.1 %>% setwd()
# chromosome <- "chr1"
chromosome <- "chrX"
search <- paste0("\\.", chromosome, ".rds$")
file <- c(list.files(path = ".", pattern = search))
variable <- file %>% gsub("\\.rds$", "", .)
loadTibbleFromRDS(variable = variable, file = file)
#NOTE Remember, "alt" is "CAST", "ref" is "129"

variable_KA <- c(
    "KA.129S1",
    "KA.CAST",
    "KA.NA",
    "KA.Ambiguous"
)
command <- paste0("<- ", variable)
operation <- makeOperation(variable_KA, command)
eval(parse(text = operation))

command <- paste0("rm(", variable, ")")
eval(parse(text = command))

#  GB assignments ---------------------
paste0(path.1, path.2) %>% setwd()
search <- paste0("^(dedup.joint.*", chromosome, ".rds$)")
file <- c(list.files(path = ".", pattern = search))
variable <- file %>% gsub("\\.rds$", "", .)
loadTibbleFromRDS(variable = variable, file = file)

variable_GB <- c(
    "GB.alt.CAST",
    "GB.ambig",
    "GB.contra",
    "GB.ref.129S1"
)
command <- paste0("<- ", variable)
operation <- makeOperation(variable_GB, command)
eval(parse(text = operation))

command <- paste0("rm(", variable, ")")
eval(parse(text = command))

rm(variable)


#  Identify intersections between "assignment" and "dedup.joint" --------------
for (i in c(1:4)) {
# for (i in c(4)) {
    command <- paste0(
        "<- ifelse(", variable_KA[i], "$coordinate %in% ", variable_GB[4], "$coordinate, '1', '0')"
    )
    operation <- makeOperation(paste0(variable_KA[i], "$in_ref_129S1"), command)
    print(operation)
    eval(parse(text = operation))
    
    command <- paste0(
        "<- ifelse(", variable_KA[i], "$coordinate %in% ", variable_GB[1], "$coordinate, '1', '0')"
    )
    operation <- makeOperation(paste0(variable_KA[i], "$in_alt_CAST"), command)
    print(operation)
    eval(parse(text = operation))
    
    command <- paste0(
        "<- ifelse(", variable_KA[i], "$coordinate %in% ", variable_GB[2], "$coordinate, '1', '0')"
    )
    operation <- makeOperation(paste0(variable_KA[i], "$in_ambig"), command)
    print(operation)
    eval(parse(text = operation))
    
    command <- paste0(
        "<- ifelse(", variable_KA[i], "$coordinate %in% ", variable_GB[3], "$coordinate, '1', '0')"
    )
    operation <- makeOperation(paste0(variable_KA[i], "$in_contra"), command)
    print(operation)
    eval(parse(text = operation))
    
    command <- paste0(
        "<- paste0(",
        variable_KA[i], "$in_ref_129S1, ",
        variable_KA[i], "$in_alt_CAST, ",
        variable_KA[i], "$in_ambig, ",
        variable_KA[i], "$in_contra)"
    )
    operation <- makeOperation(paste0(variable_KA[i], "$GB_intersection"), command)
    print(operation)
    eval(parse(text = operation))
    
    command <- paste0("<- ", variable_KA[i], "$GB_intersection %>% as.factor()")
    operation <- makeOperation(paste0(variable_KA[i], "$GB_intersection"), command)
    print(operation)
    eval(parse(text = operation))
}

#  Evaluate intersection tables
# for (i in c(1:4)) {
#     cat("# -----------------\n")
#     
#     command <- paste0(variable_KA[i], "$GB_intersection %>% table()")
#     print(command)
#     eval(parse(text = command)) %>% print()
#     
#     command <- paste0(variable_KA[i], "$GB_intersection %>% table() %>% prop.table()")
#     eval(parse(text = command)) %>% print()
#     
#     cat("\n")
# }


#  Generate table -------------------------------------------------------------
command <- paste0(
    "<- table(", variable_KA, "$GB_intersection) %>% ",
    "as.data.frame() %>% ",
    "as_tibble() %>%",
    "dplyr::rename(c(GB.assignment = Var1, ", variable_KA, " = Freq))"
)
operation <- makeOperation(paste0(variable_KA, ".table"), command)
eval(parse(text = operation))

KA.all.table <- full_join(KA.129S1.table, KA.CAST.table, by = "GB.assignment") %>%
    full_join(., KA.Ambiguous.table, by = "GB.assignment") %>%
    full_join(., KA.NA.table, by = "GB.assignment") %>%
    purrr::map_df(rev)


#  Generate table of proportions ----------------------------------------------
command <- paste0(
    "<- ", variable_KA, "$GB_intersection %>% ",
    "table() %>% ",
    "prop.table() %>% ",
    "as.data.frame() %>% ",
    "as_tibble() %>% ",
    "dplyr::rename(c(GB.assignment = ., ", variable_KA, " = Freq))"
)
operation <- makeOperation(paste0(variable_KA, ".table.prop"), command)
eval(parse(text = operation))

table.prop <- full_join(KA.129S1.table.prop, KA.CAST.table.prop, by = "GB.assignment") %>%
    full_join(., KA.Ambiguous.table.prop, by = "GB.assignment") %>%
    full_join(., KA.NA.table.prop, by = "GB.assignment") %>%
    purrr::map_df(rev)


#  Generate table of percentages ----------------------------------------------
command <- paste0("<- ", variable_KA, ".table.prop")
operation <- makeOperation(paste0(variable_KA, ".table.percent"), command)
eval(parse(text = operation))

command <- paste0(
    "<- ", variable_KA, ".table.percent$", variable_KA, " %>% convertPercent()"
)
operation <- makeOperation(paste0(variable_KA, ".table.percent$", variable_KA), command)
eval(parse(text = operation))

table.percent <- full_join(KA.129S1.table.percent, KA.CAST.table.percent, by = "GB.assignment") %>%
    full_join(., KA.Ambiguous.table.percent, by = "GB.assignment") %>%
    full_join(., KA.NA.table.percent, by = "GB.assignment") %>%
    purrr::map_df(rev)


#  Rename GB.assignment factor levels -----------------------------------------
variable_table <- c(
    "KA.all.table",
    "table.prop",
    "table.percent"
)

command <- paste0(
    "<- plyr::revalue(", variable_table, "$GB.assignment, ",
    "c('0000' = 'GB.Not_present', '0001' = 'GB.Contra', '0010' = 'GB.Ambiguous', '0100' = 'GB.CAST', '1000' = 'GB.129S1'))"
)
operation <- makeOperation(paste0(variable_table, "$GB.assignment"), command)
eval(parse(text = operation))


#  Create a full table of counts, including row and column sums ---------------
tmp.1 <- KA.all.table %>% summarise(across(2:5, ~sum(.)))

tmp.2 <- KA.all.table %>%
    select(colnames(.)[2:length(colnames(.))]) %>% 
    rowSums() %>%
    as_tibble_col()

tmp.3 <- tmp.2 %>%
    summarise(across(1, ~sum(.)))

tmp.4 <- bind_rows(tmp.2, tmp.3)

table.full <- bind_rows(KA.all.table, tmp.1)
table.full <- bind_cols(table.full, tmp.4)
table.full$GB.assignment <- table.full$GB.assignment %>% as.character()
table.full[6, 1] <- "GB.sum"
table.full$GB.assignment <- table.full$GB.assignment %>% as_factor()
table.full <- table.full %>% dplyr::rename(., sum.row = value)
colnames(table.full)[6] <- "KA.sum"

rm(tmp.1, tmp.2, tmp.3, tmp.4)


#  Create a full table of proportions -----------------------------------------
tmp.1 <- table.full$KA.sum[-6] %>% as_tibble_col()
tmp.2 <- tmp.1 / colSums(tmp.1)

table.prop.full <- cbind(table.prop, tmp.2) %>% as_tibble()
colnames(table.prop.full)[6] <- "KA.sum"

rm(tmp.1, tmp.2)


#  Create a full table of percentages -----------------------------------------
column <- c(
    "KA.129S1",
    "KA.CAST",
    "KA.Ambiguous",
    "KA.NA",
    "KA.sum"
)

table.percent.full <- list()
for (i in column) {
    command <- paste0("<- table.prop.full$", i," %>% convertPercent()")
    operation <- makeOperation(paste0("table.percent.full$", i), command)
    print(operation)
    eval(parse(text = operation))
}
table.percent.full <- table.percent.full %>%
    as_tibble() %>%
    purrr::map_df(rev)

table.percent.full <- cbind(table.percent[, 1], table.percent.full) %>%
    as_tibble()

table.percent.full$GB.assignment <- table.percent.full$GB.assignment %>%
    as.factor()


#  Print the full tables ------------------------------------------------------
table.full
table.full.df <- data.frame(table.full[, 2:5], row.names = table.full$GB.assignment)
table.full.df <- table.full.df[-6, ]

breaks <- seq(0, 50000, 100)
pheatmap::pheatmap(
    table.full.df,
    main = paste0("KA × GB assignments, ", chromosome),
    display_numbers = TRUE,
    number_format = "%.0f",
    color = colorRampPalette(c('white','red'))(length(breaks)),
    breaks = breaks,
    border_color = "#FFFFFF",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_number = 11
)

table.prop.full
table.prop.full.df <- data.frame(table.prop.full[, 2:5], row.names = table.prop.full$GB.assignment)
pheatmap::pheatmap(
    table.prop.full.df %>% as_tibble(),
    main = paste0("KA × GB assignments, ", chromosome),
    labels_row = table.prop.full$GB.assignment,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = colorRampPalette(c('white','red'))(100),
    border_color = "#FFFFFF",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_number = 11
)


#  How effective was UCSC liftOver? -------------------------------------------
setwd(path.1)

liftOver.129S1 <- loadRDS(paste0("dedup.129S1.", chromosome, ".liftOvers.rds"))
liftOver.129S1$assembly <- "129S1-SvImJ"

liftOver.CAST <- loadRDS(paste0("dedup.CAST.", chromosome, ".liftOvers.rds"))
liftOver.CAST$assembly <- "CAST-EiJ"

liftOver.129S1$liftOver_reason %>% levels()

#  Reorder liftOver_reason factor levels
variable_reason <- paste0("liftOver.", c("129S1", "CAST"), "$liftOver_reason")
order_reason <- c(
    "liftOver successful",
    "Deleted in new",
    "Partially deleted in new",
    "Split in new"
)
command <- paste0("<- forcats::fct_relevel(", variable_reason, ", order_reason)")
operation <- makeOperation(variable_reason, command)
eval(parse(text = operation))

#  Status of liftOver
#DONE Make sure the levels are ordered the same way...
liftOver.129S1 %>% ggplot(., aes(x = assembly, color = liftOver_reason, fill = liftOver_reason)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ggtitle("liftOver status:\n129S1-SvImJ to mm10") +
    ylab("") +
    xlab(chromosome) +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

liftOver.CAST %>% ggplot(., aes(x = assembly, color = liftOver_reason, fill = liftOver_reason)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ggtitle("liftOver status:\nCAST-EiJ to mm10") +
    ylab("") +
    xlab(chromosome) +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)


#  How many chromosomes changed in liftOver? ----------------------------------
variable <- c("liftOver.129S1", "liftOver.CAST")
command <- paste0("<- (", variable, "$old_rname == ", variable, "$rname)")
operation <- makeOperation(paste0(variable, "$liftOver_rname_equal"), command)
eval(parse(text = operation))

command <- paste0("<- ", variable, " %>% dplyr::relocate(liftOver_rname_equal, .after = liftOver_reason)")
operation <- makeOperation(variable, command)
eval(parse(text = operation))

#TODO Skip these analyses for now b/c need to go back and store liftOver rname,
#+ pos, etc. info for both 129S1 and CAST; right now, only have info for 129S1
# #  Neutral w/r/t 129S1
# variable <- c("tmp.final.Ambiguous")
# command <- paste0("<- (", variable, "$old_rname_129S1 == ", variable, "$rname)")
# operation <- makeOperation(paste0(variable, "$liftOver_rname_129S1_equal"), command)
# eval(parse(text = operation))
# 
# command <- paste0("<- ", variable, " %>% dplyr::relocate(liftOver_rname_129S1_equal, .after = liftOver_reason)")
# operation <- makeOperation(variable, command)
# eval(parse(text = operation))
# 
# #  Neutral w/r/t CAST
# command <- paste0("<- (", variable, "$old_rname_CAST == ", variable, "$rname)")
# operation <- makeOperation(paste0(variable, "$liftOver_rname_CAST_equal"), command)
# eval(parse(text = operation))
# 
# command <- paste0("<- ", variable, " %>% dplyr::relocate(liftOver_rname_CAST_equal, .after = liftOver_rname_129S1_equal)")
# operation <- makeOperation(variable, command)
# eval(parse(text = operation))

#  Post-liftOver on the same chromosome? TRUE, FALSE, or NA
#NOTE No need to show this
# liftOver.129S1 %>% ggplot(., aes(x = liftOver_rname_equal, color = liftOver_reason, fill = liftOver_reason)) +
#     geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
#     ggtitle("liftOver status:\n129S1-SvImJ to mm10") +
#     ylab("") +
#     xlab(chromosome) +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma)
# 
# liftOver.CAST %>% ggplot(., aes(x = liftOver_rname_equal, color = liftOver_reason, fill = liftOver_reason)) +
#     geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
#     ggtitle("liftOver status:\nCAST-EiJ to mm10") +
#     ylab("") +
#     xlab(chromosome) +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma)

#  Post-liftOver rnames per assembly
#NOTE No need to show this
# liftOver.129S1 %>% ggplot(., aes(x = assembly, color = rname, fill = rname)) +
#     geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
#     ggtitle("liftOver status:\n129S1-SvImJ to mm10") +
#     ylab("") +
#     xlab(chromosome) +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma)
# 
# liftOver.CAST %>% ggplot(., aes(x = assembly, color = rname, fill = rname)) +
#     geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
#     ggtitle("liftOver status:\nCAST-EiJ to mm10") +
#     ylab("") +
#     xlab(chromosome) +
#     theme(legend.title = element_blank()) +
#     scale_y_continuous(labels = comma)


#  Examine what alignments were given what assignments ------------------------
setwd(path.1)
joint <- bind_rows(KA.129S1, KA.CAST, KA.Ambiguous, KA.NA)

joint$assignment[is.na(joint$assignment)] <- "NA"
joint$assignment <- joint$assignment %>%
    as_factor() %>%
    plyr::revalue(
        c("Neutral" = "Ambiguous")
    )

#  Reorder trinary factor levels
order_trinary <- c(
    "010", "001", "011", "100", "110", "101", "111"
)
joint$trinary <- forcats::fct_relevel(joint$trinary, order_trinary)

joint$trinary <- joint$trinary %>%
    plyr::revalue(
        c(
            "010" = "129S1-SvImJ",
            "001" = "CAST-EiJ",
            "011" = "129S1-SvImJ, CAST-EiJ",
            "100" = "mm10-N-masked",
            "110" = "mm10-N-masked, 129S1-SvImJ",
            "101" = "mm10-N-masked, CAST-EiJ",
            "111" = "All"
        )
    )

#  What are the numbers of assignments?
ggplot(joint, aes(x = assignment)) +
    geom_bar(alpha = 0.5) +
    geom_text(stat = "count", aes(label = scales::comma(..count..)), size = 3, vjust = -0.5) +
    ggtitle("What are the numbers of assignments?") +
    ylab("") +
    xlab(paste0("assignment, ", chromosome)) +
    scale_y_continuous(labels = comma)

#  What alignments comprise each assignment?
#DONE Change the order of the levels
ggplot(joint, aes(x = assignment, color = trinary, fill = trinary)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ggtitle("What alignments comprise each assignment?") +
    ylab("") +
    xlab(paste0("assignment, ", chromosome)) +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)


#TODO Figure out where this needs to go
# liftOver.129S1$GB_intersection <- liftOver.129S1$GB_intersection %>%
#     plyr::revalue(
#         c(
#             '0000' = 'GB.Not_present',
#             '0001' = 'GB.Contra',
#             '0010' = 'GB.Ambiguous',
#             '0100' = 'GB.CAST',
#             '1000' = 'GB.129S1'
#         )
#     )
# 
# liftOver.CAST$GB_intersection <- liftOver.CAST$GB_intersection %>%
#     plyr::revalue(
#         c(
#             '0000' = 'GB.Not_present',
#             '0001' = 'GB.Contra',
#             '0010' = 'GB.Ambiguous',
#             '0100' = 'GB.CAST',
#             '1000' = 'GB.129S1'
#         )
#     )
