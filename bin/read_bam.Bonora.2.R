#!/usr/bin/env Rscript

#  Load libraries
library(dplyr)
library(forcats)
library(ggplot2)
library(magrittr)
library(parallel)
library(purrr)
library(Rsamtools)
library(scales)

options(pillar.sigfig = 8, scipen = 10000)


# -----------------------------------------------------------------------------
#  Temporary: Set up work directory (location TB∆)
"/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora" %>% setwd()
# getwd()
# list.files()
# list.dirs()


# -----------------------------------------------------------------------------
makeOperation <- function(variable, command) {
    operation <- paste(variable, command) %>% paste(., collapse = "; ")
    return(operation)
}

#  What to query from .bam files
map_info <- c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar", "mrnm", "mpos", "isize", "seq", "qual")
tag_info <- c("AS", "XS", "NM")
map_params <- Rsamtools::ScanBamParam(what = map_info, tag = tag_info)

#  Assign variables necessary for subsequent commands
# files <- list.files(pattern = "\\.dedup.bam$")
files <- list.files(pattern = "\\.MarkDuplicates.sort.bam$")
n <- length(files)
variables <- paste0("bam", 1:n)
variables_rname <- paste0("bam", 1:n, "$rname")
command_pipe <- paste0("<- bam", 1:n, " %>% ")
chromosomes <- c(paste0("chr", 1:19), "chrX")

mapply(
    assign, variables, files, MoreArgs = list(envir = parent.frame())
)

#  Assign .bam information as list
dedup.129 <- bam1 %>% Rsamtools::scanBam(., param = map_params)
dedup.CAST <- bam2 %>% Rsamtools::scanBam(., param = map_params)
predup.CAST <- bam3 %>% Rsamtools::scanBam(., param = map_params)
dedup.mm10_CAST_129_Nmasked <- bam4 %>% Rsamtools::scanBam(., param = map_params)
dedup.mm10 <- bam5 %>% Rsamtools::scanBam(., param = map_params)

#  Convert .bam information from list to dataframe to tibble
variables <- c(
    "dedup.129",
    "dedup.CAST",
    "predup.CAST",
    "dedup.mm10_CAST_129_Nmasked",
    "dedup.mm10"
)

command_pipe <- paste0("<- ", variables, " %>% ")

command <- paste0(command_pipe, "as.data.frame() %>% as_tibble()")
operation <- makeOperation(variables, command)
eval(parse(text = operation))
rm("bam1", "bam2", "bam3", "bam4", "bam5")

#  Reorder rname factor levels
variables_rname <- paste0(variables, "$rname")
command <- paste0("<- forcats::fct_relevel(", variables_rname, ", chromosomes)")
operation <- makeOperation(variables_rname, command)
eval(parse(text = operation))

#  Drop unused rname factor levels
command <- paste0("<- ", variables_rname, " %>% forcats::fct_drop()")
operation <- makeOperation(variables_rname, command)
eval(parse(text = operation))

#  Drop rows that are not chr1-19, chrX
command <- paste0(command_pipe, "filter(., rname %in% chromosomes)")
operation <- makeOperation(variables, command)
eval(parse(text = operation))

dedup.129 %>% head()
dedup.CAST %>% head()
predup.CAST %>% head()
dedup.mm10_CAST_129_Nmasked %>% head()
dedup.mm10 %>% head()


# -----------------------------------------------------------------------------
b_f.dedup.129 <- dedup.129 %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    ) %>%
    dplyr::arrange(rname, pos, tag.AS) %>% 
    dplyr::distinct(qname, seq, .keep_all = TRUE)

b_f.dedup.CAST <- dedup.CAST %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    ) %>%
    dplyr::arrange(rname, pos, tag.AS) %>% 
    dplyr::distinct(qname, seq, .keep_all = TRUE)

b_f.dedup.mm10_CAST_129_Nmasked <- dedup.mm10_CAST_129_Nmasked %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    ) %>%
    dplyr::arrange(rname, pos, tag.AS) %>% 
    dplyr::distinct(qname, seq, .keep_all = TRUE)

b_f.dedup.mm10 <- dedup.mm10 %>%
    tidyr::unite(
        b_f,
        c("qname", "seq"),
        sep = "_",
        remove = FALSE
    ) %>%
    dplyr::arrange(rname, pos, tag.AS) %>% 
    dplyr::distinct(qname, seq, .keep_all = TRUE)


# -----------------------------------------------------------------------------
d.1.b_f <- b_f.dedup.129$b_f
d.1.b_f %>% length()  # [1] 7675828
d.1.b_f %>% unique() %>% length()  # [1] 7675828

d.C.b_f <- b_f.dedup.CAST$b_f
d.C.b_f %>% length()  # [1] 7616363
d.C.b_f %>% unique() %>% length()  # [1] 7616363

d.Mn.b_f <- b_f.dedup.mm10_CAST_129_Nmasked$b_f
d.Mn.b_f %>% length()  # [1] 7574862
d.Mn.b_f %>% unique() %>% length()  # [1] 7574862

d.M.b_f <- b_f.dedup.mm10$b_f
d.M.b_f %>% length()  # [1] 7634935
d.M.b_f %>% unique() %>% length()  # [1] 7634935

# -----------------
sd.1.C <- setdiff(d.1.b_f, d.C.b_f)
sd.1.Mn <- setdiff(d.1.b_f, d.Mn.b_f)
sd.C.1 <- setdiff(d.C.b_f, d.1.b_f)
sd.C.Mn <- setdiff(d.C.b_f, d.Mn.b_f)
sd.Mn.C <- setdiff(d.Mn.b_f, d.C.b_f)
sd.Mn.1 <- setdiff(d.Mn.b_f, d.1.b_f)

sd.1.C %>% head() %>% as_tibble()
sd.1.C %>% length()  # [1] 730931
sd.1.C %>% unique() %>% length()  # [1] 730931

sd.1.Mn %>% head() %>% as_tibble()
sd.1.Mn %>% length()  # [1] 654142
sd.1.Mn %>% unique() %>% length()  # [1] 654142

sd.C.1 %>% head() %>% as_tibble()
sd.C.1 %>% length()  # [1] 671466
sd.C.1 %>% unique() %>% length()  # [1] 671466

sd.C.Mn %>% head() %>% as_tibble()
sd.C.Mn %>% length()  # [1] 707396
sd.C.Mn %>% unique() %>% length()  # [1] 707396

sd.Mn.C %>% head() %>% as_tibble()
sd.Mn.C %>% length()  # [1] 665895
sd.Mn.C %>% unique() %>% length()  # [1] 665895

sd.Mn.1 %>% head() %>% as_tibble()
sd.Mn.1 %>% length()  # [1] 553176
sd.Mn.1 %>% unique() %>% length()  # [1] 553176

# -----------------
`%notin%` <- Negate(`%in%`)
check_in.129 <- b_f.dedup.129[d.1.b_f %in% d.C.b_f, ]  # 6944897
check_notin.129 <- b_f.dedup.129[d.1.b_f %notin% d.C.b_f, ]  # 730931
6944897 + 730931  # [1] 7675828

check_in.CAST <- b_f.dedup.CAST[d.C.b_f %in% d.1.b_f, ]  # 6944897
check_notin.CAST <- b_f.dedup.CAST[d.C.b_f %notin% d.1.b_f, ]  # 671456
6944897 + 671456  # [1] 7616353


# -----------------------------------------------------------------------------
#  https://stat545.com/join-cheatsheet.html#full_joinsuperheroes-publishers

z1 <- b_f.dedup.mm10 %>% select(b_f, tag.AS)
zC <- b_f.dedup.CAST %>% select(b_f, tag.AS)
zMn <- b_f.dedup.mm10_CAST_129_Nmasked %>% select(b_f, tag.AS)

z1 <- z1 %>%
    dplyr::rename(
        oneTwentyNine.AS = tag.AS,
        # mm10.XS = tag.XS
    )
zC <- zC %>%
    dplyr::rename(
        CAST.AS = tag.AS,
        # CAST.XS = tag.XS
    )
zMn <- zMn %>%
    dplyr::rename(
        mm10_CAST_129_Nmasked.AS = tag.AS,
        # mm10.XS = tag.XS
    )

z1.zC <- full_join(z1, zC, by = "b_f")
z1.zMn <- full_join(z1, zMn, by = "b_f")
zC.zMn <- full_join(zC, zMn, by = "b_f")

z1.zC.zMn <- full_join(z1, zC, by = "b_f") %>%
    full_join(., zMn, by = "b_f")
zMn.z1.zC <- full_join(zMn, z1, by = "b_f") %>% 
    full_join(., zC, by = "b_f")

# -----------------
test1 <- zMn.z1.zC %>% 
    mutate(difference = oneTwentyNine.AS - CAST.AS)

int <- 0 %>% as.integer()  #TODO Make the integer an argument

# -----------------
test2 <- test1 %>% 
    dplyr::mutate(
        assignment.initial = case_when(
            difference >= (-1 * int) & difference <= int ~ "Neutral",
            difference > int ~ "129",
            difference < (-1 * int) ~ "CAST-EiJ"
        )
    )

test2$assignment <- ifelse(
    is.na(test2$assignment.initial),
    ifelse(
        !is.na(test2$oneTwentyNine.AS),
        "129",
        ifelse(
            !is.na(test2$CAST.AS),
            "CAST-EiJ",
            test2$assignment.initial
        )
    ),
    test2$assignment.initial
) %>% forcats::as_factor()

# -----------------
test3 <- test2 %>% 
    select(b_f, mm10_CAST_129_Nmasked.AS, oneTwentyNine.AS, CAST.AS, assignment)

test3$tmp.mm10_CAST_129_Nmasked <- ifelse(
    is.na(test3$mm10_CAST_129_Nmasked.AS), "0", "1"
)

test3$tmp.129 <- ifelse(
    is.na(test3$oneTwentyNine.AS), "0", "1"
)

test3$tmp.CAST <- ifelse(
    is.na(test3$CAST.AS), "0", "1"
)

test3$trinary <- paste0(test3$tmp.mm10_CAST_129_Nmasked, test3$tmp.129, test3$tmp.CAST) %>%
    forcats::as_factor()

test3$assignment_trinary <- paste(test3$assignment, test3$trinary) %>% as_factor()

# -----------------
test4 <- test3 %>%
    select(
        b_f,
        mm10_CAST_129_Nmasked.AS,
        oneTwentyNine.AS,
        CAST.AS,
        assignment,
        trinary,
        assignment_trinary
    )
# View(test4)


#  Make figures
# ------------------
ggplot(test4, aes(x = trinary)) +
    geom_bar(alpha = 0.5) +
    geom_text(stat = "count", aes(label = scales::comma(..count..)), size = 3, vjust = -0.5) +
    ylab("") +
    xlab("assembly-alignment combination") +
    scale_y_continuous(labels = comma)

ggplot(test4, aes(x = assignment)) +
    geom_bar(alpha = 0.5) +
    geom_text(stat = "count", aes(label = scales::comma(..count..)), size = 3, vjust = -0.5) +
    ylab("") +
    xlab("alignment assembly assignment") +
    scale_y_continuous(labels = comma)

ggplot(test4, aes(x = assignment_trinary)) +
    geom_bar(alpha = 0.5) +
    geom_text(stat = "count", aes(label = scales::comma(..count..)), size = 2.5, vjust = -0.5) +
    ylab("") +
    xlab("assembly-alignment assignment") +
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -----------------
ggplot(test4, aes(x = assignment, color = trinary, fill = trinary)) +
    geom_bar(alpha = 0.5, width = 0.9) +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

#  Cool
ggplot(test4, aes(x = assignment, color = trinary, fill = trinary)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

ggplot(test4, aes(x = trinary, color = assignment, fill = assignment)) +
    geom_bar(alpha = 0.5, width = 0.9) +
    ylab("") +
    xlab("assembly-alignment combination") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

#  Cool
ggplot(test4, aes(x = trinary, color = assignment, fill = assignment)) +
    geom_bar(alpha = 0.5, position = position_dodge2(preserve = "single"), width = 0.9) +
    ylab("") +
    xlab("assembly-alignment combination") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

# -----------------
# -------
ggplot(test4, aes(x = assignment, fill = as_factor(CAST.AS))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle("AS from alignment to CAST") +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

ggplot(test4, aes(x = assignment, fill = as_factor(CAST.AS))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle("AS from alignment to CAST") +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma) +
    facet_grid(trinary ~ .)

# -------
ggplot(test4, aes(x = assignment, fill = as_factor(oneTwentyNine.AS))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle("AS from alignment to 129") +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

ggplot(test4, aes(x = assignment, fill = as_factor(oneTwentyNine.AS))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle("AS from alignment to 129") +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma) +
    facet_grid(trinary ~ .)

# -------
ggplot(test4, aes(x = assignment, fill = as_factor(mm10_CAST_129_Nmasked.AS))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle("AS from alignment to mm10-CAST-129-N-masked") +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma)

ggplot(test4, aes(x = assignment, fill = as_factor(mm10_CAST_129_Nmasked.AS))) +
    geom_bar(alpha = 0.9, width = 0.9) +
    ggtitle("AS from alignment to mm10-CAST-129-N-masked") +
    ylab("") +
    xlab("assignment") +
    theme(legend.title = element_blank()) +
    scale_y_continuous(labels = comma) +
    facet_grid(trinary ~ .)

# -----------------
# ggplot(test4, aes(x = CAST.AS, color = assignment)) +
#     geom_histogram(binwidth = 1, fill = "white", alpha = 0.5)
# ggplot(test4, aes(x = oneTwentyNine.AS, color = assignment)) +
#     geom_histogram(binwidth = 1, fill = "white", alpha = 0.5)
# ggplot(test4, aes(x = mm10_CAST_129_Nmasked.AS, color = assignment)) +
#     geom_histogram(binwidth = 1, fill = "white", alpha = 0.5)

ggplot(test4, aes(x = CAST.AS, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to CAST")
ggplot(test4, aes(x = oneTwentyNine.AS, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to 129")
ggplot(test4, aes(x = mm10_CAST_129_Nmasked.AS, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to mm10-CAST-129-N-masked")

ggplot(test4, aes(x = CAST.AS, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to CAST") +
    facet_grid(assignment ~ .)
ggplot(test4, aes(x = oneTwentyNine.AS, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to 129") +
    facet_grid(assignment ~ .)
ggplot(test4, aes(x = mm10_CAST_129_Nmasked.AS, color = assignment, fill = assignment)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    ggtitle("AS from alignment to mm10-CAST-129-N-masked") +
    facet_grid(assignment ~ .)
