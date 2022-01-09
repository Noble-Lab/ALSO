#  Scraps ---------------------------------------------------------------------
output.1 <- vector(mode = "logical", length = length(pos))
for (i in 1:length(pos)) {
    output.1[i] <- pos[i] == mpos[i + 1] & mpos[i] == pos[i + 1]
}
output.1 %>% table()  # 54081 - 54056 = 25

output.2 <- vector(mode = "logical", length = length(output.1))
for (i in 1:length(output.1)) {
    output.2[i] <- (output.1[i] == output.1[i + 1])
}
output.2 %>% table()

test.2$discrepany <- output.2

# -------
v1 <- which(test.2$discrepany == TRUE)
v2 <- v1 - 1
v3 <- v1 + 1
v4 <- c(v1, v2, v3) %>% sort() %>% unique()
rm(v1, v2, v3)

test.3 <- test.2[v4, ]
test.3 <- test.3 %>% dplyr::relocate(discrepany, .before = qname)

test.3$isize_abs <- test.3$isize %>% abs()
test.3 <- test.3[order(test.3$isize_abs, test.3$pos), ]
test.3 <- test.3 %>% dplyr::relocate(isize_abs, .after = mpos)

# --
pos <- test.3$qpos
mpos <- test.3$qmpos

output.1 <- vector(mode = "logical", length = length(pos))
for (i in 1:length(pos)) {
    output.1[i] <- pos[i] == mpos[i + 1] & mpos[i] == pos[i + 1]
}
output.1 %>% table()  # 54081 - 54056 = 25
# .
# FALSE  TRUE 
#    32    33

output.2 <- vector(mode = "logical", length = length(output.1))
for (i in 1:length(output.1)) {
    output.2[i] <- (output.1[i] == output.1[i + 1])
}
output.2 %>% table()
# .
# FALSE 
#    64
rm(pos, mpos, output.1, output.2, i)


# -----------------
test.2$rowid <- 1:nrow(test.2)
test.4 <- test.2[test.2$rowid %notin% v4, ]
# nrow(test.2) == nrow(test.3) + nrow(test.4)  # [1] TRUE

pos <- test.4$qpos
mpos <- test.4$qmpos

pos.in.mpos <- pos %in% mpos
mpos.in.pos <- mpos %in% pos

pos.in.mpos %>% table()
mpos.in.pos %>% table()


# -----------------------------------------------------------------------------
# -----------------
# unique(pos) %>% length()
# unique(mpos) %>% length()
# 
# unique(test$qname) %>% length()

match.pos.mpos <- match(pos, mpos)
output <- vector(mode = "logical", length = length(match.pos.mpos))
for (i in 1:length(match.pos.mpos)) {
    output[i] <- ((match.pos.mpos[i] + 1) == match.pos.mpos[i + 1])
}
output %>% table()

# test.2 <- test[match.pos.mpos, ]

strand <- test.2$strand
output <- vector(mode = "logical", length = length(strand))
for (i in 1:length(strand)) {
    output[i] <- (strand[i] == strand[i + 1])
}
output %>% table()
test.2$disruption <- output
rm(strand, output, i)

# test.2 <- test.2 %>% dplyr::relocate(mpos, .after = pos)
# test.3 <- test.2[test.2$disruption == TRUE, ]

# -----------------
strand <- test$strand
output <- vector(mode = "logical", length = length(strand))
for (i in 1:length(strand)) {
    output[i] <- (strand[i] == strand[i + 1])
}
output %>% table()
test$disruption <- output
rm(strand, output, i)

test <- test %>% dplyr::relocate(mpos, .after = pos)
test.2 <- test[test$disruption == TRUE, ]
test.3 <- test[test$disruption == FALSE, ]


# -----------------
pos <- test.3$pos
mpos <- test.3$mpos
output.1 <- vector(mode = "logical", length = length(pos))
for (i in 1:length(pos)) {
    output.1[i] <- (pos[i] == mpos[i + 1])
}
output.1 %>% table()

output.2 <- vector(mode = "logical", length = length(output.1))
for (i in 1:length(output.1)) {
    output.2[i] <- (output.1[i] == output.1[i + 1])
}
output.2 %>% table()

test.3$disruption.2 <- output.2
rm(output.1, output.2, i)

test.4 <- test.3[test.3$disruption.2 == TRUE, ]
test.5 <- test.3[test.3$disruption.2 == FALSE, ]

# -----------------
pos <- test.5$pos
mpos <- test.5$mpos
output.1 <- vector(mode = "logical", length = length(pos))
for (i in 1:length(pos)) {
    output.1[i] <- (pos[i] == mpos[i + 1])
}
output.1 %>% table()

output.2 <- vector(mode = "logical", length = length(output.1))
for (i in 1:length(output.1)) {
    output.2[i] <- (output.1[i] == output.1[i + 1])
}
output.2 %>% table()


# -----------------
View(test)
rm(test)


#  Order columns by coordinate, tag.AS, and mapq ------------------------------
command <- paste0(
    "<- ", variable, "[order(",
    variable, "$qname, ",
    variable, "$rname, ",
    variable, "$pos, ",
    variable, "$mpos, -",
    variable, "$AS, -",
    variable, "$mapq), ]"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))


#  Create a "temporary tag" for use in distinguishing entries
#+ 
#+ This is necessary to filter out duplicates and join liftOver data to initial
#+ tibbles; it's named "old_coordinate", i.e., the coordinate prior to liftOver
command <- paste0(
    command_pipe,
    "tidyr::unite(",
    "old_coordinate, ",
    "c(\"qname\", \"rname\", \"pos\", \"pos_end\"), ",
    # "c(\"qname\", \"rname\", \"pos\", \"mpos\"), ",  #TODO Think about it...
    "sep = \"_\", ",
    "remove = FALSE",
    ")"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))


#  Filter out all rows without unique coordinates -----------------------------
command <- paste0(
    command_pipe, "dplyr::distinct(old_coordinate, .keep_all = TRUE)"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))

o.ambiguous$pos[o.ambiguous$pos %in% o.ambiguous$mpos] %>% length()  # [1] 108155
o.ambiguous$pos[o.ambiguous$mpos %in% o.ambiguous$pos] %>% length()  # [1] 108184

all(o.ambiguous$pos %in% o.ambiguous$mpos)  # FALSE
all(o.ambiguous$mpos %in% o.ambiguous$pos)  # FALSE
#  This means that some lines don't have mates, i.e., there are $pos associated
#+ with $mpos that aren't present, and there are $mpos associated with $pos
#+ that aren't present

o.ambiguous.bak <- o.ambiguous
o.ambiguous <- o.ambiguous.bak

o.ambiguous <- o.ambiguous %>%
    select(old_coordinate, strand, pos, mpos, isize)

o.ambiguous <- o.ambiguous[order(o.ambiguous$pos, o.ambiguous$mpos), ]


#  Test -----------------------------------------------------------------------
index <- match(o.ambiguous$pos, o.ambiguous$mpos)
# index %>% duplicated() %>% table()
# index %>% is.na() %>% table()

which(index %>% duplicated())
which(index %>% is.na())

test <- o.ambiguous[index, ]
test$index <- index

test <- test[!duplicated(test), ]
test <- test[!is.na(test$index), ]
# test <- test[test$pos != test$mpos, ]

test$abs_isize <- abs(test$isize)

test <- test[!test$pos %notin% test$mpos, ]
test <- test[!test$mpos %notin% test$pos, ]

test[test$pos %notin% test$mpos, ]
test[test$mpos %notin% test$pos, ]

test.2 <- test[order(-test$abs_isize), ]
test.2$strand
test.2$strand %>% tail(n = 1000)

index <- match(test$pos, test$mpos)
test[index, ]
test$index <- index

test$strand

is.na(test$index) %>% table()
which(test$index %>% is.na())
test[is.na(test$index), ]

duplicated(test$index) %>% table()
x1 <- which(test$index %>% duplicated())
x2 <- (x1 - 1)
x3 <- (x1 + 1)
x4 <- (x1 - 2)
x5 <- (x1 + 2)
y <- c(x4, x2, x1, x3, x5) %>% sort()

test.2 <- test[y, ] %>% select(-qname,)
rm(x1, x2, x3, x4, x5, y, test.2)


# -----------------
test.6 <- dplyr::bind_rows(test.3, test.5)
test.6$qpos_mpos <- paste0(
    test.6$qname, "_",
    test.6$pos, "_",
    test.6$mpos
)
test.6 <- test.6 %>% dplyr::relocate(qpos_mpos, .before = discrepancy)

output.1 <- testMatesPaired(pos = test.6$qpos, mpos = test.6$qmpos)
output.2 <- testLogicalAlternating(logical = output.1)
output.2 %>% table()

test.6$discrepancy <- output.2

discrepancy <- createVectorDiscrepancy(discrepancy = test.6$discrepancy)

d1 <- test.6[discrepancy, ]
d2 <- test.6[duplicated(test.6$qpos_mpos), ]

# plus_minus <- createVectorPlusMinus(which(test.6$qpos_mpos %>% duplicated()))
# d.3 <- test.6[plus_minus, ]
# View(d.3)


# -----------------
test.7 <- test.6[!duplicated(test.6$qpos_mpos), ]

output.1 <- testMatesPaired(pos = test.7$qpos, mpos = test.7$qmpos)
output.2 <- testLogicalAlternating(logical = output.1)
output.2 %>% table()

test.7$discrepancy <- output.2

discrepancy <- createVectorDiscrepancy(discrepancy = test.7$discrepancy)

d1 <- test.7[discrepancy, ]
d2 <- test.7[duplicated(test.7$qpos_mpos), ]
