
#  Need to have loaded in dedup.129S1 from the following script:
#+ test.PE-processing.2021-1220.R

o.ambiguous <- dedup.129S1[dedup.129S1$mate_status == "ambiguous", ]
o.ambiguous$pos %>% duplicated() %>% table()
# .
#  FALSE   TRUE 
# 111462 148792

var <- "o.ambiguous"
command <- paste0("<- forcats::fct_relevel(", var, "$rname, chromosome)")
operation <- makeOperation(paste0(var, "$rname"), command)
eval(parse(text = operation))


#  Drop unused rname factor levels --------------------------------------------
command <- paste0("<- ", var, "$rname %>% forcats::fct_drop()")
operation <- makeOperation(paste0(var, "$rname"), command)
eval(parse(text = operation))


#  Drop rows that are not chr1-19, chrX ---------------------------------------
chromosomes <- c(paste0("chr", 1:19), "chrX")
command_pipe <- paste0("<- ", var, " %>% ")
command <- paste0(command_pipe, "filter(., rname %in% chromosomes)")
operation <- makeOperation(var, command)
eval(parse(text = operation))


#  Create and append pos_end, mpos_end columns --------------------------------
command <- paste0("<- ", var, "$pos + 49")
operation <- makeOperation(paste0(var, "$pos_end"), command)
eval(parse(text = operation))

command <- paste0("<- ", var, "$mpos + 49")
operation <- makeOperation(paste0(var, "$mpos_end"), command)
eval(parse(text = operation))

command <- paste0(
    command_pipe, "dplyr::relocate(pos_end, .after = pos) %>% dplyr::relocate(mpos_end, .after = mpos)"
)
operation <- makeOperation(var, command)
eval(parse(text = operation))

var <- "o.ambiguous"

#  Filter out rows with mapq less than 30 -------------------------------------
command <- paste0("<- ", var, "[", var, "$mapq >= 30, ]")
operation <- makeOperation(var, command)
eval(parse(text = operation))

#  Test: Deduplicating prior to sorting, etc. ---------------------------------
# o.ambiguous$pos_mpos <- paste0(o.ambiguous$pos, "_", o.ambiguous$mpos)
o.ambiguous$criteria <- paste0(
    o.ambiguous$qname, "_",
    o.ambiguous$flag, "_",
    o.ambiguous$pos, "_",
    o.ambiguous$mpos
)

# o.ambiguous$criteria %>% head(n = 100)
# o.ambiguous$criteria %>% duplicated() %>% head(n = 100)

# Create a back-up for the test var
o.ambiguous <- o.ambiguous.pre
o.ambiguous.pre <- o.ambiguous

# Get mpos close to pos
o.ambiguous <- o.ambiguous %>% dplyr::relocate(mpos, .after = pos)

#MAYBE
o.ambiguous$groupid <- o.ambiguous$groupid %>% as_factor()
# length_groupid <- o.ambiguous %>%
#     dplyr::group_by(groupid) %>%
#     summarize(no_rows = length(groupid))

o.ambiguous$qpos <- paste0(o.ambiguous$qname, "_", o.ambiguous$pos)
o.ambiguous$qmpos <- paste0(o.ambiguous$qname, "_", o.ambiguous$mpos)
o.ambiguous.same <- o.ambiguous[o.ambiguous$qpos == o.ambiguous$qmpos, ]
o.ambiguous <- o.ambiguous[!(o.ambiguous$qpos == o.ambiguous$qmpos), ]
#HERE - 1: all(o.ambiguous$qpos == ambiguous.dedup.129S1$qpos)  #[1] TRUE

test <- o.ambiguous[!duplicated(o.ambiguous$criteria), ]
# nrow(o.ambiguous) - nrow(test)  # [1] 140199

#HERE - 2: all(test$qpos == ambiguous.dedup.129S1$qpos)  #[1] TRUE
#HERE - 3: all(test$qpos == ambiguous.dedup.129S1$qpos)  #[1] TRUE

# -----------------
pos <- paste0(test$qname, "_", test$pos)
mpos <- paste0(test$qname, "_", test$mpos)

pos.in.mpos <- pos %in% mpos
mpos.in.pos <- mpos %in% pos

pos.in.mpos %>% table()
mpos.in.pos %>% table()
#HERE - 4: all(vPIM.ambiguous.dedup.129S1 == pos.in.mpos)  #[1] TRUE

# -------
test.2 <- test[pos.in.mpos, ]
#HERE - 5: all(test.ambiguous.dedup.129S1$qpos == test.2$qpos)  #[1] TRUE

pos <- paste0(test.2$qname, "_", test.2$pos)
mpos <- paste0(test.2$qname, "_", test.2$mpos)

pos.in.mpos <- pos %in% mpos
mpos.in.pos <- mpos %in% pos

pos.in.mpos %>% table()
# .
#   TRUE
# 108138
mpos.in.pos %>% table()
# .
#   TRUE
# 108138
#HERE - 6: tables are the same

test.2$isize_abs <- test.2$isize %>% abs()
test.2 <- test.2 %>% dplyr::relocate(isize_abs, .after = mpos)
#HERE - 6.5: all(test.2$qpos == test.ambiguous.dedup.129S1$qpos)  #[1] TRUE

# -------
test.3 <- test.2[order(test.2$isize_abs, test.2$qpos, test.2$qmpos), ]
#HERE - 7: all(test.3$qpos == test.ambiguous.dedup.129S1$qpos)  #[1] TRUE

pos <- test.3$qpos
mpos <- test.3$qmpos

output.1 <- vector(mode = "logical", length = length(pos))
for (i in 1:length(pos)) {
    output.1[i] <- pos[i] == mpos[i + 1] & mpos[i] == pos[i + 1]
}
output.1 %>% table()  # 54070 - 54067 = 3
# .
# FALSE  TRUE 
# 54070 54067
#HERE - 8: Tables same

output.2 <- vector(mode = "logical", length = length(output.1))
for (i in 1:length(output.1)) {
    output.2[i] <- (output.1[i] == output.1[i + 1])
}
output.2 %>% table()
# .
# FALSE   TRUE
# 108132     4
#HERE - 9: Tables same

test.3$discrepancy <- output.2
v1 <- which(test.3$discrepancy == TRUE)
v2 <- v1 - 1
v3 <- v1 + 1
discrepancy_plus_minus <- c(v1, v2, v3) %>% sort() %>% unique()
rm(v1, v2, v3)
#HERE - 10: all(discrepancy_plus_minus == vector_discrepancy)  #[1] TRUE


# -----------------
test.4 <- test.3[discrepancy_plus_minus, ]
test.4 <- test.4[order(test.4$isize_abs, test.4$qpos, test.4$qmpos, test.4$strand), ]
#HERE - 11: all(test.4$qpos == discrepancy.ambiguous.dedup.129S1$qpos)  #[1] TRUE


# -----------------
test.5 <- test.3[!duplicated(test.3$pos) & !duplicated(test.3$mpos), ]
#HERE - 12: all(test.5$qpos == test.ambiguous.dedup.129S1$qpos)  #[1] TRUE

pos <- test.5$qpos
mpos <- test.5$qmpos

output.1 <- vector(mode = "logical", length = length(pos))
for (i in 1:length(pos)) {
    output.1[i] <- pos[i] == mpos[i + 1] & mpos[i] == pos[i + 1]
}
output.1 %>% table()  # 53520 - 53521 = -1
# .
# FALSE  TRUE
# 53520 53521
#HERE - 13: OK; tables the same

output.2 <- vector(mode = "logical", length = length(output.1))
for (i in 1:length(output.1)) {
    output.2[i] <- (output.1[i] == output.1[i + 1])
}
output.2 %>% table()
# .
#  FALSE
# 107040
#HERE - 14: OK; tables the same

rm(pos, mpos, output.1, output.2, i)
