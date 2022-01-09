# -----------------------------------------------------------------------------
#  Scraps, TBD

#  List .bams in directory
# list.files(pattern = "\\.bam$")

# x <- as.list(rnorm(10))
# names(x) <- paste("bam", 1:length(x), sep = "")
# list2env(x , envir = .GlobalEnv)

# x <- list.files(pattern = "\\.bam$") %>% as.list()
# names(x) <- paste(x[1:length(x)])
# list2env(x , envir = .GlobalEnv)
# x %>% rm()

#  Read in an entire .bam
# `mm10-CAST-129-Nmasked.F121-6-CASTx129.undifferentiated.dedup.bam` <-
#     `mm10-CAST-129-Nmasked.F121-6-CASTx129.undifferentiated.dedup.bam` %>% scanBam()
# 
# `mm10.F121-6-CASTx129.undifferentiated.dedup.bam` <-
#     `mm10.F121-6-CASTx129.undifferentiated.dedup.bam` %>% scanBam()

# test <- scanBam(list.files(pattern = "\\.bam$")[1], param = map_params)

# -------------------------------------
#  Scraps: Learning to work with mapply()
# Vars <- c("_car","_bike","_lorry")
# Dat <- c(10,20,22)
# 
# for (i in 1:length(Vars)) { assign(Vars[i], Dat[i]) }
# rm(i, Vars, Dat, "_car", "_bike", "_lorry")
# 
# mapply(assign, Vars, Dat, MoreArgs=list(envir=parent.frame()))
# mapply(assign, Vars, Dat)
# rm(i, Vars, Dat, "_car", "_bike", "_lorry")

# -------------------------------------


# -------------------------------------
# bam1$qname %>% unique() %>% length()  # [1] 2972
# bam2$qname %>% unique() %>% length()  # [1] 2961
# bam3$qname %>% unique() %>% length()  # [1] 2945
# bam4$qname %>% unique() %>% length()  # [1] 2959

# bam1$seq %>% unique() %>% length()  # [1] 6326297
# bam2$seq %>% unique() %>% length()  # [1] 6314847
# bam3$seq %>% unique() %>% length()  # [1] 6337920
# bam4$seq %>% unique() %>% length()  # [1] 6389895


# -------------------------------------
# bam1$qname[!(bam1$qname %in% bam2$qname)] %>% length()  # [1] 28
# bam2$qname[!(bam2$qname %in% bam1$qname)] %>% length()  # [1] 58

# test <- bam1 %>%
#     group_by(qname) %>%
#     summarize(
#         mean = mean(tag.AS),
#         median = median(tag.AS),
#         min = min(tag.AS),
#         max = max(tag.AS),
#         sd = sd(tag.AS),
#         iqr = IQR(tag.AS),
#         # n = n(tag.AS),
#         n_distinct = n_distinct(tag.AS),
#         
#     )
# 
# test %>% head()


# -----------------------------------------------------------------------------
#  Working with *.fragments.txt.gz
files_gz <- list.files(pattern = "\\.gz$")
n_gz <- files_gz %>% length()
variables_gz <- paste0("gz", 1:n_gz)
command_pipe <- paste0("<- gz", 1:n_gz, " %>% ")

#  Read in .txt.gz filenames
command <- paste0("<- files_gz[", 1:n_gz, "]")
operation <- make_operation(variables_gz, command)
eval(parse(text = operation))

command <- paste0(command_pipe, "gzfile(., \"rt\")")
operation <- make_operation(variables_gz, command)
eval(parse(text = operation))

command <- paste0(command_pipe, "read.delim(., header = FALSE)")
operation <- make_operation(variables_gz, command)
eval(parse(text = operation))

#  Sort rows by rname, pos, mpos, and qual
command <- paste0(command_pipe, "arrange(., V1, V2, V3, V4, V5)")
operation <- make_operation(variables_gz, command)
eval(parse(text = operation))

gz1 %>% head()
gz2 %>% head()
gz3 %>% head()
gz4 %>% head()

gz1$V4 %>% unique() %>% length()  # [1] 382
gz2$V4 %>% unique() %>% length()  # [1] 382
gz3$V4 %>% unique() %>% length()  # [1] 382
gz4$V4 %>% unique() %>% length()  # [1] 382

# -------------------------------------
test <- bam1$qname %>% head()
test_qname <- gsub(':.*', '', test)

bam1_qname <- gsub(':.*', '', bam1$qname)
bam1_qname[!(bam1_qname %in% gz1$V4)]  # character(0)
bam1_qname %>% unique() %>% length()  # [1] 382

bam2_qname <- gsub(':.*', '', bam2$qname)
bam2_qname[!(bam2_qname %in% gz2$V4)]  # character(0)
bam2_qname %>% unique() %>% length()  # [1] 382

bam3_qname <- gsub(':.*', '', bam3$qname)
bam3_qname[!(bam3_qname %in% gz3$V4)]  # character(0)
bam3_qname %>% unique() %>% length()  # [1] 382

bam4_qname <- gsub(':.*', '', bam4$qname)
bam4_qname[!(bam4_qname %in% gz4$V4)]  # character(0)
bam4_qname %>% unique() %>% length()  # [1] 382


# -----------------------------------------------------------------------------
#  View names of .bam fields
# bam[[1]] %>% names()
# [1] "qname"  "flag"   "rname"  "strand" "pos"    "qwidth" "mapq"   "cigar"  "mrnm"   "mpos"  "isize"  "seq"    "qual"

#  Check distribution of .bam flags
bam[[1]]$flag %>% table()
# .
# 83      99     147     163 
# 2478882 2474233 2474238 2478868
#  URL to understand bit-wise flags: samformat.info/sam-format-flag

.unlist <- function (x) {
  #  Function to collapse list of lists into a single list (as per the
  #+ Rsamtools vignette)
  #+
  #+ do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#  Store names of .bam fields
bam_field <- bam[[1]] %>% names()

#  Go through each .bam field and unlist it
list <- lapply(
  bam_field,
  function(y) .unlist(lapply(bam, "[[", y))
)

#  Store the list as a data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

bam_df %>% dim()
# [1] 9906221      13

bam_df$qname %>% length()
# [1] 9906221

bam_df$qname %>% head()
# [1] "CGGCTATGAGTTCTCATGCAGCCGGCTTGTACTGAC:0" "CGGCTATGAGTTCTCATGCAGCCGGCTTGTACTGAC:0" "CGGCTATGAGTTCTCATGCAGCCGGCTTGTACTGAC:0"
# [4] "CGGCTATGAGTTCTCATGCAGCCGGCTTGTACTGAC:0" "CGGCTATGAGTTCTCATGCAGCCGGCTTGTACTGAC:0" "CGGCTATGAGTTCTCATGCAGCCGGCTTGTACTGAC:0"

bam_df$flag %>% head()
# [1]  99  99  99  99 147 147

bam_df$flag %>% table()
# .
# 83      99     147     163 
# 2478882 2474233 2474238 2478868

bam_df$rname %>% head()
# [1] chr10 chr10 chr10 chr10 chr10 chr10
# Levels: chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chrM chrX chrY

bam_df$strand %>% head()
# [1] + + + + - -
#    Levels: + - *

bam_df$pos %>% head()
# > bam_df$pos %>% head()
# [1] 3117273 3117273 3117273 3117273 3117434 3117434

bam_df$qwidth %>% head()
# [1] 50 50 50 50 50 50

bam_df$mapq %>% head()
# [1] 42 42 42 42 42 42

bam_df$mapq %>% table()
# .
# 11      12      14      15      16      17      18      21      22      23      24      25      26      27      30      31      32      33      34 
# 21118   48453    3598   12170    4092   45422   46562   36328   18092  145919  147274   30626   22560   11444   53193   47392   41958    6348   39349 
# 35      36      37      38      39      40      42 
# 43543   11728   16750   11586   10080  404924 8625712

bam_df$cigar %>% head()
# [1] "50M"      "50M"      "50M"      "50M"      "12M2I36M" "12M2I36M"

# bam_df$cigar %>% table()

bam_df$mrnm %>% head()
# [1] chr10 chr10 chr10 chr10 chr10 chr10
# 22 Levels: chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr1 chr2 chr3 ... chrY

bam_df$mpos %>% head()
# [1] 3117434 3117434 3117434 3117434 3117273 3117273

bam_df$isize %>% head()
# [1]  209  209  209  209 -209 -209

bam_df$seq %>% head()
# DNAStringSet object of length 6:
#     width seq
# [1]    50 TGCCAAATGTCTAAACCACAACTCAGGAATTTTCTCAGAACTGAGCAGTA
# [2]    50 TGCCAAATGTCTAAACCACAACTCAGGAATTTTCTCAGAACTGAGCAGTA
# [3]    50 TGCCAAATGTCTAAACCACAACTCAGGAATTTTCTCAGAACTGAGCAGTA
# [4]    50 TGCCAAATGTCTAAACCACAACTCAGGAATTTTCTCAGAACTGAGCAGTA
# [5]    50 CACTCTGTAGACCAGGCTGGCCTTGAACTCAGGATCTGCCTGCTGCTGCC
# [6]    50 CACTCTGTAGACCAGGCTGGCCTTGAACTCAGGATCTGCCTGCTGCTGCC

bam_df$qual %>% head()
# PhredQuality object of length 6:
#     width seq
# [1]    50 AAAAAEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# [2]    50 A/AAAEEEEEEEAEEEEEEEEEEEEEEEE/EAAEEEEAEEEEE6EEAE/A
# [3]    50 AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# [4]    50 AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# [5]    50 EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA
# [6]    50 A/EEA/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEAAAAA


# -----------------------------------------------------------------------------
#  - stackoverflow.com/questions/44652234/read-bam-with-tags-to-tbl-recursive-dplyrbind-cols-for-list-of-lists
#+ - support.bioconductor.org/p/49521/

# map_info <- c("rname", "strand", "pos")
map_info <- c("qname", "flag", "rname", "pos", "cigar", "mpos", "isize", "seq", "qual")
# tag_info <- c("AS", "XS", "MD", "NM")
tag_info <- "AS"

map_params <- ScanBamParam(what = map_info, tag = tag_info)

test <- scanBam(list.files(pattern = "\\.bam$")[1], param = map_params)
# test_tbl <- bind_cols(
#     do.call(bind_cols, bam[[1]][map_info]),
#     do.call(bind_cols, bam[[1]]$tag)
# )

test_df <- test %>% as.data.frame()
test_df %>% dim()
# [1] 9906221       3
