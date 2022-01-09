#!/usr/bin/env Rscript

library(parallel)
library(Rsamtools)
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


# Set up functions ------------------------------------------------------------
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
# chromosome <- "chr11"  #TODO No successful liftOvers: why? 
# chromosome <- "chr15"  #TODO No successful liftOvers: why?
chromosome <- "chrX"
file <- list.files(pattern = paste0("\\", chromosome, ".bam$")) %>%
    stringr::str_subset("^dedup\\.", negate = TRUE) %>% 
    stringr::str_subset("^mm10\\.", negate = TRUE)
variable <- file %>%
    stringr::str_split(., "-") %>%
    lapply(., `[[`, 1) %>%
    unlist() %>%
    paste0("dedup.", .)
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)


#  Assign .bam information as list --------------------------------------------

#  What to query from .bam files
map_info <- c(
    "qname", "flag", "rname", "strand", "pos", "mapq",
    "cigar", "mrnm", "mpos", "isize", "seq"
)
tag_info <- c("AS", "XS", "NM")
map_params <- Rsamtools::ScanBamParam(what = map_info, tag = tag_info)

command <- paste0(
    "<- ", variable, " %>% ",
        "Rsamtools::scanBam(., param = map_params)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

rm(map_info, tag_info, map_params)


#  Convert .bam information from list to dataframe to tibble ------------------
command <- paste0("<- ", variable, " %>% ", "as.data.frame() %>% as_tibble()")
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Removal of Picard MarkDuplicates flags -------------------------------------
dedup.129S1$flag %>% table()
duplicate_flags <- c(1107, 1123, 1171, 1187)
command <- paste0(
    "<- ", variable, " %>% dplyr::filter(flag %notin% duplicate_flags)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Reorder rname factor levels ------------------------------------------------
command <- paste0("<- forcats::fct_relevel(", variable, "$rname, chromosome)")
operation <- makeOperation(paste0(variable, "$rname"), command)
evaluateOperation(operation)


#  Drop unused rname factor levels --------------------------------------------
command <- paste0("<- ", variable, "$rname %>% forcats::fct_drop()")
operation <- makeOperation(paste0(variable, "$rname"), command)
evaluateOperation(operation)


#  Drop rows that are not chr1-19, chrX ---------------------------------------
chromosomes <- c(paste0("chr", 1:19), "chrX")
command <- paste0(
    "<- ", variable, " %>% ",
        "filter(., rname %in% chromosomes)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Create and append pos_end, mpos_end columns --------------------------------
command <- paste0("<- ", variable, "$pos + 49")
operation <- makeOperation(paste0(variable, "$pos_end"), command)
evaluateOperation(operation)

command <- paste0("<- ", variable, "$mpos + 49")
operation <- makeOperation(paste0(variable, "$mpos_end"), command)
evaluateOperation(operation)

command <- paste0(
    "<- ", variable, " %>% ",
        "dplyr::relocate(pos_end, .after = pos) %>% ",
        "dplyr::relocate(mpos_end, .after = mpos)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Order columns by coordinate, tag.AS, and mapq ------------------------------
command <- paste0(
    "<- ", variable, "[",
        "order(",
            variable, "$qname, ",
            variable, "$rname, ",
            variable, "$pos, ",
            variable, "$mpos, ",
            "-", variable, "$tag.AS, ",
            "-", variable, "$mapq",
        "), ",
    "]"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Based on coordinate value, identify if a given row is a duplicate ----------
# command <- paste0(
#     "<- ave(",
#         variable, "$coordinate, ",
#         variable, "$coordinate, ",
#         "FUN = length",
#     ") > 1L"
# )
# operation <- makeOperation(paste0(variable, "$duplicated"), command)
# evaluateOperation(operation)


#  Create a "temporary tag" for use in distinguishing entries
#+ 
#+ This is necessary to filter out duplicates and join liftOver data to initial
#+ tibbles; it's named "old_coordinate", i.e., the coordinate prior to liftOver
command <- paste0(
    "<- ", variable, " %>% ",
        "tidyr::unite(",
            "old_coordinate, ",
            "c(\"qname\", \"rname\", \"pos\", \"pos_end\"), ",
            "sep = \"_\", ",
            "remove = FALSE",
    ")"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Filter out all rows without unique coordinates -----------------------------
command <- paste0(
    "<- ", variable, " %>% ",
        "dplyr::distinct(old_coordinate, .keep_all = TRUE)"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)
dedup.129S1 %>% nrow()  # [1] 494694 (chr1); [1] 274025 (chrX, w/o rm MarkDuplicates); [1] 243239 (chrX)
dedup.CAST %>% nrow()  # [1] 492644 (chr1); [1] 281850 (chrX, w/o rm MarkDuplicates); [1] 249783 (chrX)
dedup.mm10 %>% nrow()  # [1] 491416 (chr1); [1] 283022 (chrX, w/o rm MarkDuplicates); [1] 252171 (chrX)

# -------
# STOP, GO TO SHORTCUT
# -------


#  Perform liftOvers of the 129S1 and CAST tibbles ----------------------------

#  Create and write out 129S1 and CAST .bed file(s) for liftOver
#+ pos
variable_bed <- paste0(
    variable[1:2], ".pos.", chromosome, ".bed"
)

command <- paste0(
    "<- tibble::tibble(",
        "\"rname\" = ", variable[1:2], "$rname, ",
        "\"pos\" = ", variable[1:2], "$pos, ",
        "\"pos_end\" = ", variable[1:2], "$pos_end, ",
        "\"old_coordinate\" = ", variable[1:2], "$old_coordinate",
    ")"
)
operation <- makeOperation(variable_bed, command)
evaluateOperation(operation)

operation <- paste0(
    "readr::write_tsv(",
        variable_bed, ", ",
        "file = \"", variable_bed, "\", ",
        "col_names = FALSE",
    ")"
)
evaluateOperation(operation)

#+ mpos
# variable_bed <- paste0(variable[1:2], ".mpos.", chromosome, ".bed")
# 
# command <- paste0(
#     "<- tibble::tibble(",
#         "\"rname\" = ", variable[1:2], "$rname, ",
#         "\"mpos\" = ", variable[1:2], "$mpos, ",
#         "\"mpos_end\" = ", variable[1:2], "$mpos_end, ",
#         "\"old_coordinate\" = ", variable[1:2], "$old_coordinate",
#     ")"
# )
# operation <- makeOperation(variable_bed, command)
# evaluateOperation(operation)
# 
# command <- paste0(
#     "readr::write_tsv(",
#         variable_bed, ", ",
#         "file = \"", variable_bed, "\", ",
#         "col_names = FALSE",
#     ")"
# )
# evaluateOperation(command)

# -------
#  Run shell code to perform the liftOver for 129S1
shell_code <- paste0(
"
# -------
#  Start recording time
start=$(date +%s)

# -------
#  Set up functions
renameChrCommon129S1() {
    #  Convert 129S1-SvImJ chromosomes names from common to official
    cat ${1} | \
    sed 's/chr1/CM003934.1/g' | \
    sed 's/chr2/CM003935.1/g' | \
    sed 's/chr3/CM003936.1/g' | \
    sed 's/chr4/CM003937.1/g' | \
    sed 's/chr5/CM003938.1/g' | \
    sed 's/chr6/CM003939.1/g' | \
    sed 's/chr7/CM003940.1/g' | \
    sed 's/chr8/CM003941.1/g' | \
    sed 's/chr9/CM003942.1/g' | \
    sed 's/chr10/CM003943.1/g' | \
    sed 's/chr11/CM003944.1/g' | \
    sed 's/chr12/CM003945.1/g' | \
    sed 's/chr13/CM003946.1/g' | \
    sed 's/chr14/CM003947.1/g' | \
    sed 's/chr15/CM003948.1/g' | \
    sed 's/chr16/CM003949.1/g' | \
    sed 's/chr17/CM003950.1/g' | \
    sed 's/chr18/CM003951.1/g' | \
    sed 's/chr19/CM003952.1/g' | \
    sed 's/chrX/CM003953.1/g' > ${2}
}


renameChr129S1Common() {
    #  Convert 129S1-SvImJ chromosome names from official to common
    cat ${1} | \
    sed 's/CM003934.1/chr1/g' | \
    sed 's/CM003935.1/chr2/g' | \
    sed 's/CM003936.1/chr3/g' | \
    sed 's/CM003937.1/chr4/g' | \
    sed 's/CM003938.1/chr5/g' | \
    sed 's/CM003939.1/chr6/g' | \
    sed 's/CM003940.1/chr7/g' | \
    sed 's/CM003941.1/chr8/g' | \
    sed 's/CM003942.1/chr9/g' | \
    sed 's/CM003943.1/chr10/g' | \
    sed 's/CM003944.1/chr11/g' | \
    sed 's/CM003945.1/chr12/g' | \
    sed 's/CM003946.1/chr13/g' | \
    sed 's/CM003947.1/chr14/g' | \
    sed 's/CM003948.1/chr15/g' | \
    sed 's/CM003949.1/chr16/g' | \
    sed 's/CM003950.1/chr17/g' | \
    sed 's/CM003951.1/chr18/g' | \
    sed 's/CM003952.1/chr19/g' | \
    sed 's/CM003953.1/chrX/g' > ${2}
}


displaySpinningIcon() {
    #  Display \"spinning icon\" while background process runs
    spin='-\\|/'
    i=0
    while kill -0 ${1} 2> /dev/null; do
        i=$(( (i + 1) % 4 ))
        printf \"\r${spin:$i:1} %s \" \"${2}\"
        sleep .1
    done
}


# -------
#  Set up paths, files
path_in=/Users/kalavattam/Downloads/to-do/get_unique_fragments/Bonora
path_out=/Users/kalavattam/Downloads/to-do/get_unique_fragments/Bonora

path_liftOver=/Users/kalavattam/Downloads/to-do/2021-1105-1107/liftOver

from_R=$(find . -maxdepth 1 -type f -name \"dedup.129S1.pos.", chromosome, ".bed\" -printf '%f\n')
prefix_from_R=${from_R%????}

# -------
#  Run liftOver loop
for file_prefix in \"${prefix_from_R[@]}\"; do
    #  Declare variables
    file_in=${file_prefix}.bed
    file_in_rename=${file_prefix}.rename.bed

    file_out_lifted=${file_prefix}.lifted.bed
    file_out_lifted_rename=${file_prefix}.lifted.rename.bed

    file_out_unlifted=${file_prefix}.unlifted.bed
    file_out_unlifted_rename=${file_prefix}.unlifted.rename.bed

    file_liftOver=129S1-SvImJ-to-mm10.over.chain.gz

    #  Variables for liftOver
    in=${path_in}/${file_in_rename}
    chain=${path_liftOver}/${file_liftOver}
    out_lift=${path_out}/${file_out_lifted}
    out_unlift=${path_out}/${file_out_unlifted}

    # -------
    #  Convert infile chromosomes names from common to official
    renameChrCommon129S1 ${file_in} ${file_in_rename}
    # head ${file_in_rename}

    # -------
    #  Do the liftOver in the background, saving those regions that are and are
    #+ not lifted
    liftOver -bedPlus=3 -tab ${in} ${chain} ${out_lift} ${out_unlift} &

    #  Display a spinning icon while liftOver is taking place
    displaySpinningIcon $! \"Lifting over $(basename ${in})\"

    # -------
    #  Convert the outfile chromosome names from official to common
    renameChr129S1Common ${file_out_lifted} ${file_out_lifted_rename} &
    displaySpinningIcon $! \"Renaming $(basename ${file_out_lifted})\"

    renameChr129S1Common ${file_out_unlifted} ${file_out_unlifted_rename} &
    displaySpinningIcon $! \"Renaming $(basename ${file_out_unlifted})\"

    # -------
    #  Clean up
    # rm ${file_in}
    rm ${file_in_rename}
    rm ${file_out_lifted}
    rm ${file_out_unlifted}

    mv ${file_out_lifted_rename} ${file_out_lifted}
    mv ${file_out_unlifted_rename} ${file_out_unlifted}
done

# -------
#  Record end time
end=$(date +%s)

#  Echo run time
run_time=$(echo ${end} - ${start} | bc -l)
echo ''
echo \"Run time: ${run_time} seconds\"
echo ''
"
)
system(shell_code)

# -------
#  Run shell code to perform the liftOver for CAST
shell_code <- paste0(
"
# -------
#  Start recording time
start=$(date +%s)

# -------
#  Set up functions
renameChrCommonCAST() {
    #  Convert CAST-EiJ chromosomes names from common to official
    cat ${1} | \
    sed 's/chr1/CM003994.1/g' | \
    sed 's/chr2/CM003995.1/g' | \
    sed 's/chr3/CM003996.1/g' | \
    sed 's/chr4/CM003997.1/g' | \
    sed 's/chr5/CM003998.1/g' | \
    sed 's/chr6/CM003999.1/g' | \
    sed 's/chr7/CM004000.1/g' | \
    sed 's/chr8/CM004001.1/g' | \
    sed 's/chr9/CM004002.1/g' | \
    sed 's/chr10/CM004003.1/g' | \
    sed 's/chr11/CM004004.1/g' | \
    sed 's/chr12/CM004005.1/g' | \
    sed 's/chr13/CM004006.1/g' | \
    sed 's/chr14/CM004007.1/g' | \
    sed 's/chr15/CM004008.1/g' | \
    sed 's/chr16/CM004009.1/g' | \
    sed 's/chr17/CM004010.1/g' | \
    sed 's/chr18/CM004011.1/g' | \
    sed 's/chr19/CM004012.1/g' | \
    sed 's/chrX/CM004013.1/g' > ${2}
}


renameChrCASTCommon() {
    #  Convert CAST-EiJ chromosome names from official to common
    cat ${1} | \
    sed 's/CM003994.1/chr1/g' | \
    sed 's/CM003995.1/chr2/g' | \
    sed 's/CM003996.1/chr3/g' | \
    sed 's/CM003997.1/chr4/g' | \
    sed 's/CM003998.1/chr5/g' | \
    sed 's/CM003999.1/chr6/g' | \
    sed 's/CM004000.1/chr7/g' | \
    sed 's/CM004001.1/chr8/g' | \
    sed 's/CM004002.1/chr9/g' | \
    sed 's/CM004003.1/chr10/g' | \
    sed 's/CM004004.1/chr11/g' | \
    sed 's/CM004005.1/chr12/g' | \
    sed 's/CM004006.1/chr13/g' | \
    sed 's/CM004007.1/chr14/g' | \
    sed 's/CM004008.1/chr15/g' | \
    sed 's/CM004009.1/chr16/g' | \
    sed 's/CM004010.1/chr17/g' | \
    sed 's/CM004011.1/chr18/g' | \
    sed 's/CM004012.1/chr19/g' | \
    sed 's/CM004013.1/chrX/g' > ${2}
}


displaySpinningIcon() {
    #  Display \"spinning icon\" while background process runs
    spin='-\\|/'
    i=0
    while kill -0 ${1} 2> /dev/null; do
        i=$(( (i + 1) % 4 ))
        printf \"\r${spin:$i:1} %s \" \"${2}\"
        sleep .1
    done
}


# -------
#  Set up paths, files
path_in=/Users/kalavattam/Downloads/to-do/get_unique_fragments/Bonora
path_out=/Users/kalavattam/Downloads/to-do/get_unique_fragments/Bonora

path_liftOver=/Users/kalavattam/Downloads/to-do/2021-1105-1107/liftOver

from_R=$(find . -maxdepth 1 -type f -name \"dedup.CAST.pos.", chromosome, ".bed\" -printf '%f\n')
prefix_from_R=${from_R%????}

# -------
#  Run liftOver loop
for file_prefix in \"${prefix_from_R[@]}\"; do
    #  Declare variables
    file_in=${file_prefix}.bed
    file_in_rename=${file_prefix}.rename.bed

    file_out_lifted=${file_prefix}.lifted.bed
    file_out_lifted_rename=${file_prefix}.lifted.rename.bed

    file_out_unlifted=${file_prefix}.unlifted.bed
    file_out_unlifted_rename=${file_prefix}.unlifted.rename.bed

    file_liftOver=CAST-EiJ-to-mm10.over.chain.gz

    #  Variables for liftOver
    in=${path_in}/${file_in_rename}
    chain=${path_liftOver}/${file_liftOver}
    out_lift=${path_out}/${file_out_lifted}
    out_unlift=${path_out}/${file_out_unlifted}

    # -------
    #  Convert infile chromosomes names from common to official
    renameChrCommonCAST ${file_in} ${file_in_rename}
    # head ${file_in_rename}

    # -------
    #  Do the liftOver in the background, saving those regions that are and are
    #+ not lifted
    liftOver -bedPlus=3 -tab ${in} ${chain} ${out_lift} ${out_unlift} &

    #  Display a spinning icon while liftOver is taking place
    displaySpinningIcon $! \"Lifting over $(basename ${in})\"

    # -------
    #  Convert the outfile chromosome names from official to common
    renameChrCASTCommon ${file_out_lifted} ${file_out_lifted_rename} &
    displaySpinningIcon $! \"Renaming $(basename ${file_out_lifted})\"

    renameChrCASTCommon ${file_out_unlifted} ${file_out_unlifted_rename} &
    displaySpinningIcon $! \"Renaming $(basename ${file_out_unlifted})\"

    # -------
    #  Clean up
    # rm ${file_in}
    rm ${file_in_rename}
    rm ${file_out_lifted}
    rm ${file_out_unlifted}

    mv ${file_out_lifted_rename} ${file_out_lifted}
    mv ${file_out_unlifted_rename} ${file_out_unlifted}
done

# -------
#  Record end time
end=$(date +%s)

#  Echo run time
run_time=$(echo ${end} - ${start} | bc -l)
echo ''
echo \"Run time: ${run_time} seconds\"
echo ''
"
)
system(shell_code)

# -------
# SHORTCUT
# -------
#  Read in the "lifted" data
file <- list.files(pattern = paste0("\\.", chromosome,".lifted.bed$"))
variable_lifted <- file
mapply(
    assign, variable_lifted, file, MoreArgs = list(envir = parent.frame())
)

command <- paste0(
    "<- readr::read_tsv(",
        "file = \"", file, "\", ",
        "col_names = c(",
            "\"liftOver_rname\", ",
            "\"liftOver_pos\", ",
            "\"liftOver_pos_end\", ",
            "\"old_coordinate\"",
        ")",
    ")"
)
operation <- makeOperation(variable_lifted, command)
evaluateOperation(operation)

command <- paste0("<- \"#liftOver successful\"")
operation <- makeOperation(paste0(file, "$liftOver_reason"), command)
evaluateOperation(operation)

# -------
# STOP, GO TO SHORTCUT
# -------


# -------
#  Munge the "unlifted" data
shell_code <- paste0(
"
set +o posix

# -------
#  Start recording time
start=$(date +%s)

# -------
#  Set up files
typeset -a from_R
while IFS=  read -r -d $'\\0'; do
    from_R+=(\"$REPLY\")
done < <(find . -maxdepth 1 -type f -name \"dedup.*.pos.", chromosome, ".bed\" -print0)

typeset -a prefixes_from_R
for i in \"${from_R[@]}\"; do prefixes_from_R+=( \"${i%????}\" ); done

# -------
#  Loop for munging files
for file_prefix in \"${prefixes_from_R[@]}\"; do
    file_out_unlifted=${file_prefix}.unlifted.bed

    # -------
    #  Append tab to the end of each line
    cat ${file_out_unlifted} | sed 's/$/\t/' > ${file_out_unlifted/.bed/.sed.bed}

    #  Move every odd row down to the beginning of the subsequent even row
    awk 'NR%2==0{print p,$0}{p=$0}' ${file_out_unlifted/.bed/.sed.bed} > ${file_out_unlifted/.bed/.munged.bed}

    # -------
    #  Clean up
    rm ${file_out_unlifted}
    rm ${file_out_unlifted/.bed/.sed.bed}
    mv ${file_out_unlifted/.bed/.munged.bed} ${file_out_unlifted}
done

# -------
#  Record end time
end=$(date +%s)

#  Echo run time
run_time=$(echo ${end} - ${start} | bc -l)
echo ''
echo \"Run time: ${run_time} seconds\"
echo ''
"
)
system(shell_code)

# -------
# SHORTCUT
# -------
#  Read in the "unlifted" data
file <- list.files(pattern = paste0("\\.", chromosome,".unlifted.bed$"))
variable_unlifted <- file[1:length(file)]
mapply(
    assign, variable_unlifted, file, MoreArgs = list(envir = parent.frame())
)

command <- paste0(
    "<- readr::read_tsv(",
        "file = \"", file, "\", ",
        "col_names = c(",
            "\"liftOver_reason\", ",
            "\"liftOver_rname\", ",
            "\"liftOver_pos\", ",
            "\"liftOver_pos_end\", ",
            "\"old_coordinate\"",
        ")",
    ")"
)
operation <- makeOperation(variable_unlifted, command)
evaluateOperation(operation)


#  Munge the information in the "unlifted" 129S1 and CAST variables
command <- paste0("<- NA_character_")
operation <- makeOperation(
    paste0(variable_unlifted, "$liftOver_rname"), command
)
evaluateOperation(operation)

command <- paste0("<- NA_integer_")
operation <- makeOperation(
    paste0(variable_unlifted, "$liftOver_pos"), command
)
evaluateOperation(operation)

command <- paste0("<- NA_integer_")
operation <- makeOperation(
    paste0(variable_unlifted, "$liftOver_pos_end"), command
)
evaluateOperation(operation)

command <- paste0(
    "<- dplyr::bind_rows(", variable_unlifted, ", ", variable_lifted, ") %>% ",
        "full_join(., ", variable[1:2], ", by = \"old_coordinate\") %>% ",
        "dplyr::arrange(\"old_coordinate\")"
)
operation <- makeOperation(variable[1:2], command)
evaluateOperation(operation)

command <- paste0("<- ", variable[1:2], "$liftOver_reason %>% as_factor()")
operation <- makeOperation(paste0(variable[1:2], "$liftOver_reason"), command)
evaluateOperation(operation)

command <- paste0(
    "<- ", variable[1:2], " %>% ",
        "dplyr::relocate(",
            "c(liftOver_rname, liftOver_pos, liftOver_pos_end, liftOver_reason), ",
            ".after = pos_end",
        ")"
)
operation <- makeOperation(variable[1:2], command)
evaluateOperation(operation)

dedup.129S1$liftOver_reason <- plyr::revalue(
    dedup.129S1$liftOver_reason,
    c(
        "#Deleted in new" = "Deleted in new",
        "#Partially deleted in new" = "Partially deleted in new",
        "#Split in new" = "Split in new",
        "#liftOver successful" = "liftOver successful"
    )
)

dedup.CAST$liftOver_reason <- plyr::revalue(
    dedup.CAST$liftOver_reason,
    c(
        "#Deleted in new" = "Deleted in new",
        "#Partially deleted in new" = "Partially deleted in new",
        "#Split in new" = "Split in new",
        "#liftOver successful" = "liftOver successful"
    )
)


#  Create back-ups ------------------------------------------------------------
# bak.dedup.129S1 <- dedup.129S1
# bak.dedup.CAST <- dedup.CAST
# bak.dedup.mm10 <- dedup.mm10

# dedup.129S1 <- bak.dedup.129S1
# dedup.CAST <- bak.dedup.CAST
# dedup.mm10 <- bak.dedup.mm10


#  Check on NAs and duplicates... ---------------------------------------------
test <- dedup.129S1[is.na(dedup.129S1$liftOver_reason), ]  # 9275 first time, 9349 the next, 9386, 9404; 0 once the bug is fixed...
test.dup <- test %>%
    dplyr::filter(old_coordinate %in% unique(.[["old_coordinate"]][
        duplicated(.[["old_coordinate"]])
    ]))  # 0

test <- dedup.CAST[is.na(dedup.CAST$liftOver_reason), ]  # 17809 first time, 17950 the next, 18021, 18056; 0 once the bug is fixed...
test.dup <- test %>%
    dplyr::filter(old_coordinate %in% unique(.[["old_coordinate"]][
        duplicated(.[["old_coordinate"]])
    ]))  # 0

rm(test, test.dup)
#NOTE Bug arises from running the '#  Munge the "unlifted" data' code chunk more than once
#TODO If not 0, then stop the program and print a warning


#  Clean up both CAST and 129S1 -----------------------------------------------
rm(
    command,
    file,
    operation,
    shell_code
)

operation <- paste0(
    "rm(",
        variable_bed, ", ",
        variable_lifted, ", ",
        variable_unlifted,
    ")"
)
evaluateOperation(operation)

rm(
    variable,
    variable_bed,
    variable_lifted,
    variable_unlifted
)

# system("rm *.bed")


#  Munging and creation of "read" and coordinate "columns" --------------------
#+ 
#+ 1. Reorder the pos_end and mpos_end columns
#+ 2. Create a column of concatenated barcode_sequence strings (i.e., "reads")  #NOTE No longer doing this
#+ 3. And then create a column of concatenated barcode_fragment strings (i.e.,
#+    "coordinates," where a single coordinate is rname, pos, and pos_end)
#+ 
#+ #NOTE
#+ - Terminology change: "b_s" is now "read"
#+ - Terminology change: "b_f" is now "coordinate"
variable <- c("dedup.129S1", "dedup.CAST")
command <- paste0(
    "<- ", variable, " %>% ",
        "dplyr::rename(",
            "old_rname = rname, ",
            "old_pos = pos, ",
            "old_pos_end = pos_end, ",
            "old_mpos = mpos, ",
            "old_mpos_end = mpos_end, ",
            "rname = liftOver_rname, ",
            "pos = liftOver_pos, pos_end = liftOver_pos_end",
        ")"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)

variable <- c("dedup.129S1", "dedup.CAST", "dedup.mm10")
command <- paste0(
    "<- ", variable, " %>% ",
        "tidyr::unite(",
            "coordinate, ",
            "c(\"qname\", \"rname\", \"pos\", \"pos_end\"), ",
            "sep = \"_\", ",
            "remove = FALSE",
        ")"
)
operation <- makeOperation(variable, command)
evaluateOperation(operation)


#  Save the tibbles containing all liftOver assignments as .rds files ---------
#+ 
#+ Save the tibble for alignments to the N-masked mm10 genome too, even though
#+ no liftOvers were performed on it
operation <- paste0(
    "saveRDS(",
        variable, ", ",
        "file = '", variable, ".", chromosome, ".post-liftOvers.rds'",
    ")"
)
evaluateOperation(operation)
#TODO 1/2 Example of what these output files look like and/or some systematic
#TODO 2/2 way to tell they're from this script/step in the pipeline

#  Output is, e.g., "dedup.129S1.chrX.post-liftOvers.rds"

rm(list = ls())
