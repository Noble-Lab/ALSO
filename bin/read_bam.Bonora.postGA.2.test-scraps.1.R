#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
#  ...from read_bam.Bonora.3.R

#  Load libraries
library(dplyr)
library(forcats)
library(magrittr)
library(parallel)
library(purrr)
library(Rsamtools)

options(pillar.sigfig = 8, scipen = 10000)


# -----------------------------------------------------------------------------
#  Temporary: Set up work directory (location TB∆)
"/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora" %>% setwd()


# -----------------------------------------------------------------------------
#  Set up functions
`%notin%` <- Negate(`%in%`)

makeOperation <- function(variable, command) {
    operation <- paste(variable, command) %>% paste(., collapse = "; ")
    return(operation)
}


#  Set up what to query from .bam files ---------------------------------------
map_info <- c("qname", "flag", "rname", "strand", "pos", "mapq", "cigar", "mrnm", "mpos", "isize", "seq")
tag_info <- c("AS", "XS", "NM")
map_params <- Rsamtools::ScanBamParam(what = map_info, tag = tag_info)


#  Assign variables necessary for subsequent commands -------------------------
chromosome <- "chr1"
# chromosome <- "chrX"
# chromosome <- "chr11"  #TODO No successful liftOvers: why? 
# chromosome <- "chr15"  #TODO No successful liftOvers: why?
file <- list.files(pattern = paste0("\\", chromosome, ".bam$"))
file <- file[c(1, 2, 3)]
variable <- file %>% strsplit(., "-") %>% lapply(., `[[`, 1) %>% unlist() %>% paste0("dedup.", .)

mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)


#  Assign .bam information as list --------------------------------------------
command_pipe <- paste0("<- ", variable, " %>% ")
command <- paste0(command_pipe, "Rsamtools::scanBam(., param = map_params)")
operation <- makeOperation(variable, command)
eval(parse(text = operation))

rm(map_info, tag_info, map_params)


#  Convert .bam information from list to dataframe to tibble ------------------
command <- paste0(command_pipe, "as.data.frame() %>% as_tibble()")
operation <- makeOperation(variable, command)
eval(parse(text = operation))


#  Reorder rname factor levels ------------------------------------------------
command <- paste0("<- forcats::fct_relevel(", variable, "$rname, chromosome)")
operation <- makeOperation(paste0(variable, "$rname"), command)
eval(parse(text = operation))


#  Drop unused rname factor levels --------------------------------------------
command <- paste0("<- ", variable, "$rname %>% forcats::fct_drop()")
operation <- makeOperation(paste0(variable, "$rname"), command)
eval(parse(text = operation))


#  Drop rows that are not chr1-19, chrX ---------------------------------------
chromosomes <- c(paste0("chr", 1:19), "chrX")
command <- paste0(command_pipe, "filter(., rname %in% chromosomes)")
operation <- makeOperation(variable, command)
eval(parse(text = operation))


#  Create and append pos_end, mpos_end columns --------------------------------
command <- paste0("<- ", variable, "$pos + 49")
operation <- makeOperation(paste0(variable, "$pos_end"), command)
eval(parse(text = operation))

command <- paste0("<- ", variable, "$mpos + 49")
operation <- makeOperation(paste0(variable, "$mpos_end"), command)
eval(parse(text = operation))

command <- paste0(
    command_pipe, "dplyr::relocate(pos_end, .after = pos) %>% dplyr::relocate(mpos_end, .after = mpos)"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))


#  Order columns by coordinate, tag.AS, and mapq ------------------------------
command <- paste0(
    "<- ", variable, "[order(", variable, "$qname, ", variable, "$rname, ", variable, "$pos, ", variable, "$mpos, -", variable, "$tag.AS, -", variable, "$mapq), ]"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))


#  Based on coordinate value, identify if a given row is a duplicate ----------
# command <- paste0("<- ave(", variable, "$coordinate, ", variable, "$coordinate, FUN = length) > 1L")
# operation <- makeOperation(paste0(variable, "$duplicated"), command)
# eval(parse(text = operation))


#  Create separate j's for duplicated and non-duplicated observations ---------
# command <- paste0(
#     command_pipe,
#     "dplyr::group_by(duplicated) %>% ",
#     "dplyr::group_split()"
# )
# operation <- makeOperation(paste0(variable, ".split"), command)
# eval(parse(text = operation))
# 
# command <- paste0("<- ", variable, ".split[[1]]")
# operation <- makeOperation(paste0(variable, ".nondup"), command)
# eval(parse(text = operation))
# 
# command <- paste0("<- ", variable, ".split[[2]]")
# operation <- makeOperation(paste0(variable, ".dup"), command)
# eval(parse(text = operation))
# 
# command <- paste0("rm(", variable, ".split)")
# eval(parse(text = command))


#  Create categories based on the number of rows comprising a "duplicate group"
# variable <- paste0(variable, ".dup")
# command_pipe <- paste0("<- ", variable, " %>% ")
# command <- paste0(
#     command_pipe,
#     "group_by(coordinate) %>% ",
#     "dplyr::mutate(group_n_members = n()) %>% ",
#     "dplyr::mutate_at(vars(group_n_members), dplyr::funs(factor)) %>% ",
#     "dplyr::mutate(group_n_members = recode(group_n_members, \"1\" = \"one\", \"2\" = \"two\", \"3\" = \"three\"))"
# )
# operation <- makeOperation(variable, command)
# eval(parse(text = operation))

#NOTE The above hashed-out code is not an appropriate strategy since duplicate groups can have up to 13 rows...


#  Create a "temporary tag" for use in distinguishing entries
#+ 
#+ This is necessary to filter out duplicates and join liftOver data to initial
#+ tibbles
command <- paste0(
    command_pipe, "tidyr::unite(old_coordinate, c(\"qname\", \"rname\", \"pos\", \"pos_end\"), sep = \"_\", remove = FALSE)"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))


#  Filter out all rows without unique coordinates -----------------------------
command <- paste0(command_pipe, "dplyr::distinct(old_coordinate, .keep_all = TRUE)")
operation <- makeOperation(variable, command)
eval(parse(text = operation))
dedup.129S1 %>% nrow()  # [1] 494694 (chr1); [1] 274025 (chrX); [1] 260377
dedup.CAST %>% nrow()  # [1] 492644 (chr1); [1] 281850 (chrX); [1] 261514
dedup.mm10 %>% nrow()  # [1] 491416 (chr1); [1] 283022 (chrX); [1] 262658


#  Perform liftOvers of the 129S1 and CAST tibbles ----------------------------

#  Create and write out 129S1 and CAST .bed file(s) for liftOver
#+ pos
variable_bed <- paste0(
    variable[1:2], ".pos.", chromosome, ".bed"
)

command <- paste0(
    "<- tibble::tibble(\"rname\" = ", variable[1:2], "$rname, ",
    "\"pos\" = ", variable[1:2], "$pos, ",
    "\"pos_end\" = ", variable[1:2], "$pos_end, ",
    "\"old_coordinate\" = ", variable[1:2], "$old_coordinate)"
)
operation <- makeOperation(variable_bed, command)
eval(parse(text = operation))

command <- paste0(
    "readr::write_tsv(", variable_bed, ", file = \"", variable_bed, "\", col_names = FALSE)"
)
eval(parse(text = command))

#+ mpos
# variable_bed <- paste0(
#     variable[1:2], ".mpos.", chromosome, ".bed"
# )
#
# command <- paste0(
#     "<- tibble::tibble(\"rname\" = ", variable[1:2], "$rname, ",
#     "\"mpos\" = ", variable[1:2], "$mpos, ",
#     "\"mpos_end\" = ", variable[1:2], "$mpos_end, ",
#     "\"old_coordinate\" = ", variable[1:2], "$old_coordinate)"
# )
# operation <- makeOperation(variable_bed, command)
# eval(parse(text = operation))
#
# command <- paste0(
#     "readr::write_tsv(", variable_bed, ", file = \"", variable_bed, "\", col_names = FALSE)"
# )
# eval(parse(text = command))

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
    "<- readr::read_tsv(file = \"", file, "\", col_names = c(\"liftOver_rname\", \"liftOver_pos\", \"liftOver_pos_end\", \"old_coordinate\"))"
)
operation <- makeOperation(variable_lifted, command)
eval(parse(text = operation))

command <- paste0("<- \"#liftOver successful\"")
operation <- makeOperation(paste0(file, "$liftOver_reason"), command)
eval(parse(text = operation))

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
# Read in the "unlifted" data
file <- list.files(pattern = paste0("\\.", chromosome,".unlifted.bed$"))
variable_unlifted <- file[1:length(file)]
mapply(
    assign, variable_unlifted, file, MoreArgs = list(envir = parent.frame())
)

command <- paste0(
    "<- readr::read_tsv(file = \"", file, "\", col_names = c(\"liftOver_reason\", \"liftOver_rname\", \"liftOver_pos\", \"liftOver_pos_end\", \"old_coordinate\"))"
)
operation <- makeOperation(variable_unlifted, command)
eval(parse(text = operation))


#  Munge the information in the "unlifted" 129S1 and CAST variables
command <- paste0("<- NA_character_")
operation <- makeOperation(paste0(variable_unlifted, "$liftOver_rname"), command)
eval(parse(text = operation))

command <- paste0("<- NA_integer_")
operation <- makeOperation(paste0(variable_unlifted, "$liftOver_pos"), command)
eval(parse(text = operation))

command <- paste0("<- NA_integer_")
operation <- makeOperation(paste0(variable_unlifted, "$liftOver_pos_end"), command)
eval(parse(text = operation))

command <- paste0(
    "<- dplyr::bind_rows(", variable_unlifted, ", ", variable_lifted, ") %>% ",
    "full_join(., ", variable[1:2], ", by = \"old_coordinate\") %>% ",
    "dplyr::arrange(\"old_coordinate\")"
)
operation <- makeOperation(variable[1:2], command)
eval(parse(text = operation))

command <- paste0("<- ", variable[1:2], "$liftOver_reason %>% as_factor()")
operation <- makeOperation(paste0(variable[1:2], "$liftOver_reason"), command)
eval(parse(text = operation))

command <- paste0(
    command_pipe[1:2],
    "dplyr::relocate(c(liftOver_rname, liftOver_pos, liftOver_pos_end, liftOver_reason), .after = pos_end)"
)
operation <- makeOperation(variable[1:2], command)
eval(parse(text = operation))

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
bak.dedup.129S1 <- dedup.129S1
bak.dedup.CAST <- dedup.CAST
bak.dedup.mm10 <- dedup.mm10

dedup.129S1 <- bak.dedup.129S1
dedup.CAST <- bak.dedup.CAST
dedup.mm10 <- bak.dedup.mm10

#  Checking on NAs and duplicates... ------------------------------------------
#TODO
test <- dedup.129S1[is.na(dedup.129S1$liftOver_reason), ]  # 9275 first time, 9349 the next, 9386, 9404; 0 once the bug is fixed...
test.dup <- test %>% dplyr::filter(old_coordinate %in% unique(.[["old_coordinate"]][duplicated(.[["old_coordinate"]])]))  # 0

test <- dedup.CAST[is.na(dedup.CAST$liftOver_reason), ]  # 17809 first time, 17950 the next, 18021, 18056; 0 once the bug is fixed...
test.dup <- test %>% dplyr::filter(old_coordinate %in% unique(.[["old_coordinate"]][duplicated(.[["old_coordinate"]])]))  # 0


# -------
#  Clean up both CAST and 129S1
rm(
    command,
    command_pipe,
    file,
    dedup.129S1.pos.bed,
    dedup.129S1.pos.lifted.bed,
    dedup.129S1.pos.unlifted.bed,
    dedup.CAST.pos.bed,
    dedup.CAST.pos.lifted.bed,
    dedup.CAST.pos.unlifted.bed,
    operation,
    shell_code,
)

rm(
    test,
    test.dup
)

# system("rm *.bed")


#  Munging and creation of "read" and coordinate "columns" --------------------
#+ 
#+ Reorder the pos_end and mpos_end columns; create a column of concatenated
#+ barcode_sequence strings (i.e., "reads"); and then create a column of
#+ concatenated barcode_fragment strings (i.e., "coordinates," where a single
#+ coordinate is rname, pos, and pos_end)
#+ 
#+ #NOTE
#+ - Terminology change: "b_s" is now "read"
#+ - Terminology change: "b_f" is now "coordinate"
variable <- c("dedup.129S1", "dedup.CAST")
command_pipe <- paste0("<- ", variable, " %>% ")
command <- paste0(
    command_pipe, "dplyr::rename(old_rname = rname, old_pos = pos, old_pos_end = pos_end, old_mpos = mpos, old_mpos_end = mpos_end, rname = liftOver_rname, pos = liftOver_pos, pos_end = liftOver_pos_end)"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))

variable <- c("dedup.129S1", "dedup.CAST", "dedup.mm10")
command_pipe <- paste0("<- ", variable, " %>% ")
command <- paste0(
    command_pipe,
    "tidyr::unite(coordinate, c(\"qname\", \"rname\", \"pos\", \"pos_end\"), sep = \"_\", remove = FALSE)"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))


# # Checking the liftOvers ------------------------------------------------------
# test.1 <- dedup.129S1 %>% group_by(liftOver_reason) %>% group_split()
# 
# test.1.deleted_in_new <- test.1[[1]]
# test.1.deleted_in_new.dup <- test.1.deleted_in_new %>% dplyr::filter(coordinate %in% unique(.[["coordinate"]][duplicated(.[["coordinate"]])]))
# 
# test.1.split_in_new <- test.1[[2]]
# test.1.split_in_new.dup <- test.1.split_in_new %>% dplyr::filter(coordinate %in% unique(.[["coordinate"]][duplicated(.[["coordinate"]])]))
# 
# test.1.partially_deleted <- test.1[[3]]
# test.1.partially_deleted.dup <- test.1.partially_deleted %>% dplyr::filter(coordinate %in% unique(.[["coordinate"]][duplicated(.[["coordinate"]])]))
# 
# test.1.liftOver_successful <- test.1[[4]]
# test.1.liftOver_successful.dup <- test.1.liftOver_successful %>% dplyr::filter(coordinate %in% unique(.[["coordinate"]][duplicated(.[["coordinate"]])]))
# 
# # --
# test.C <- dedup.CAST %>% group_by(liftOver_reason) %>% group_split()
# 
# test.C.deleted_in_new <- test.C[[1]]
# test.C.deleted_in_new.dup <- test.C.deleted_in_new %>% dplyr::filter(coordinate %in% unique(.[["coordinate"]][duplicated(.[["coordinate"]])]))
# 
# test.C.split_in_new <- test.C[[2]]
# test.C.split_in_new.dup <- test.C.split_in_new %>% dplyr::filter(coordinate %in% unique(.[["coordinate"]][duplicated(.[["coordinate"]])]))
# 
# test.C.partially_deleted <- test.C[[3]]
# test.C.partially_deleted.dup <- test.C.partially_deleted %>% dplyr::filter(coordinate %in% unique(.[["coordinate"]][duplicated(.[["coordinate"]])]))
# 
# test.C.liftOver_successful <- test.C[[4]]
# test.C.liftOver_successful.dup <- test.C.liftOver_successful %>% dplyr::filter(coordinate %in% unique(.[["coordinate"]][duplicated(.[["coordinate"]])]))
# 
# rm(
#     test.1,
#     test.1.deleted_in_new,
#     test.1.deleted_in_new.dup,
#     test.1.liftOver_successful,
#     test.1.liftOver_successful.dup,
#     test.1.split_in_new,
#     test.1.split_in_new.dup,
#     test.1.partially_deleted,
#     test.1.partially_deleted.dup,
#     test.C,
#     test.C.deleted_in_new,
#     test.C.deleted_in_new.dup,
#     test.C.liftOver_successful,
#     test.C.liftOver_successful.dup,
#     test.C.partially_deleted,
#     test.C.partially_deleted.dup,
#     test.C.split_in_new,
#     test.C.split_in_new.dup
# )


#  Work with only the "liftOver successful" observations ----------------------
l.variable <- paste0("l.", variable)[c(1, 2)]
command <- paste0("<- ", variable[c(1, 2)], " %>% group_by(liftOver_reason) %>% group_split()")
operation <- makeOperation(l.variable, command)
eval(parse(text = operation))


#  Create "minimal tibbles" ("m") of only ... and AS values -------------------
#+ 
#+ Rename the AS columns to identify their alignment of origin
s.variable <- c(
    paste0(l.variable, "[[4]]"),
    paste0(variable[3])
)
m.variable <- paste0("m.", variable)
suffix <- m.variable %>% strsplit(., "\\.") %>% lapply(., `[[`, 3) %>% unlist()
command <- paste0("<- ", s.variable, " %>% dplyr::select(coordinate, tag.AS) %>% dplyr::rename(AS.", suffix, " = tag.AS)")  #TODO read or coordinate
operation <- makeOperation(m.variable, command)
eval(parse(text = operation))

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
m.full$tmp.mm10 <- ifelse(
    is.na(m.full$AS.mm10), "0", "1"
)

m.full$dedup.129S1 <- ifelse(
    is.na(m.full$AS.129S1), "0", "1"
)

m.full$dedup.CAST <- ifelse(
    is.na(m.full$AS.CAST), "0", "1"
)

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
    "<- m.full %>% dplyr::filter_at(vars(assignment), any_vars(. %in% \"",
    lapply(variable, `[[`, 2), "\"))"
)
operation <- makeOperation(lapply(variable, `[[`, 1), command)
eval(parse(text = operation))


#  Make $assignment factors again ---------------------------------------------
command <- paste0("<- ", lapply(variable, `[[`, 1), "$assignment %>% as.factor()")
operation <- makeOperation(paste0(lapply(variable, `[[`, 1), "$assignment"), command)
eval(parse(text = operation))

assignment.NA$assignment <- NA
assignment.NA$assignment <- NA_character_

rm(variable, m.dedup.129S1, m.dedup.CAST, m.dedup.mm10, m.full)

# bak.assignment.129S1 <- assignment.129S1
# bak.assignment.CAST <- assignment.CAST
# bak.assignment.NA <- assignment.NA
# bak.assignment.Neutral <- assignment.Neutral

assignment.129S1 <- bak.assignment.129S1
assignment.CAST <- bak.assignment.CAST
assignment.NA <- bak.assignment.NA
assignment.Neutral <- bak.assignment.Neutral


#  Full joins for 129S1 and CAST ----------------------------------------------

#  Assignment.129S1
assignment.129S1 <- full_join(dedup.129S1, assignment.129S1, by = "coordinate") %>%
    tidyr::drop_na("AS.129S1") %>% 
    dplyr::select(-dplyr::one_of("AS.129S1")) %>%
    dplyr::rename(
        AS.129S1 = tag.AS,
        XS.129S1 = tag.XS,
        NM.129S1 = tag.NM
    ) %>%
    dplyr::select(-XS.129S1, -NM.129S1) %>%
    dplyr::relocate(AS.mm10, .after = seq)

#  Assignment.CAST
assignment.CAST <- full_join(dedup.CAST, assignment.CAST, by = "coordinate") %>%
    tidyr::drop_na("AS.CAST") %>% 
    dplyr::select(-dplyr::one_of("AS.CAST")) %>%
    dplyr::rename(
        AS.CAST = tag.AS,
        XS.CAST = tag.XS,
        NM.CAST = tag.NM
    ) %>%
    dplyr::select(-XS.CAST, -NM.CAST) %>%
    dplyr::relocate(AS.CAST, .after = AS.129S1)

#  Assignment.Neutral
assignment.Neutral <- full_join(dedup.129S1, assignment.Neutral, by = "coordinate") %>%
    tidyr::drop_na("AS.129S1") %>% 
    dplyr::select(-dplyr::one_of("AS.129S1")) %>%
    dplyr::rename(
        AS.129S1 = tag.AS,
        XS.129S1 = tag.XS,
        NM.129S1 = tag.NM
    ) %>%
    dplyr::select(-XS.129S1, -NM.129S1) %>%
    dplyr::relocate(AS.mm10, .before = AS.129S1) %>%
    dplyr::rename(
        c(
            old_coordinate_129S1 = old_coordinate,
            old_rname_129S1 = old_rname,
            old_pos_129S1 = old_pos,
            old_pos_end_129S1 = old_pos_end,
            old_mpos_129S1 = old_mpos,
            old_mpos_end_129S1 = old_mpos_end
        )
    )

#  Munge assignment.Neutral
assignment.Neutral.CAST <- full_join(dedup.CAST, assignment.Neutral, by = "coordinate") %>%
    tidyr::drop_na("AS.CAST")

assignment.Neutral$old_coordinate_CAST <- assignment.Neutral.CAST$old_coordinate
assignment.Neutral$old_rname_CAST <- assignment.Neutral.CAST$old_rname
assignment.Neutral$old_pos_CAST <- assignment.Neutral.CAST$old_pos
assignment.Neutral$old_pos_end_CAST <- assignment.Neutral.CAST$old_pos_end
assignment.Neutral$old_mpos_CAST <- assignment.Neutral.CAST$old_mpos
assignment.Neutral$old_mpos_end_CAST <- assignment.Neutral.CAST$old_mpos_end

assignment.Neutral <- assignment.Neutral %>%
    dplyr::relocate(old_coordinate_CAST, .before = coordinate) %>% 
    dplyr::relocate(old_rname_CAST, .after = old_rname_129S1) %>%
    dplyr::relocate(old_pos_CAST, .after = old_pos_end_129S1) %>%
    dplyr::relocate(old_pos_end_CAST, .after = old_pos_CAST) %>%
    dplyr::relocate(old_mpos_CAST, .after = old_mpos_end_129S1) %>%
    dplyr::relocate(old_mpos_end_CAST, .after = old_mpos_CAST)

rm(assignment.Neutral.CAST)

#  Assignment.NA
assignment.NA <- full_join(dedup.mm10, assignment.NA, by = "coordinate") %>%
    tidyr::drop_na("AS.mm10") %>% 
    dplyr::select(-dplyr::one_of("AS.mm10")) %>%
    dplyr::rename(
        AS.mm10 = tag.AS,
        XS.mm10 = tag.XS,
        NM.mm10 = tag.NM
    ) %>%
    dplyr::select(-old_coordinate, -XS.mm10, -NM.mm10)


#  Save the "assignment tibbles" as RDS
variable <- paste0("assignment.", c("129S1", "CAST", "NA", "Neutral"))
command <- paste0("saveRDS(", variable, ", file = ", "\"", variable, ".", chromosome, ".rds\")")
eval(parse(text = command))
