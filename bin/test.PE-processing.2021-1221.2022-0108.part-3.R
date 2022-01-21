#!/usr/bin/env Rscript

library(stringr)
library(tidyverse)

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


#  Having run part 2 first, load part-2 .Rdata into environment ---------------
load("test.PE-processing.2021-1221.2022-0106.part-2.Rdata")
#  Doing so loads in appropriate functions, variables, etc.


#  Script name ----------------------------------------------------------------
script <- "test.PE-processing.2021-1221.2022-0108.part-3.R"


# -------------------------------------
# STOP, GO TO SHORTCUT
# -------------------------------------
#  Perform liftOvers of the 129S1 and CAST tibbles ----------------------------

#  Create and write out 129S1 and CAST .bed file(s) for liftOver
chromosome <- "chrX"

writeBedOperation <- function(variable_bed) {
    paste0(
        "readr::write_tsv(",
            variable_bed, ", ",
            "file = \"", variable_bed, "\", ",
            "col_names = FALSE",
        ")"
    )
}

#  pos
variable_bed <- paste0(variable_tbl, ".pos.", chromosome, ".bed")

command <- paste0(
    "<- tibble::tibble(",
        "\"rname\" = ", variable_tbl, "$rname, ",
        "\"pos\" = ", variable_tbl, "$pos, ",
        "\"pos_end\" = ", variable_tbl, "$pos_end, ",
        "\"criteria\" = ", variable_tbl, "$criteria",
    ")"
)
operation <- makeOperation(variable_bed, command)
evaluateOperation(operation)

operation <- writeBedOperation(variable_bed)
evaluateOperation(operation)

#  mpos
variable_bed <- paste0(variable_tbl, ".mpos.", chromosome, ".bed")

command <- paste0(
    "<- tibble::tibble(",
        "\"mrnm\" = ", variable_tbl, "$mrnm, ",
        "\"mpos\" = ", variable_tbl, "$mpos, ",
        "\"mpos_end\" = ", variable_tbl, "$mpos_end, ",
        "\"criteria\" = ", variable_tbl, "$criteria",
    ")"
)
operation <- makeOperation(variable_bed, command)
evaluateOperation(operation)

operation <- writeBedOperation(variable_bed)
evaluateOperation(operation)


#  Run shell code to perform the liftOvers ------------------------------------
mate <- c("pos", "mpos")

#  129S1 ------------------------------
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

from_R=$(find . -maxdepth 1 -type f -name \"tbl.129S1.", mate, ".", chromosome, ".bed\" -printf '%f\n')
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
for (i in 1:length(shell_code)) {
    system(shell_code[i])
}

#  CAST -------------------------------
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

from_R=$(find . -maxdepth 1 -type f -name \"tbl.CAST.", mate, ".", chromosome, ".bed\" -printf '%f\n')
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
for (i in 1:length(shell_code)) {
    system(shell_code[i])
}


# -------------------------------------
# SHORTCUT
# -------------------------------------
#  Read in the "lifted" data --------------------------------------------------
file <- list.files(pattern = paste0("\\.", chromosome,".lifted.bed$")) %>%
    stringr::str_subset("^tbl\\.")
variable_lifted <- file
mapply(
    assign, variable_lifted, file, MoreArgs = list(envir = parent.frame())
)

#  pos
command <- paste0(
    "<- readr::read_tsv(",
        "file = \"", stringr::str_subset(file, "\\.pos\\."), "\", ",
        "col_names = c(",
            "\"liftOver_rname\", ",
            "\"liftOver_pos\", ",
            "\"liftOver_pos_end\", ",
            "\"criteria\"",
        ")",
    ")"
)
operation <- makeOperation(
    stringr::str_subset(variable_lifted, "\\.pos\\."), command
)
evaluateOperation(operation)

#  mpos
command <- paste0(
    "<- readr::read_tsv(",
        "file = \"", stringr::str_subset(file, "\\.mpos\\."), "\", ",
        "col_names = c(",
            "\"liftOver_mrnm\", ",
            "\"liftOver_mpos\", ",
            "\"liftOver_mpos_end\", ",
            "\"criteria\"",
        ")",
    ")"
)
operation <- makeOperation(
    stringr::str_subset(variable_lifted, "\\.mpos\\."), command
)
evaluateOperation(operation)

#  Note liftOver status
command <- paste0("<- \"#liftOver successful\"")
operation <- makeOperation(
    paste0(variable_lifted, "$liftOver_reason"), command
)
evaluateOperation(operation)


variable_bed_mpos <- variable_bed
variable_bed_pos <- variable_bed %>% stringr::str_replace("mpos", "pos")

operation <- paste0(
    "rm(", variable_bed_mpos, ", ", variable_bed_pos, ")"
)
evaluateOperation(operation)


# -------------------------------------
# STOP, GO TO SHORTCUT
# -------------------------------------
#  Munge the "unlifted" data --------------------------------------------------
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
done < <(find . -maxdepth 1 -type f -name \"tbl.*.", mate, ".", chromosome, ".bed\" -print0)

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
for (i in 1:length(shell_code)) {
    system(shell_code[i])
}

#  Don't mind the shell errors for tbl.mm10.*
system("rm tbl.mm10.pos.chrX.unlifted.bed tbl.mm10.mpos.chrX.unlifted.bed")


# -------------------------------------
# SHORTCUT
# -------------------------------------
#  Read in the "unlifted" data ------------------------------------------------
file <- list.files(pattern = paste0("\\.", chromosome,".unlifted.bed$")) %>%
    stringr::str_subset("^tbl\\.")
variable_unlifted <- file
mapply(
    assign, variable_unlifted, file, MoreArgs = list(envir = parent.frame())
)

#  pos
command <- paste0(
    "<- readr::read_tsv(",
        "file = \"", stringr::str_subset(file, "\\.pos\\."), "\", ",
        "col_names = c(",
            "\"liftOver_reason\", ",
            "\"liftOver_rname\", ",
            "\"liftOver_pos\", ",
            "\"liftOver_pos_end\", ",
            "\"criteria\"",
        ")",
    ")"
)
operation <- makeOperation(
    stringr::str_subset(variable_unlifted, "\\.pos\\."), command
)
evaluateOperation(operation)

#  mpos
command <- paste0(
    "<- readr::read_tsv(",
        "file = \"", stringr::str_subset(file, "\\.mpos\\."), "\", ",
        "col_names = c(",
            "\"liftOver_reason\", ",
            "\"liftOver_mrnm\", ",
            "\"liftOver_mpos\", ",
            "\"liftOver_mpos_end\", ",
            "\"criteria\"",
        ")",
    ")"
)
operation <- makeOperation(
    stringr::str_subset(variable_unlifted, "\\.mpos\\."), command
)
evaluateOperation(operation)


#  Munge the information in the "unlifted" 129S1 and CAST variables -----------
# rm(
#     tbl.129S1.mpos.chrX.unlifted.bed,
#     tbl.129S1.pos.chrX.unlifted.bed,
#     tbl.CAST.mpos.chrX.unlifted.bed,
#     tbl.CAST.pos.chrX.unlifted.bed,
# )

#  pos
command <- paste0("<- NA_character_")
operation <- makeOperation(
    paste0(stringr::str_subset(variable_unlifted, "\\.pos\\."), "$liftOver_rname"), command
)
evaluateOperation(operation)

command <- paste0("<- NA_integer_")
operation <- makeOperation(
    paste0(stringr::str_subset(variable_unlifted, "\\.pos\\."), "$liftOver_pos"), command
)
evaluateOperation(operation)

command <- paste0("<- NA_integer_")
operation <- makeOperation(
    paste0(stringr::str_subset(variable_unlifted, "\\.pos\\."), "$liftOver_pos_end"), command
)
evaluateOperation(operation)

#  mpos
command <- paste0("<- NA_character_")
operation <- makeOperation(
    paste0(stringr::str_subset(variable_unlifted, "\\.mpos\\."), "$liftOver_mrnm"), command
)
evaluateOperation(operation)

command <- paste0("<- NA_integer_")
operation <- makeOperation(
    paste0(stringr::str_subset(variable_unlifted, "\\.mpos\\."), "$liftOver_mpos"), command
)
evaluateOperation(operation)

command <- paste0("<- NA_integer_")
operation <- makeOperation(
    paste0(stringr::str_subset(variable_unlifted, "\\.mpos\\."), "$liftOver_mpos_end"), command
)
evaluateOperation(operation)


#  Temporary: save.image() ----------------------------------------------------
# save.image("test.Rdata")
# load("test.Rdata")


#  Perform a full join of the unlifted and lifted .bed files ------------------
variable_bed <- c(
    "tbl.129S1.mpos.chrX.bed",
    "tbl.129S1.pos.chrX.bed",
    "tbl.CAST.mpos.chrX.bed",
    "tbl.CAST.pos.chrX.bed"
)
command <- paste0(
    "<- ", "dplyr::bind_rows(", variable_unlifted, ", ", variable_lifted, ")"
)
operation <- makeOperation(variable_bed, command)
evaluateOperation(operation)

#  Remove unneeded .bed files (cleaning up environment)
operation <- paste0("rm(", variable_lifted, ", ", variable_unlifted, ")")
evaluateOperation(operation)


#  Perform a full join of liftOver information to main tibbles ----------------

#  Munge files for full join of liftOver information
command <- paste0(
    "<- ", stringr::str_subset(variable_bed_mpos, "mm10", negate = TRUE), " %>% ",
        "dplyr::select(-liftOver_reason)"
)
operation <- makeOperation(
    stringr::str_subset(variable_bed_mpos, "mm10", negate = TRUE), command
)
evaluateOperation(operation)

#  Full join of liftOver information
command <- paste0(
    "<- ", "dplyr::full_join(",
        variable_tbl[1:2], ", ",
        variable_bed_pos[1:2], ", ",
        "by = \"criteria\"",
    ")", " %>% ",
        "dplyr::full_join(",
            ".", ", ",
            variable_bed_mpos[1:2], ", ",
            "by = \"criteria\"",
        ")", " %>% ",
        "dplyr::relocate(liftOver_mpos, .after = liftOver_pos) %>%",
        "dplyr::relocate(liftOver_mpos_end, .after = liftOver_pos_end) %>%",
        "dplyr::relocate(liftOver_rname, .before = rname) %>%",
        "dplyr::relocate(liftOver_mrnm, .after = liftOver_rname) %>%",
        "dplyr::relocate(c(liftOver_pos:liftOver_mpos_end), .before = pos) %>%",
        "dplyr::relocate(liftOver_reason, .before = mate_status)"
)
operation <- makeOperation(variable_tbl[1:2], command)
evaluateOperation(operation)

#  Factorize $liftOver_reason
command <- paste0(
    "<- ", stringr::str_subset(variable_tbl, "mm10", negate = TRUE), "$liftOver_reason %>% ",
        "as_factor()"
)
operation <- makeOperation(
    paste0(stringr::str_subset(variable_tbl, "mm10", negate = TRUE), "$liftOver_reason"),
    command
)
evaluateOperation(operation)

#  Remove "#" from $liftOver_reason levels
tbl.129S1$liftOver_reason <- plyr::revalue(
    tbl.129S1$liftOver_reason,
    c(
        "#Deleted in new" = "Deleted in new",
        "#Partially deleted in new" = "Partially deleted in new",
        "#Split in new" = "Split in new",
        "#liftOver successful" = "liftOver successful"
    )
)

tbl.CAST$liftOver_reason <- plyr::revalue(
    tbl.CAST$liftOver_reason,
    c(
        "#Deleted in new" = "Deleted in new",
        "#Partially deleted in new" = "Partially deleted in new",
        "#Split in new" = "Split in new",
        "#liftOver successful" = "liftOver successful"
    )
)

#  For easier reading, replace "liftOver_*" with "lO_*" in column names
operation <- paste0(
    "colnames(", variable_tbl[1:2], ")", " <- ",
        "gsub(",
            "\"liftOver_\", \"lO_\", colnames(", variable_tbl[1:2], ")",
        ")"
)
evaluateOperation(operation)


#  Check on NAs and duplicates... ---------------------------------------------
test <- tbl.129S1[is.na(tbl.129S1$lO_reason), ]
test.dup <- test %>%
    dplyr::filter(criteria %in% unique(.[["criteria"]][
        duplicated(.[["criteria"]])
    ]))  # 0

test <- tbl.CAST[is.na(tbl.CAST$lO_reason), ]
test.dup <- test %>%
    dplyr::filter(criteria %in% unique(.[["criteria"]][
        duplicated(.[["criteria"]])
    ]))  # 0

rm(test, test.dup)
#NOTE Bug arises from running the '#  Munge the "unlifted" data' code chunk more than once
#TODO If not 0, then stop the program and print a warning


#  Clean up both CAST and 129S1 -----------------------------------------------
rm(
    chromosome,
    command,
    file,
    i,
    mate,
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
    variable_bed_mpos,
    variable_bed_pos,
    variable_lifted,
    variable_tbl,
    variable_uniline,
    variable_unlifted
)

rm(operation)

# system("rm *.bed")

# save.image("test.Rdata")
# load("test.Rdata")


#  Sort the tibbles while maintaining mate pairs ------------------------------
variable_tbl <- paste0("tbl.", c("129S1", "CAST"))

for (i in 1:length(variable_tbl)) {
    #  Sort the tibble of interest by pos while maintaining proper mate pairs
    df <- eval(parse(text = variable_tbl[i]))
    
    #  Create tibbles from odd and even rows
    odd.seq <- seq(1, nrow(df), 2)
    even.seq <- (seq(1, nrow(df), 2)) + 1
    
    odd <- df[odd.seq, ]
    even <- df[even.seq, ]
    
    #  Create a tibble with two mates per row (joined by "groupid_3"); then,
    #+ sort by pos.x (the odd-row pos value)
    df <- full_join(odd, even, by = "groupid") %>%
        dplyr::arrange(lO_pos.x) %>%  # Sort by pos.x
        dplyr::rename(groupid.x = groupid)  # new_name = old_name
    df$groupid.y <- df$groupid.x
    df <- df %>% dplyr::relocate(groupid.y, .after = qmpos.y)
    
    #  Rename column names to denote odd/even status
    colnames(df) <- str_replace_all(colnames(df), "\\.x", "\\.odd")
    colnames(df) <- str_replace_all(colnames(df), "\\.y", "\\.even")
    
    #  Split the tibble by odd or even status
    odd <- df[stringr::str_subset(colnames(df), "\\.odd")]
    even <- df[stringr::str_subset(colnames(df), "\\.even")]
    
    #  Strip suffixes from column names
    colnames(odd) <- str_replace_all(colnames(odd), "\\.odd", "")
    colnames(even) <- str_replace_all(colnames(even), "\\.even", "")
    
    #  Odd tibble is tibble #1, even tibble is tibble #2
    odd$tibble <- "1"
    even$tibble <- "2"
    
    #  Again, interleave the odd/even rows, then arrange them by group ID and
    #+ tibble number
    command <- paste0(
        "<- odd %>%",
            "dplyr::mutate(groupid = row_number()) %>% ",
            "dplyr::bind_rows(even %>% mutate(groupid = row_number())) %>% ",
            "dplyr::arrange(groupid, tibble) %>% ",
            "dplyr::relocate(groupid, .after = qmpos) %>% ",
            "dplyr::select(-tibble)"
    )
    operation <- makeOperation(variable_tbl[i], command)
    evaluateOperation(operation)
    
    #  Again, check to make sure that there are no more than two entries per
    #+ group ID
    n_occur <- data.frame(table(eval(parse(text = variable_tbl[i]))$groupid))
    print(n_occur[n_occur$Freq > 2, ])
    rm(n_occur)
}

rm(df, i, even, even.seq, odd, odd.seq, variable_tbl)


#  Save environment image -----------------------------------------------------
save.image(stringr::str_remove(script, ".R") %>% paste0(., ".Rdata"))


#  Script is completed; clean up environment ----------------------------------
rm(list = ls())

#  Now, go to part 4
