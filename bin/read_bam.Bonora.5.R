#!/usr/bin/env Rscript

#  Load libraries
library(BSgenome)
library(dplyr)
library(forcats)
library(ggplot2)
library(magrittr)
library(parallel)
library(purrr)
library(Rsamtools)
library(rtracklayer)
library(scales)

options(pillar.sigfig = 8, scipen = 10000)

setwd("/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora")


# -----------------------------------------------------------------------------
# rm(
#     dedup.129S1,
#     dedup.CAST,
#     dedup.mm10
# )
# 
# rm(
#     chromosomes,
#     command_pipe,
#     command,
#     file,
#     i,
#     int,
#     map_info,
#     mungeTmpGB,
#     n,
#     operation,
#     prefix_assembly,
#     prefix_GB,
#     prefix_KA,
#     prefix_tmp.df,
#     shell_code_129S1,
#     shell_code_CAST,
#     shell_code,
#     spin,
#     suffix_GB,
#     suffix_in_GB,
#     suffix_not_in_GB,
#     tag_info,
#     tmp.CAST,
#     tmp.CAST.bed,
#     tmp.CAST.lifted.bed,
#     tmp.CAST.mpos.bed,
#     tmp.CAST.mpos.bed,
#     tmp.CAST.mpos.lifted.bed,
#     tmp.CAST.mpos.unlifted.bed,
#     tmp.CAST.pos.bed,
#     tmp.CAST.pos.bed,
#     tmp.CAST.pos.lifted.bed,
#     tmp.CAST.pos.unlifted.bed,
#     tmp.CAST.unlifted.bed,
#     tmp.df.129S1.all_dups,
#     tmp.df.129S1.each_dup,
#     tmp.df.129S1.in_all,
#     tmp.df.129S1.in_any,
#     tmp.df.129S1.not_in_GB,
#     tmp.df.129S1.in_alt.CAST,
#     tmp.df.129S1.in_ambig,
#     tmp.df.129S1.in_contra,
#     tmp.df.129S1.in_ref.129S1,
#     tmp.df.CAST.all_dups,
#     tmp.df.CAST.each_dup,
#     tmp.df.CAST.in_all,
#     tmp.df.CAST.in_any,
#     tmp.df.CAST.not_in_GB,
#     tmp.df.CAST.in_alt.CAST,
#     tmp.df.CAST.in_ambig,
#     tmp.df.CAST.in_contra,
#     tmp.df.CAST.in_ref.129S1,
#     tmp.df.CAST.not_in_alt.CAST,
#     tmp.df.CAST.not_in_ambig,
#     tmp.df.CAST.not_in_contra,
#     tmp.df.CAST.not_in_ref.129S1,
#     tmp.df.in_alt,
#     tmp.df.not_in_alt,
#     tmp.final.129S1,
#     tmp.final.CAST,
#     tmp.join,
#     tmp.join.129S1,
#     tmp.join.CAST,
#     tmp.join.CAST.not_in_GB,
#     tmp.method.GB,
#     tmp.method.GB.129S1.alt.CAST,
#     tmp.method.GB.129S1.ambig,
#     tmp.method.GB.129S1.contra,
#     tmp.method.GB.129S1.ref.129S1,
#     tmp.method.GB.CAST.alt.CAST,
#     tmp.method.GB.CAST.ambig,
#     tmp.method.GB.CAST.contra,
#     tmp.method.GB.CAST.ref.129S1,
#     tmp.method.KA,
#     tmp.method.KA.129S1.alt.CAST,
#     tmp.method.KA.129S1.ambig,
#     tmp.method.KA.129S1.contra,
#     tmp.method.KA.129S1.ref.129S1,
#     tmp.method.KA.CAST.alt.CAST,
#     tmp.method.KA.CAST.ambig,
#     tmp.method.KA.CAST.contra,
#     tmp.method.KA.CAST.ref.129S1,
#     tmp.rbind,
#     variable_all,
#     variable_any,
#     variable_assignment,
#     variable_column,
#     variable_join_1,
#     variable_join_2,
#     variable_tibble,
#     variable_rname,
#     variable
# )

# system("rm *.bed")

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
            b_s,
            GB_rname,
            GB_strand,
            GB_pos,
            GB_mrnm,
            GB_mpos
        )
}


# -----------------------------------------------------------------------------
#  ...need to have run read_bam.Bonora.3_4.R first


# -------
#  Load in data
# n <- 50000
# tmp.CAST <- assignment.CAST[1:n, ]
# tmp.129S1 <- assignment.129S1[1:n, ]

tmp.CAST <- assignment.CAST
tmp.129S1 <- assignment.129S1


# -------
#  CAST
tmp.CAST$in_alt_CAST <- ifelse(
    tmp.CAST$b_s %in% GB.alt.CAST$b_s, "1", "0"  #TODO b_f or b_s
)
tmp.CAST$in_ref_129S1 <- ifelse(
    tmp.CAST$b_s %in% GB.ref.129S1$b_s, "1", "0"  #TODO b_f or b_s
)
tmp.CAST$in_ambig <- ifelse(
    tmp.CAST$b_s %in% GB.ambig$b_s, "1", "0"  #TODO b_f or b_s
)
tmp.CAST$in_contra <- ifelse(
    tmp.CAST$b_s %in% GB.contra$b_s, "1", "0"  #TODO b_f or b_s
)
tmp.CAST$intersection_GB <- paste0(
    tmp.CAST$in_alt_CAST,
    tmp.CAST$in_ref_129S1,
    tmp.CAST$in_ambig,
    tmp.CAST$in_contra
)
tmp.CAST$intersection_GB <- tmp.CAST$intersection_GB %>% as.factor()

# --
#  129S1
tmp.129S1$in_alt_CAST <- ifelse(
    tmp.129S1$b_s %in% GB.alt.CAST$b_s, "1", "0"  #TODO b_f or b_s
)
tmp.129S1$in_ref_129S1 <- ifelse(
    tmp.129S1$b_s %in% GB.ref.129S1$b_s, "1", "0"  #TODO b_f or b_s
)
tmp.129S1$in_ambig <- ifelse(
    tmp.129S1$b_s %in% GB.ambig$b_s, "1", "0"  #TODO b_f or b_s
)
tmp.129S1$in_contra <- ifelse(
    tmp.129S1$b_s %in% GB.contra$b_s, "1", "0"  #TODO b_f or b_s
)
tmp.129S1$intersection_GB <- paste0(
    tmp.129S1$in_alt_CAST,
    tmp.129S1$in_ref_129S1,
    tmp.129S1$in_ambig,
    tmp.129S1$in_contra
)
tmp.129S1$intersection_GB <- tmp.129S1$intersection_GB %>% as.factor()


# -------
#  CAST
tmp.CAST.pos.bed <- tibble::tibble(
    "rname" = tmp.CAST$rname,
    "pos" = tmp.CAST$pos,
    "pos_end" = tmp.CAST$pos_end,
    "b_s" = tmp.CAST$b_s  #TODO b_f or b_s
)

# tmp.CAST.mpos.bed <- tibble::tibble(
#     "mrnm" = tmp.CAST$mrnm,
#     "mpos" = tmp.CAST$mpos,
#     "mpos_end" = tmp.CAST$mpos_end,
#     "b_s" = tmp.CAST$b_s  #TODO b_f or b_s
# )

readr::write_tsv(
    tmp.CAST.pos.bed,
    file = "tmp.CAST.pos.bed",
    col_names = FALSE
)

# readr::write_tsv(
#     tmp.CAST.mpos.bed,
#     file = "tmp.CAST.mpos.bed",
#     col_names = FALSE
# )

# --
#  129S1
tmp.129S1.pos.bed <- tibble::tibble(
    "rname" = tmp.129S1$rname,
    "pos" = tmp.129S1$pos,
    "pos_end" = tmp.129S1$pos_end,
    "b_s" = tmp.129S1$b_s  #TODO b_f or b_s
)

# tmp.129S1.mpos.bed <- tibble::tibble(
#     "mrnm" = tmp.129S1$mrnm,
#     "mpos" = tmp.129S1$mpos,
#     "mpos_end" = tmp.129S1$mpos_end,
#     "b_s" = tmp.129S1$b_s  #TODO b_f or b_s
# )

readr::write_tsv(
    tmp.129S1.pos.bed,
    file = "tmp.129S1.pos.bed",
    col_names = FALSE
)

# readr::write_tsv(
#     tmp.129S1.mpos.bed,
#     file = "tmp.129S1.mpos.bed",
#     col_names = FALSE
# )


# -------
#  CAST
shell_code <- 
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

prefixes_from_R=\"tmp.CAST.pos\"
# prefixes_from_R=(
#     tmp.CAST.pos
#     tmp.CAST.mpos
# )

# -------
#  Run liftOver loop
for file_prefix in \"${prefixes_from_R[@]}\"; do
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
system(shell_code)

# --
#  129S1
shell_code <- 
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

prefixes_from_R=\"tmp.129S1.pos\"
# prefixes_from_R=(
#     tmp.129S1.pos
#     tmp.129S1.mpos
# )

# -------
#  Run liftOver loop
for file_prefix in \"${prefixes_from_R[@]}\"; do
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
system(shell_code)


# -------
#  CAST and 129S1
file <- list.files(pattern = "\\.lifted.bed$")
variable <- file
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

command <- paste0(
    "<- readr::read_tsv(file = \"", file, "\", col_names = c(\"liftOver_rname\", \"liftOver_pos\", \"liftOver_pos_end\", \"b_s\"))"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))

shell_code <-
"
# -------
#  Start recording time
start=$(date +%s)

# -------
#  Set up files
prefixes_from_R=(
    tmp.CAST.pos
    tmp.129S1.pos
)
# prefixes_from_R=(
#     tmp.CAST.pos
#     tmp.CAST.mpos
#     tmp.129S1.pos
#     tmp.129S1.mpos
# )

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
system(shell_code)


# -------
#  CAST and 129S1
file <- list.files(pattern = "\\.unlifted.bed$")
variable <- file[1:length(file)]
mapply(
    assign, variable, file, MoreArgs = list(envir = parent.frame())
)

command <- paste0(
    "<- readr::read_tsv(file = \"", file, "\", col_names = c(\"liftOver_reason\", \"liftOver_rname\", \"liftOver_pos\", \"liftOver_pos_end\", \"b_s\"))"
)
operation <- makeOperation(variable, command)
eval(parse(text = operation))

# --
#  CAST
tmp.CAST.pos.unlifted.bed$liftOver_rname <- NA_character_
tmp.CAST.pos.unlifted.bed$liftOver_pos <- NA_integer_
tmp.CAST.pos.unlifted.bed$liftOver_pos_end <- NA_integer_

tmp.join.CAST <- dplyr::bind_rows(tmp.CAST.pos.unlifted.bed, tmp.CAST.pos.lifted.bed) %>%
    full_join(., tmp.CAST, by = "b_s") %>%
    select(-liftOver_pos_end) %>%
    dplyr::arrange("b_s") 

tmp.join.CAST$liftOver_reason[is.na(tmp.join.CAST$liftOver_reason)] <- "#liftOver successful"
tmp.join.CAST$liftOver_reason <- tmp.join.CAST$liftOver_reason %>% as_factor()
tmp.join.CAST <- tmp.join.CAST %>% dplyr::relocate(c(liftOver_rname, liftOver_pos, liftOver_reason), .after = pos_end)

# --
#  129S1
tmp.129S1.pos.unlifted.bed$liftOver_rname <- NA_character_
tmp.129S1.pos.unlifted.bed$liftOver_pos <- NA_integer_
tmp.129S1.pos.unlifted.bed$liftOver_pos_end <- NA_integer_

tmp.join.129S1 <- dplyr::bind_rows(tmp.129S1.pos.unlifted.bed, tmp.129S1.pos.lifted.bed) %>%
    full_join(., tmp.129S1, by = "b_s") %>%
    select(-liftOver_pos_end) %>%
    dplyr::arrange("b_s") 

tmp.join.129S1$liftOver_reason[is.na(tmp.join.129S1$liftOver_reason)] <- "#liftOver successful"
tmp.join.129S1$liftOver_reason <- tmp.join.129S1$liftOver_reason %>% as_factor()
tmp.join.129S1 <- tmp.join.129S1 %>% dplyr::relocate(c(liftOver_rname, liftOver_pos, liftOver_reason), .after = pos_end)


# -------
#  Clean up both CAST and 129S1
rm(
    file,
    n,
    shell_code,
    tmp.129S1,
    tmp.129S1.mpos.bed,
    tmp.129S1.mpos.lifted.bed,
    tmp.129S1.mpos.unlifted.bed,
    tmp.129S1.pos.bed,
    tmp.129S1.pos.lifted.bed,
    tmp.129S1.pos.unlifted.bed,
    tmp.CAST,
    tmp.CAST.mpos.bed,
    tmp.CAST.mpos.lifted.bed,
    tmp.CAST.mpos.unlifted.bed,
    tmp.CAST.pos.bed,
    tmp.CAST.pos.lifted.bed,
    tmp.CAST.pos.unlifted.bed
)

system("rm *.bed")


# -------
#  CAST
tmp.method.KA.CAST.alt.CAST <- tmp.join.CAST[tmp.join.CAST$b_s %in% GB.alt.CAST$b_s, ]
tmp.method.KA.CAST.ref.129S1 <- tmp.join.CAST[tmp.join.CAST$b_s %in% GB.ref.129S1$b_s, ]
tmp.method.KA.CAST.ambig <- tmp.join.CAST[tmp.join.CAST$b_s %in% GB.ambig$b_s, ]
tmp.method.KA.CAST.contra <- tmp.join.CAST[tmp.join.CAST$b_s %in% GB.contra$b_s, ]

tmp.method.GB.CAST.alt.CAST <- GB.alt.CAST[GB.alt.CAST$b_s %in% tmp.join.CAST$b_s, ]
tmp.method.GB.CAST.ref.129S1 <- GB.ref.129S1[GB.ref.129S1$b_s %in% tmp.join.CAST$b_s, ]
tmp.method.GB.CAST.ambig <- GB.ambig[GB.ambig$b_s %in% tmp.join.CAST$b_s, ]
tmp.method.GB.CAST.contra <- GB.contra[GB.contra$b_s %in% tmp.join.CAST$b_s, ]

# --
#  129S1
tmp.method.KA.129S1.alt.CAST <- tmp.join.129S1[tmp.join.129S1$b_s %in% GB.alt.CAST$b_s, ]
tmp.method.KA.129S1.ref.129S1 <- tmp.join.129S1[tmp.join.129S1$b_s %in% GB.ref.129S1$b_s, ]
tmp.method.KA.129S1.ambig <- tmp.join.129S1[tmp.join.129S1$b_s %in% GB.ambig$b_s, ]
tmp.method.KA.129S1.contra <- tmp.join.129S1[tmp.join.129S1$b_s %in% GB.contra$b_s, ]

tmp.method.GB.129S1.alt.CAST <- GB.alt.CAST[GB.alt.CAST$b_s %in% tmp.join.129S1$b_s, ]
tmp.method.GB.129S1.ref.129S1 <- GB.ref.129S1[GB.ref.129S1$b_s %in% tmp.join.129S1$b_s, ]
tmp.method.GB.129S1.ambig <- GB.ambig[GB.ambig$b_s %in% tmp.join.129S1$b_s, ]
tmp.method.GB.129S1.contra <- GB.contra[GB.contra$b_s %in% tmp.join.129S1$b_s, ]


# -----
#  CAST and 129S1
variable <- c(
    "tmp.method.GB.CAST.alt.CAST",
    "tmp.method.GB.CAST.ref.129S1",
    "tmp.method.GB.CAST.ambig",
    "tmp.method.GB.CAST.contra",
    "tmp.method.GB.129S1.alt.CAST",
    "tmp.method.GB.129S1.ref.129S1",
    "tmp.method.GB.129S1.ambig",
    "tmp.method.GB.129S1.contra"
)
command_pipe <- paste0("<- ", variable, " %>% ")
command <- paste0(command_pipe, "mungeTmpGB()")
operation <- makeOperation(variable, command)
eval(parse(text = operation))

rm(
    variable,
    command_pipe,
    command,
    operation
)


# -------
#  Analyses to determine if liftOver_rname == GB_rname, liftOver_pos == GB_pos
prefix_KA <- "tmp.method.KA"
prefix_GB <- "tmp.method.GB"
prefix_tmp.df <- "tmp.df"
prefix_assembly <- c("CAST", "129S1")
suffix_GB <- c(
    "alt.CAST",
    "ref.129S1",
    "ambig",
    "contra"
)
suffix_in_GB <- c(paste0("in_", suffix_GB))
suffix_not_in_GB <- c(paste0("not_", suffix_in_GB))

# variable_assignment <- paste0(prefix_tmp.df, ".", prefix_assembly, ".", suffix_in_GB)
# variable_join_1 <- paste0(prefix_KA, ".", prefix_assembly, ".", suffix_GB)
# variable_join_2 <- paste0(prefix_GB, ".", prefix_assembly, ".", suffix_GB)

#  See stackoverflow.com/questions/16143700/pasting-two-vectors-with-combinations-of-all-vectors-elements
variable_assignment <- apply(expand.grid(prefix_tmp.df, prefix_assembly, suffix_in_GB), 1, paste0, collapse = ".")
variable_join_1 <- apply(expand.grid(prefix_KA, prefix_assembly, suffix_GB), 1, paste0, collapse = ".")
variable_join_2 <- apply(expand.grid(prefix_GB, prefix_assembly, suffix_GB), 1, paste0, collapse = ".")

command <- paste0("<- dplyr::full_join(", variable_join_1, ", ", variable_join_2, ", ", "by = \"b_s\") %>% dplyr::arrange(\"b_s\")")  #TODO b_f or b_s
operation <- makeOperation(variable_assignment, command)
eval(parse(text = operation))

command <- paste0("<- (", variable_assignment, "$liftOver_rname == ", variable_assignment, "$GB_rname)")
operation <- makeOperation(paste0(variable_assignment,"$equal_rname_GB_liftOver"), command)
eval(parse(text = operation))

command <- paste0("<- (", variable_assignment, "$liftOver_pos == ", variable_assignment, "$GB_pos)")
operation <- makeOperation(paste0(variable_assignment, "$equal_pos_GB_liftOver"), command)
eval(parse(text = operation))

rm(
    prefix_GB,
    prefix_KA,
    tmp.method.GB.CAST.alt.CAST,
    tmp.method.GB.CAST.ref.129S1,
    tmp.method.GB.CAST.ambig,
    tmp.method.GB.CAST.contra,    
    tmp.method.KA.CAST.alt.CAST,
    tmp.method.KA.CAST.ref.129S1,
    tmp.method.KA.CAST.ambig,
    tmp.method.KA.CAST.contra
)


# -------
#  Row-binding and munging...
variable_all <- c("tmp.df.CAST.in_all", "tmp.df.129S1.in_all")
variable_any <- c("tmp.df.CAST.in_any", "tmp.df.129S1.in_any")

# --
#  CAST
command <- paste0("<- rbind(", variable_assignment[1], ", ", variable_assignment[3], ", ", variable_assignment[5], ", ", variable_assignment[7], ")")
operation <- makeOperation(variable_all[1], command)
eval(parse(text = operation))

# --
#  129S1
command <- paste0("<- rbind(", variable_assignment[2], ", ", variable_assignment[4], ", ", variable_assignment[6], ", ", variable_assignment[8], ")")
operation <- makeOperation(variable_all[2], command)
eval(parse(text = operation))


# -------
#  Filter out all rows without unique b_s
#+ CAST and 129S1
command <- paste0("<- ", variable_all, " %>% dplyr::arrange(rname, pos, flag) %>% dplyr::distinct(b_s, .keep_all = TRUE)")  #TODO b_f or b_s
operation <- makeOperation(variable_any, command)
eval(parse(text = operation))

#  Test: Get one of each duplicated element
#+ CAST and 129S1
command <- paste0("<- ", variable_all, "%>% dplyr::filter(duplicated(.[[\"b_s\"]]))")
operation <- makeOperation(paste0("tmp.df.", prefix_assembly,".each_dup"), command)
eval(parse(text = operation))

#  Test: Get all of the duplicated elements
#+ CAST and 129S1
command <- paste0("<- ", variable_all, "%>% dplyr::filter(b_s %in% unique(.[[\"b_s\"]][duplicated(.[[\"b_s\"]])]))")
operation <- makeOperation(paste0("tmp.df.", prefix_assembly,".all_dups"), command)
eval(parse(text = operation))


# -------
#  Create a tibble of rows not found in the following: in_alt.CAST,
#+ in_ref.129S1, in_ambig, in_contra

# --
#  CAST
tmp.df.CAST.not_in_GB <- dplyr::anti_join(tmp.join.CAST, tmp.df.CAST.in_any, by = "b_s")  # [1] 36008
nrow(tmp.df.CAST.in_all) + nrow(tmp.df.CAST.not_in_GB)  # [1] 124145
nrow(tmp.df.CAST.in_any) + nrow(tmp.df.CAST.not_in_GB)  # [1] 122045

# --
#  129S1
tmp.df.129S1.not_in_GB <- dplyr::anti_join(tmp.join.129S1, tmp.df.129S1.in_any, by = "b_s")  # [1] 28708
nrow(tmp.df.129S1.in_all) + nrow(tmp.df.129S1.not_in_GB)  # [1] 123629
nrow(tmp.df.129S1.in_any) + nrow(tmp.df.129S1.not_in_GB)  # [1] 121090


# -------
tmp.final.CAST <- dplyr::bind_rows(tmp.df.CAST.in_any, tmp.df.CAST.not_in_GB)
tmp.final.129S1 <- dplyr::bind_rows(tmp.df.129S1.in_any, tmp.df.129S1.not_in_GB)

variable <- c("tmp.final.CAST", "tmp.final.129S1")
command <- paste0("<- (", variable, "$rname == ", variable, "$liftOver_rname)")
operation <- makeOperation(paste0(variable, "$liftOver_rname_equal"), command)
eval(parse(text = operation))

command <- paste0("<- ", variable, " %>% dplyr::relocate(liftOver_rname_equal, .after = liftOver_reason)")
operation <- makeOperation(variable, command)
eval(parse(text = operation))

saveRDS(tmp.final.CAST, file = "tmp.final.CAST.rds", compress = FALSE)
saveRDS(tmp.final.129S1, file = "tmp.final.129S1.rds", compress = FALSE)

View(tmp.final.CAST)
View(tmp.final.129S1)


# -------
#  Clean up
rm(
    command,
    operation,
    prefix_assembly,
    prefix_tmp.df,
    suffix_GB,
    suffix_in_GB,
    suffix_not_in_GB,
    variable,
    variable_all,
    variable_any,
    variable_assignment,
    variable_join_1,
    variable_join_2
)
