
#  - https://lokraj.me/post/orf-genome/
#+ - https://stackoverflow.com/questions/55248624/passing-variable-in-r-system-command
#+ - https://stackoverflow.com/questions/12226846/count-letter-differences-of-two-strings

library(AnnotationDbi)
library(BSgenome)
library(GenomicFeatures)


# -----------------------------------------------------------------------------
writeSeed <- function(
    package,
    title,
    description,
    version,
    author,
    maintainer,
    organism,
    common_name,
    provider,
    provider_version,
    release_date,
    release_name,
    BS_genome_object_name,
    source_url,
    organism_biocview,
    src_data_files,
    seqs_srcdir,
    seqfile_name,
    circ_seqs
) {
    paste0(
    "Package: ", package, "\n",
    "Title: ", title, "\n",
    "Description: ", description, "\n",
    "Version: ", version, "\n",
    "Author: ", author, "\n",
    "Maintainer: ", maintainer, "\n",
    "organism: ", organism, "\n",
    "common_name: ", common_name, "\n",
    "provider: ", provider, "\n",
    "provider_version: ", provider_version, "\n",
    "release_date: ", release_date, "\n",
    "release_name: ", release_name, "\n",
    "BSgenomeObjname: ", BS_genome_object_name, "\n",
    "source_url: ", source_url, "\n",
    "organism_biocview: ", organism_biocview, "\n",
    "SrcDataFiles: ", src_data_files, "\n",
    "seqs_srcdir: ", seqs_srcdir, "\n",
    "seqfile_name: ", seqfile_name, "\n",
    "circ_seqs: ", circ_seqs, "\n"
    )
}


convertSeedTxtDcf <- function(seed_name) {
    seed <- read.dcf(
        seed_name,
        fields = NULL,
        all = FALSE,
        keep.white = NULL
    )
    write.dcf(
        seed,
        file = gsub("\\.txt$", ".dcf", seed_name),
        append = FALSE,
        useBytes = FALSE,
        indent = 0.1 * getOption("width"),
        width = 0.9 * getOption("width"),
        keep.white = NULL
    )
}


setwd("/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/2021-1105-1107/")
# getwd()
# list.files()
# list.dirs()

# #  Load gff3 as TxDb, then save as .txdb file
# gff3 <- GenomicFeatures::makeTxDbFromGFF("Mus_musculus_129s1svimj.129S1_SvImJ_v1.104.gff3.gz", format = "auto")
# AnnotationDbi::saveDb(gff3, "129S1-SvImJ.txdb")
# 
# #  Practice loading .db file
# gff3 <- AnnotationDbi::loadDb(file = "129S1-SvImJ.txdb")
# 
# #  Check TxDb object column names
# gff3 %>% AnnotationDbi::columns()
# 
# #  Extract genes
# genes <- GenomicFeatures::genes(gff3, "GENEID")
# genes %>% head()
# 
# #  Extract intergenic regions
# intergenic <- GenomicRanges::gaps(genes)
# intergenic %>% head()

#  The TxDb (taxonomic database) object shows the positions where genes and/or
#+ features are within the genome. However, it does not actually contain DNA
#+ sequence information. To extract DNA sequences, we need to use the BSgenome
#+ package.
# system("mkdir -p 129S1-SvImJ/src_seqdir")
# system("mv 129S1-SvImJ.fa 129S1-SvImJ/src_seqdir")

# (...) #TODO Make sure my names are consistent across local and remote
# (...) #TODO From GB's directories, get GRCm38.SNPs-inserted.CAST-EiJ.129S1-SvImJ
# (   ) #TODO From GB's directories, get GRCm38.SNPs-N-masked.129S1-SvImJ
directories <- c(
    "129S1-SvImJ",
    "CAST-EiJ",
    "GRCm38.SNPs-inserted.129S1-SvImJ",
    "GRCm38.SNPs-inserted.CAST-EiJ",
    "GRCm38.SNPs-inserted.CAST-EiJ.129S1-SvImJ",
    "GRCm38.SNPs-N-masked.129S1-SvImJ",
    "GRCm38.SNPs-N-masked.CAST-EiJ",
    "GRCm38.SNPs-N-masked.CAST-EiJ.129S1-SvImJ"
)

directories_list <- list(
    "129S1-SvImJ" = c(
        "Mmusculus129S1SvImJ",
        "Mus musculus strain 129S1/SvImJ"
    ),
    "CAST-EiJ" = c(
        "MmusculusCASTEiJ",
        "Mus musculus castaneus"
    ),
    "GRCm38.SNPs-inserted.129S1-SvImJ" = c(
        "Mmusculus129Inserted",
        "Mus musculus with 129S1/SvImJ SNPs inserted"
    ),
    "GRCm38.SNPs-inserted.CAST-EiJ" = c(
        "MmusculusCASTInserted",
        "Mus musculus with CAST-EiJ SNPs inserted"
    ),
    "GRCm38.SNPs-inserted.CAST-EiJ.129S1-SvImJ" = c(
        "MmusculusCAST129Inserted",
        "Mus musculus with CAST-EiJ and 129S1/SvImJ SNPs inserted"
    ),
    "GRCm38.SNPs-N-masked.129S1-SvImJ" = c(
        "Mmusculus129Nmasked",
        "Mus musculus with 129S1-SvImJ SNPs N-masked"
    ),
    "GRCm38.SNPs-N-masked.CAST-EiJ" = c(
        "MmusculusCASTNmasked",
        "Mus musculus with CAST-EiJ SNPs N-masked"
    ),
    "GRCm38.SNPs-N-masked.CAST-EiJ.129S1-SvImJ" = c(
        "MmusculusCAST129Nmasked",
        "Mus musculus with CAST-EiJ and 129S1/SvImJ SNPs N-masked"
    )
)

sapply(
    directories,
    function(directories) {
        system(paste0("mkdir -p ", directories, "/src_seqdir"))
    }
)

# sapply(
#     directories,
#     function(directories) {
#         system(paste0("touch ", directories, ".txt "))
#     }
# )

sapply(
    directories,
    function(directories) {
        # system(paste0("mv ", directories, ".txt ", directories, "/src_seqdir"))
        system(paste0("mv ", directories, ".2bit ", directories, "/src_seqdir"))
    }
)


#  129S1-SvImJ ----------------------------------------------------------------
package <- directories_list[["129S1-SvImJ"]][1]
title <- "Full autosome and X sequences for Mus musculus strain 129S1/SvImJ (ENA GCA_001624185)"
description <- "Full autosome and X sequences for Mus musculus strain 129S1/SvImJ as provided by ENA (GCA_001624185; 26 April, 2016) and stored in Biostrings objects."
version <- "0.1"
author <- "Kris Alavattam via Doran et al., Genome Biol 2016, DOI: 10.1186/s13059-016-1024-y"
maintainer <- "Kris Alavattam <kga0@uw.edu>"
organism <- directories_list[["129S1-SvImJ"]][2]
common_name <- "129"
provider <- "ENA"
provider_version <- "0.1"
release_date <- "26 April, 2016"
release_name <- directories_list[["129S1-SvImJ"]][2]
BS_genome_object_name <- "Mus musculus 129S1-SvImJ"
source_url <- "https://www.ebi.ac.uk/ena/browser/view/GCA_001624185.1"
organism_biocview <- "AnnotationData, BSgenome"
src_data_files <- "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_001624185.1?download=true&gzip=true"
seqs_srcdir <- "/Users/kalavattam/Downloads/to-do/2021-1105-1107/129S1-SvImJ/src_seqdir"
seqfile_name <- "129S1-SvImJ.2bit"
circ_seqs <- "character(0)"

seed_name <- paste0("129S1-SvImJ", ".seed.txt")
seed_contents <- writeSeed(
    package = package,
    title = title,
    description = description,
    version = version,
    author = author,
    maintainer = maintainer,
    organism = organism,
    common_name = common_name,
    provider = provider,
    provider_version = provider_version,
    release_date = release_date,
    release_name = release_name,
    BS_genome_object_name = BS_genome_object_name,
    source_url = source_url,
    organism_biocview = organism_biocview,
    src_data_files = src_data_files,
    seqs_srcdir = seqs_srcdir,
    seqfile_name = seqfile_name,
    circ_seqs = circ_seqs
)
readr::write_file(seed_contents, file = seed_name)

#  Create a plain-text seed file
# readr::read_file(seed_name)
readLines(seed_name)

#  Convert .txt file to .dcf format
convertSeedTxtDcf(seed_name = seed_name)

#  Remove output directory if it exists
unlink(
    directories_list[["129S1-SvImJ"]][1],
    recursive = TRUE,
    force = TRUE
)

#  Now, create the BSgenome package
BSgenome::forgeBSgenomeDataPkg(
    gsub("\\.txt$", ".dcf", seed_name)
)

#  Create the package
system(
    paste0("R CMD build ", directories_list[["129S1-SvImJ"]][1])
)
# * checking for file ‘129S1SvImJ/DESCRIPTION’ ...Error : file '/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/2021-1105-1107/129S1SvImJ/DESCRIPTION' is not in valid DCF format
# EXISTS but not correct format

#  Solution: Need to add four spaces to begin line two of field Description

# #  Check if everything is OK with the package using following command in the
# #+ terminal:
# system(
#     paste0("R CMD check ", directories_list[["129S1-SvImJ"]][1], "_0.1.tar.gz")
# )
# # * using log directory ‘/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/2021-1105-1107/Mmusculus129S1SvImJ.Rcheck’
# # * using R version 4.1.0 (2021-05-18)
# # * using platform: x86_64-apple-darwin13.4.0 (64-bit)
# # * using session charset: UTF-8
# # * checking for file ‘Mmusculus129S1SvImJ/DESCRIPTION’ ... OK
# # * this is package ‘Mmusculus129S1SvImJ’ version ‘1.0’
# # * checking package namespace information ... OK
# # * checking package dependencies ... OK
# # * checking if this is a source package ... OK
# # * checking if there is a namespace ... OK
# # * checking for executable files ... OK
# # * checking for hidden files and directories ... OK
# # * checking for portable file names ... OK
# # * checking for sufficient/correct file permissions ... OK
# # * checking whether package ‘Mmusculus129S1SvImJ’ can be installed ... OK
# # * checking installed package size ... NOTE
# # installed size is 652.2Mb
# # sub-directories of 1Mb or more:
# #   extdata  652.1Mb
# # * checking package directory ... OK
# # * checking DESCRIPTION meta-information ... NOTE
# # Package listed in more than one of Depends, Imports, Suggests, Enhances:
# #   ‘BSgenome’
# # A package should be listed in only one of these fields.
# # * checking top-level files ... OK
# # * checking for left-over files ... OK
# # * checking index information ... OK
# # * checking package subdirectories ... OK
# # * checking R files for non-ASCII characters ... OK
# # * checking R files for syntax errors ... OK
# # * checking whether the package can be loaded ... OK
# # * checking whether the package can be loaded with stated dependencies ... OK
# # * checking whether the package can be unloaded cleanly ... OK
# # * checking whether the namespace can be loaded with stated dependencies ... OK
# # * checking whether the namespace can be unloaded cleanly ... OK
# # * checking loading without being on the library search path ... OK
# # * checking dependencies in R code ... OK
# # * checking S3 generic/method consistency ... OK
# # * checking replacement functions ... OK
# # * checking foreign function calls ... OK
# # * checking R code for possible problems ... OK
# # * checking Rd files ... NOTE
# # prepare_Rd: package.Rd:15-17: Dropping empty section \details
# # * checking Rd metadata ... OK
# # * checking Rd cross-references ... OK
# # * checking for missing documentation entries ... OK
# # * checking for code/documentation mismatches ... OK
# # * checking Rd \usage sections ... OK
# # * checking Rd contents ... OK
# # * checking for unstated dependencies in examples ... OK
# # * checking examples ... OK
# # * checking PDF version of manual ... WARNING
# # LaTeX errors when creating PDF version.
# # This typically indicates Rd problems.
# # * checking PDF version of manual without hyperrefs or index ...sh: -c: line 0: syntax error near unexpected token `('
# # sh: -c: line 0: ` '/Users/kalavattam/miniconda3/envs/r41_env/lib/R/bin/R' CMD Rd2pdf  --batch --no-preview --build-dir='/var/folders/bb/bw7y2j396vd5mg4jgddw7ssw0000gn/T//Rtmpjv5Q24/Rd2pdf5958796d92fa' --no-clean --no-index -o  Mmusculus129S1SvImJ-manual.pdf  /Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/2021-1105-1107/Mmusculus129S1SvImJ.Rcheck/Mmusculus129S1SvImJ > 'Rdlatex.log' 2>&1'
# # errorCondition()
# # Re-running with no redirection of stdout/stderr.
# # sh: -c: line 0: syntax error near unexpected token `('
# # sh: -c: line 0: ` '/Users/kalavattam/miniconda3/envs/r41_env/lib/R/bin/R' CMD Rd2pdf  --batch --no-preview --build-dir='/var/folders/bb/bw7y2j396vd5mg4jgddw7ssw0000gn/T//Rtmpjv5Q24/Rd2pdf59588cc995e' --no-clean --no-index -o  Mmusculus129S1SvImJ-manual.pdf  /Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/2021-1105-1107/Mmusculus129S1SvImJ.Rcheck/Mmusculus129S1SvImJ'
# # * DONE
# # 
# # Status: 1 ERROR, 1 WARNING, 3 NOTEs
# # See ‘/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/2021-1105-1107/Mmusculus129S1SvImJ.Rcheck/00check.log’ for details.

#  Install the package using the following command:
system(
    paste0("R CMD install ", directories_list[["129S1-SvImJ"]][1], "_0.1.tar.gz")
)


#  CAST-EiJ -------------------------------------------------------------------
package <- directories_list[["CAST-EiJ"]][1]
title <- "Full autosome and X sequences for Mus musculus castaneus (ENA GCA_001632555)"
description <- "Full autosome and X sequences for Mus musculus castaneus as provided by ENA (GCA_001632555; 26 April, 2016) and stored in Biostrings objects."
version <- "0.1"
author <- "Kris Alavattam via Doran et al., Genome Biol 2016, DOI: 10.1186/s13059-016-1024-y"
maintainer <- "Kris Alavattam <kga0@uw.edu>"
organism <- directories_list[["CAST-EiJ"]][2]
common_name <- "CAST"
provider <- "ENA"
provider_version <- "0.1"
release_date <- "26 April, 2016"
release_name <- directories_list[["CAST-EiJ"]][2]
BS_genome_object_name <- "Mus musculus CAST-EiJ"
source_url <- "https://www.ebi.ac.uk/ena/browser/view/GCA_001632555.1"
organism_biocview <- "AnnotationData, BSgenome"
src_data_files <- "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_001632555.1?download=true&gzip=true"
seqs_srcdir <- "/Users/kalavattam/Downloads/to-do/2021-1105-1107/CAST-EiJ/src_seqdir"
seqfile_name <- "CAST-EiJ.2bit"
circ_seqs <- "character(0)"

seed_name <- paste0("CAST-EiJ", ".seed.txt")
seed_contents <- writeSeed(
    package = package,
    title = title,
    description = description,
    version = version,
    author = author,
    maintainer = maintainer,
    organism = organism,
    common_name = common_name,
    provider = provider,
    provider_version = provider_version,
    release_date = release_date,
    release_name = release_name,
    BS_genome_object_name = BS_genome_object_name,
    source_url = source_url,
    organism_biocview = organism_biocview,
    src_data_files = src_data_files,
    seqs_srcdir = seqs_srcdir,
    seqfile_name = seqfile_name,
    circ_seqs = circ_seqs
)
readr::write_file(seed_contents, file = seed_name)

#  Create a plain-text seed file
# readr::read_file(seed_name)
readLines(seed_name)

#  Convert .txt file to .dcf format
convertSeedTxtDcf(seed_name = seed_name)

#  Remove output directory if it exists
unlink(
    directories_list[["CAST-EiJ"]][1],
    recursive = TRUE,
    force = TRUE
)

#  Now, create the BSgenome package
BSgenome::forgeBSgenomeDataPkg(
    gsub("\\.txt$", ".dcf", seed_name)
)

#  Create the package
system(
    paste0("R CMD build ", directories_list[["CAST-EiJ"]][1])
)

# #  Check if everything is OK with the package using following command in the
# #+ terminal:
# system(
#     paste0("R CMD check ", directories_list[["CAST-EiJ"]][1], "_0.1.tar.gz")
# )

#  Install the package using the following command:
system(
    paste0("R CMD install ", directories_list[["CAST-EiJ"]][1], "_0.1.tar.gz")
)


#  GRCm38.SNPs-inserted.129S1-SvImJ -------------------------------------------
package <- directories_list[["GRCm38.SNPs-inserted.129S1-SvImJ"]][1]
title <- "Full autosome, X, Y, and MT sequences for Mus musculus with 129S1-SvImJ SNPs inserted"
description <- "Full autosome, X, Y, and MT sequences for Mus musculus with 129S1-SvImJ SNPs inserted; stored in Biostrings objects."
version <- "0.1"
author <- "Kris Alavattam"
maintainer <- "Kris Alavattam <kga0@uw.edu>"
organism <- directories_list[["GRCm38.SNPs-inserted.129S1-SvImJ"]][2]
common_name <- "mm10 with 129 SNPs inserted"
provider <- NA
provider_version <- "0.1"
release_date <- NA
release_name <- directories_list[["GRCm38.SNPs-inserted.129S1-SvImJ"]][2]
BS_genome_object_name <- "Mus musculus with 129S1-SvImJ SNPs inserted"
source_url <- NA
organism_biocview <- "AnnotationData, BSgenome"
src_data_files <- NA
seqs_srcdir <- "/Users/kalavattam/Downloads/to-do/2021-1105-1107/GRCm38.SNPs-inserted.129S1-SvImJ/src_seqdir"
seqfile_name <- "GRCm38.SNPs-inserted.129S1-SvImJ.2bit"
circ_seqs <- "character(0)"

seed_name <- paste0("GRCm38.SNPs-inserted.129S1-SvImJ", ".seed.txt")
seed_contents <- writeSeed(
    package = package,
    title = title,
    description = description,
    version = version,
    author = author,
    maintainer = maintainer,
    organism = organism,
    common_name = common_name,
    provider = provider,
    provider_version = provider_version,
    release_date = release_date,
    release_name = release_name,
    BS_genome_object_name = BS_genome_object_name,
    source_url = source_url,
    organism_biocview = organism_biocview,
    src_data_files = src_data_files,
    seqs_srcdir = seqs_srcdir,
    seqfile_name = seqfile_name,
    circ_seqs = circ_seqs
)
readr::write_file(seed_contents, file = seed_name)

#  Create a plain-text seed file
# readr::read_file(seed_name)
readLines(seed_name)

#  Convert .txt file to .dcf format
convertSeedTxtDcf(seed_name = seed_name)

#  Remove output directory if it exists
unlink(
    directories_list[["GRCm38.SNPs-inserted.129S1-SvImJ"]][1],
    recursive = TRUE,
    force = TRUE
)

#  Now, create the BSgenome package
BSgenome::forgeBSgenomeDataPkg(
    gsub("\\.txt$", ".dcf", seed_name)
)

#  Create the package
system(
    paste0("R CMD build ", directories_list[["GRCm38.SNPs-inserted.129S1-SvImJ"]][1])
)

# #  Check if everything is OK with the package using following command in the
# #+ terminal:
# system(
#     paste0("R CMD check ", directories_list[["GRCm38.SNPs-inserted.129S1-SvImJ"]][1], "_0.1.tar.gz")
# )

#  Install the package using the following command:
system(
    paste0("R CMD install ", directories_list[["GRCm38.SNPs-inserted.129S1-SvImJ"]][1], "_0.1.tar.gz")
)


#  GRCm38.SNPs-inserted.CAST-EiJ ----------------------------------------------
package <- directories_list[["GRCm38.SNPs-inserted.CAST-EiJ"]][1]
title <- "Full autosome, X, Y, and MT sequences for Mus musculus with CAST-EiJ SNPs inserted"
description <- "Full autosome, X, Y, and MT sequences for Mus musculus with CAST-EiJ SNPs inserted; stored in Biostrings objects."
version <- "0.1"
author <- "Kris Alavattam"
maintainer <- "Kris Alavattam <kga0@uw.edu>"
organism <- directories_list[["GRCm38.SNPs-inserted.CAST-EiJ"]][2]
common_name <- "mm10 with CAST SNPs inserted"
provider <- NA
provider_version <- "0.1"
release_date <- NA
release_name <- directories_list[["GRCm38.SNPs-inserted.CAST-EiJ"]][2]
BS_genome_object_name <- "Mus musculus with CAST-EiJ SNPs inserted"
source_url <- NA
organism_biocview <- "AnnotationData, BSgenome"
src_data_files <- NA
seqs_srcdir <- "/Users/kalavattam/Downloads/to-do/2021-1105-1107/GRCm38.SNPs-inserted.CAST-EiJ/src_seqdir"
seqfile_name <- "GRCm38.SNPs-inserted.CAST-EiJ.2bit"
circ_seqs <- "character(0)"

seed_name <- paste0("GRCm38.SNPs-inserted.CAST-EiJ", ".seed.txt")
seed_contents <- writeSeed(
    package = package,
    title = title,
    description = description,
    version = version,
    author = author,
    maintainer = maintainer,
    organism = organism,
    common_name = common_name,
    provider = provider,
    provider_version = provider_version,
    release_date = release_date,
    release_name = release_name,
    BS_genome_object_name = BS_genome_object_name,
    source_url = source_url,
    organism_biocview = organism_biocview,
    src_data_files = src_data_files,
    seqs_srcdir = seqs_srcdir,
    seqfile_name = seqfile_name,
    circ_seqs = circ_seqs
)
readr::write_file(seed_contents, file = seed_name)

#  Create a plain-text seed file
# readr::read_file(seed_name)
readLines(seed_name)

#  Convert .txt file to .dcf format
convertSeedTxtDcf(seed_name = seed_name)

#  Remove output directory if it exists
unlink(
    directories_list[["GRCm38.SNPs-inserted.CAST-EiJ"]][1],
    recursive = TRUE,
    force = TRUE
)

#  Now, create the BSgenome package
BSgenome::forgeBSgenomeDataPkg(
    gsub("\\.txt$", ".dcf", seed_name)
)

#  Create the package
system(
    paste0("R CMD build ", directories_list[["GRCm38.SNPs-inserted.CAST-EiJ"]][1])
)

# #  Check if everything is OK with the package using following command in the
# #+ terminal:
# system(
#     paste0("R CMD check ", directories_list[["GRCm38.SNPs-inserted.CAST-EiJ"]][1], "_0.1.tar.gz")
# )

#  Install the package using the following command:
system(
    paste0("R CMD install ", directories_list[["GRCm38.SNPs-inserted.CAST-EiJ"]][1], "_0.1.tar.gz")
)

#  GRCm38.SNPs-inserted.CAST-EiJ.129S1-SvImJ ----------------------------------
package <- directories_list[["GRCm38.SNPs-inserted.CAST-EiJ.129S1-SvImJ"]][1]
title <- "Full autosome, X, Y, and MT sequences for Mus musculus with CAST-EiJ and 129S1-SvImJ SNPs inserted"
description <- "Full autosome, X, Y, and MT sequences for Mus musculus with CAST-EiJ and 129S1-SvImJ SNPs inserted; stored in Biostrings objects."
version <- "0.1"
author <- "Kris Alavattam"
maintainer <- "Kris Alavattam <kga0@uw.edu>"
organism <- directories_list[["GRCm38.SNPs-inserted.CAST-EiJ.129S1-SvImJ"]][2]
common_name <- "mm10 with CAST and 129 SNPs inserted"
provider <- NA
provider_version <- "0.1"
release_date <- NA
release_name <- directories_list[["GRCm38.SNPs-inserted.CAST-EiJ.129S1-SvImJ"]][2]
BS_genome_object_name <- "Mus musculus with CAST-EiJ and 129S1-SvImJ SNPs inserted"
source_url <- NA
organism_biocview <- "AnnotationData, BSgenome"
src_data_files <- NA
seqs_srcdir <- "/Users/kalavattam/Downloads/to-do/2021-1105-1107/GRCm38.SNPs-inserted.CAST-EiJ.129S1-SvImJ/src_seqdir"
seqfile_name <- "GRCm38.SNPs-inserted.CAST-EiJ.129S1-SvImJ.2bit"
circ_seqs <- "character(0)"

seed_name <- paste0("GRCm38.SNPs-inserted.CAST-EiJ.129S1-SvImJ", ".seed.txt")
seed_contents <- writeSeed(
    package = package,
    title = title,
    description = description,
    version = version,
    author = author,
    maintainer = maintainer,
    organism = organism,
    common_name = common_name,
    provider = provider,
    provider_version = provider_version,
    release_date = release_date,
    release_name = release_name,
    BS_genome_object_name = BS_genome_object_name,
    source_url = source_url,
    organism_biocview = organism_biocview,
    src_data_files = src_data_files,
    seqs_srcdir = seqs_srcdir,
    seqfile_name = seqfile_name,
    circ_seqs = circ_seqs
)
readr::write_file(seed_contents, file = seed_name)

#  Create a plain-text seed file
# readr::read_file(seed_name)
readLines(seed_name)

#  Convert .txt file to .dcf format
convertSeedTxtDcf(seed_name = seed_name)

#  Remove output directory if it exists
unlink(
    directories_list[["GRCm38.SNPs-inserted.CAST-EiJ.129S1-SvImJ"]][1],
    recursive = TRUE,
    force = TRUE
)

#  Now, create the BSgenome package
BSgenome::forgeBSgenomeDataPkg(
    gsub("\\.txt$", ".dcf", seed_name)
)

#  Create the package
system(
    paste0("R CMD build ", directories_list[["GRCm38.SNPs-inserted.CAST-EiJ.129S1-SvImJ"]][1])
)

# #  Check if everything is OK with the package using following command in the
# #+ terminal:
# system(
#     paste0("R CMD check ", directories_list[["GRCm38.SNPs-inserted.CAST-EiJ.129S1-SvImJ"]][1], "_0.1.tar.gz")
# )

#  Install the package using the following command:
system(
    paste0("R CMD install ", directories_list[["GRCm38.SNPs-inserted.CAST-EiJ.129S1-SvImJ"]][1], "_0.1.tar.gz")
)


#  GRCm38.SNPs-N-masked.129S1-SvImJ ----------------------------------
package <- directories_list[["GRCm38.SNPs-N-masked.129S1-SvImJ"]][1]
title <- "Full autosome, X, Y, and MT sequences for Mus musculus with 129S1-SvImJ SNPs N-masked"
description <- "Full autosome, X, Y, and MT sequences for Mus musculus with 129S1-SvImJ SNPs N-masked; stored in Biostrings objects."
version <- "0.1"
author <- "Kris Alavattam"
maintainer <- "Kris Alavattam <kga0@uw.edu>"
organism <- directories_list[["GRCm38.SNPs-N-masked.129S1-SvImJ"]][2]
common_name <- "mm10 with 129 SNPs N-masked"
provider <- NA
provider_version <- "0.1"
release_date <- NA
release_name <- directories_list[["GRCm38.SNPs-N-masked.129S1-SvImJ"]][2]
BS_genome_object_name <- "Mus musculus with 129S1-SvImJ SNPs N-masked"
source_url <- NA
organism_biocview <- "AnnotationData, BSgenome"
src_data_files <- NA
seqs_srcdir <- "/Users/kalavattam/Downloads/to-do/2021-1105-1107/GRCm38.SNPs-N-masked.129S1-SvImJ/src_seqdir"
seqfile_name <- "GRCm38.SNPs-N-masked.129S1-SvImJ.2bit"
circ_seqs <- "character(0)"

seed_name <- paste0("GRCm38.SNPs-N-masked.129S1-SvImJ", ".seed.txt")
seed_contents <- writeSeed(
    package = package,
    title = title,
    description = description,
    version = version,
    author = author,
    maintainer = maintainer,
    organism = organism,
    common_name = common_name,
    provider = provider,
    provider_version = provider_version,
    release_date = release_date,
    release_name = release_name,
    BS_genome_object_name = BS_genome_object_name,
    source_url = source_url,
    organism_biocview = organism_biocview,
    src_data_files = src_data_files,
    seqs_srcdir = seqs_srcdir,
    seqfile_name = seqfile_name,
    circ_seqs = circ_seqs
)
readr::write_file(seed_contents, file = seed_name)

#  Create a plain-text seed file
# readr::read_file(seed_name)
readLines(seed_name)

#  Convert .txt file to .dcf format
convertSeedTxtDcf(seed_name = seed_name)

#  Remove output directory if it exists
unlink(
    directories_list[["GRCm38.SNPs-N-masked.129S1-SvImJ"]][1],
    recursive = TRUE,
    force = TRUE
)

#  Now, create the BSgenome package
BSgenome::forgeBSgenomeDataPkg(
    gsub("\\.txt$", ".dcf", seed_name)
)

#  Create the package
system(
    paste0("R CMD build ", directories_list[["GRCm38.SNPs-N-masked.129S1-SvImJ"]][1])
)

# #  Check if everything is OK with the package using following command in the
# #+ terminal:
# system(
#     paste0("R CMD check ", directories_list[["GRCm38.SNPs-N-masked.129S1-SvImJ"]][1], "_0.1.tar.gz")
# )

#  Install the package using the following command:
system(
    paste0("R CMD install ", directories_list[["GRCm38.SNPs-N-masked.129S1-SvImJ"]][1], "_0.1.tar.gz")
)


#  GRCm38.SNPs-N-masked.CAST-EiJ ----------------------------------------------
package <- directories_list[["GRCm38.SNPs-N-masked.CAST-EiJ"]][1]
title <- "Full autosome, X, Y, and MT sequences for Mus musculus N-masked at CAST-EiJ SNPs"
description <- "Full autosome, X, Y, and MT sequences for Mus musculus N-masked at CAST-EiJ SNPs; stored in Biostrings objects."
version <- "0.1"
author <- "Kris Alavattam"
maintainer <- "Kris Alavattam <kga0@uw.edu>"
organism <- directories_list[["GRCm38.SNPs-N-masked.CAST-EiJ"]][2]
common_name <- "mm10 with CAST SNPs N-masked"
provider <- NA
provider_version <- "0.1"
release_date <- NA
release_name <- directories_list[["GRCm38.SNPs-N-masked.CAST-EiJ"]][2]
BS_genome_object_name <- "Mus musculus with CAST-EiJ SNPs N-masked"
source_url <- NA
organism_biocview <- "AnnotationData, BSgenome"
src_data_files <- NA
seqs_srcdir <- "/Users/kalavattam/Downloads/to-do/2021-1105-1107/GRCm38.SNPs-N-masked.CAST-EiJ/src_seqdir"
seqfile_name <- "GRCm38.SNPs-N-masked.CAST-EiJ.2bit"
circ_seqs <- "character(0)"

seed_name <- paste0("GRCm38.SNPs-N-masked.CAST-EiJ", ".seed.txt")
seed_contents <- writeSeed(
    package = package,
    title = title,
    description = description,
    version = version,
    author = author,
    maintainer = maintainer,
    organism = organism,
    common_name = common_name,
    provider = provider,
    provider_version = provider_version,
    release_date = release_date,
    release_name = release_name,
    BS_genome_object_name = BS_genome_object_name,
    source_url = source_url,
    organism_biocview = organism_biocview,
    src_data_files = src_data_files,
    seqs_srcdir = seqs_srcdir,
    seqfile_name = seqfile_name,
    circ_seqs = circ_seqs
)
readr::write_file(seed_contents, file = seed_name)

#  Create a plain-text seed file
# readr::read_file(seed_name)
readLines(seed_name)

#  Convert .txt file to .dcf format
convertSeedTxtDcf(seed_name = seed_name)

#  Remove output directory if it exists
unlink(
    directories_list[["GRCm38.SNPs-N-masked.CAST-EiJ"]][1],
    recursive = TRUE,
    force = TRUE
)

#  Now, create the BSgenome package
BSgenome::forgeBSgenomeDataPkg(
    gsub("\\.txt$", ".dcf", seed_name)
)

#  Create the package
system(
    paste0("R CMD build ", directories_list[["GRCm38.SNPs-N-masked.CAST-EiJ"]][1])
)

# #  Check if everything is OK with the package using following command in the
# #+ terminal:
# system(
#     paste0("R CMD check ", directories_list[["GRCm38.SNPs-N-masked.CAST-EiJ"]][1], "_0.1.tar.gz")
# )

#  Install the package using the following command:
system(
    paste0("R CMD install ", directories_list[["GRCm38.SNPs-N-masked.CAST-EiJ"]][1], "_0.1.tar.gz")
)


#  GRCm38.SNPs-N-masked.CAST-EiJ.129S1-SvImJ ----------------------------------
package <- directories_list[["GRCm38.SNPs-N-masked.CAST-EiJ.129S1-SvImJ"]][1]
title <- "Full autosome, X, Y, and MT sequences for Mus musculus with CAST-EiJ and 129S1-SvImJ SNPs N-masked"
description <- "Full autosome, X, Y, and MT sequences for Mus musculus with CAST-EiJ and 129S1-SvImJ SNPs N-masked; stored in Biostrings objects."
version <- "0.1"
author <- "Kris Alavattam"
maintainer <- "Kris Alavattam <kga0@uw.edu>"
organism <- directories_list[["GRCm38.SNPs-N-masked.CAST-EiJ.129S1-SvImJ"]][2]
common_name <- "mm10 with CAST and 129 SNPs N-masked"
provider <- NA
provider_version <- "0.1"
release_date <- NA
release_name <- directories_list[["GRCm38.SNPs-N-masked.CAST-EiJ.129S1-SvImJ"]][2]
BS_genome_object_name <- "Mus musculus with CAST-EiJ and 129S1-SvImJ SNPs N-masked"
source_url <- NA
organism_biocview <- "AnnotationData, BSgenome"
src_data_files <- NA
seqs_srcdir <- "/Users/kalavattam/Downloads/to-do/2021-1105-1107/GRCm38.SNPs-N-masked.CAST-EiJ.129S1-SvImJ/src_seqdir"
seqfile_name <- "GRCm38.SNPs-N-masked.CAST-EiJ.129S1-SvImJ.2bit"
circ_seqs <- "character(0)"

seed_name <- paste0("GRCm38.SNPs-N-masked.CAST-EiJ.129S1-SvImJ", ".seed.txt")
seed_contents <- writeSeed(
    package = package,
    title = title,
    description = description,
    version = version,
    author = author,
    maintainer = maintainer,
    organism = organism,
    common_name = common_name,
    provider = provider,
    provider_version = provider_version,
    release_date = release_date,
    release_name = release_name,
    BS_genome_object_name = BS_genome_object_name,
    source_url = source_url,
    organism_biocview = organism_biocview,
    src_data_files = src_data_files,
    seqs_srcdir = seqs_srcdir,
    seqfile_name = seqfile_name,
    circ_seqs = circ_seqs
)
readr::write_file(seed_contents, file = seed_name)

#  Create a plain-text seed file
# readr::read_file(seed_name)
readLines(seed_name)

#  Convert .txt file to .dcf format
convertSeedTxtDcf(seed_name = seed_name)

#  Remove output directory if it exists
unlink(
    directories_list[["GRCm38.SNPs-N-masked.CAST-EiJ.129S1-SvImJ"]][1],
    recursive = TRUE,
    force = TRUE
)

#  Now, create the BSgenome package
BSgenome::forgeBSgenomeDataPkg(
    gsub("\\.txt$", ".dcf", seed_name)
)

#  Create the package
system(
    paste0("R CMD build ", directories_list[["GRCm38.SNPs-N-masked.CAST-EiJ.129S1-SvImJ"]][1])
)

# #  Check if everything is OK with the package using following command in the
# #+ terminal:
# system(
#     paste0("R CMD check ", directories_list[["GRCm38.SNPs-N-masked.CAST-EiJ.129S1-SvImJ"]][1], "_0.1.tar.gz")
# )

#  Install the package using the following command:
system(
    paste0("R CMD install ", directories_list[["GRCm38.SNPs-N-masked.CAST-EiJ.129S1-SvImJ"]][1], "_0.1.tar.gz")
)
