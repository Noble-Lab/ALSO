#!/bin/bash

#  4dn_calls.sh

# -----------------------------------------------------------------------------
#  Called 2021-1003, run "C57BL-6NJ", run "CAST-EiJ"
./analyze_sciatac \
--demux_output /net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run1 /net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run2 /net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run3 \
--outdir /net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1003_Disteche_C57BL-6NJ \
--genome C57BL-6NJ \
--samples Disteche_sample_1 Disteche_sample_2 Disteche_sample_3 Disteche_sample_4 Disteche_sample_5 Disteche_sample_6 Disteche_sample_7 Disteche_sample_8 Disteche_sample_9 Disteche_sample_10 Disteche_sample_11 Disteche_sample_12 Disteche_sample_13 Disteche_sample_14 Disteche_sample_15 Disteche_sample_16 Disteche_sample_17 Disteche_sample_18 Disteche_sample_19 Disteche_sample_20 Disteche_sample_21 Disteche_sample_22

./analyze_sciatac \
--demux_output /net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run1 /net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run2 /net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run3 \
--outdir /net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1003_Disteche_CAST-EiJ \
--genome CAST-EiJ \
--samples Disteche_sample_1 Disteche_sample_2 Disteche_sample_3 Disteche_sample_4 Disteche_sample_5 Disteche_sample_6 Disteche_sample_7 Disteche_sample_8 Disteche_sample_9 Disteche_sample_10 Disteche_sample_11 Disteche_sample_12 Disteche_sample_13 Disteche_sample_14 Disteche_sample_15 Disteche_sample_16 Disteche_sample_17 Disteche_sample_18 Disteche_sample_19 Disteche_sample_20 Disteche_sample_21 Disteche_sample_22

#  Experiment failure b/c I mapped using indices built against GCA .fa files.
#+ This is a problem b/c these contain duplicate chromosomes that, in turn,
#+ throw errors in various steps of the pipeline (e.g., with samtools); see...
#+ (reference job out files)


# -----------------------------------------------------------------------------
#  Called 2021-1007, run "mm10_Cast_129_Nmasked"
./analyze_sciatac \
--demux_output "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1007_Bonora-et-al_mm10-CAST-129-Nmasked" \
--genome "mm10_Cast_129_Nmasked" \
--samples \
"F123-CASTx129.embryoid-body-day-3" \
"F123-CASTx129.undifferentiated" \
"F121-6-CASTx129.embryoid-body-day-3" \
"F121-6-CASTx129.embryoid-body-day-11" \
"F123-CASTx129.embryoid-body-day-11" \
"F121-6-CASTx129.embryoid-body-day-7" \
"F121-6-CASTx129.undifferentiated" \
"F121-6-CASTx129.npcs" \
"F123-CASTx129.embryoid-body-day-7"

#  Strange errors encountered; see strange-errors.txt for details


# -----------------------------------------------------------------------------
#  Called 2021-1008
./analyze_sciatac \
--demux_output "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1007_Bonora-et-al_mm10-CAST-129-Nmasked" \
--genome "mm10_Cast_129_Nmasked" \
--samples \
"F123-CASTx129.embryoid-body-day-3" \
"F123-CASTx129.undifferentiated" \
"F121-6-CASTx129.embryoid-body-day-3" \
"F121-6-CASTx129.embryoid-body-day-11" \
"F123-CASTx129.embryoid-body-day-11" \
"F121-6-CASTx129.embryoid-body-day-7" \
"F121-6-CASTx129.undifferentiated" \
"F121-6-CASTx129.npcs" \
"F123-CASTx129.embryoid-body-day-7" \
--dry

#  Manually add what's missing to see if the pipeline continues on its own
touch results/2021-1007_Bonora-et-al_mm10-CAST-129-Nmasked/align_reads/.run.0.lane1.F121-6-CASTx129.embryoid-body-day-7.bam.done
touch results/2021-1007_Bonora-et-al_mm10-CAST-129-Nmasked/align_reads/.run.0.lane1.F123-CASTx129.undifferentiated.bam.done


#  Re-running "mm10_Cast_129_Nmasked"; yes, the pipeline does recognize what has
#+ been done and continues to run forward
./analyze_sciatac \
--demux_output "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1007_Bonora-et-al_mm10-CAST-129-Nmasked" \
--genome "mm10_Cast_129_Nmasked" \
--samples \
"F123-CASTx129.embryoid-body-day-3" \
"F123-CASTx129.undifferentiated" \
"F121-6-CASTx129.embryoid-body-day-3" \
"F121-6-CASTx129.embryoid-body-day-11" \
"F123-CASTx129.embryoid-body-day-11" \
"F121-6-CASTx129.embryoid-body-day-7" \
"F121-6-CASTx129.undifferentiated" \
"F121-6-CASTx129.npcs" \
"F123-CASTx129.embryoid-body-day-7"


#  Run CAST-EiJ using appropriate .fa, Bowtie2 indices
./analyze_sciatac \
--demux_output "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1008_Bonora-et-al_CAST-EiJ" \
--genome "CAST-EiJ" \
--samples \
"F123-CASTx129.embryoid-body-day-3" \
"F123-CASTx129.undifferentiated" \
"F121-6-CASTx129.embryoid-body-day-3" \
"F121-6-CASTx129.embryoid-body-day-11" \
"F123-CASTx129.embryoid-body-day-11" \
"F121-6-CASTx129.embryoid-body-day-7" \
"F121-6-CASTx129.undifferentiated" \
"F121-6-CASTx129.npcs" \
"F123-CASTx129.embryoid-body-day-7"


# -----------------------------------------------------------------------------
#  Called 2021-1023-1024, mm10_Cast_129_Nmasked

# conda activate pipeline-test_env  # Always run w/this environment activated

# conda install -c conda-forge \
# r-readr r-stringr r-irlba \
# r-varhandle r-rtsne r-plyr \
# r-ggplot2 r-patchwork r-pdist \
# r-heatmaply r-rcolorbrewer r-mass \
# r-devtools

# conda install -c bioconda \
# r-seurat bioconductor-rgreat bioconductor-rcistarget \
# bioconductor-monocle r-seqminer bioconductor-ensdb.hsapiens.v75 \
# bioconductor-ggbio bioconductor-biomart r-gplots

# $ R
# R> devtools::install_github("aertslab/RcisTarget")
# R> devtools::install_github("aertslab/AUCell")
# R> devtools::install_github("aertslab/cisTopic")

cd /net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/src/sciatac_pipeline || exit 1

#  Re-running "mm10_Cast_129_Nmasked"
./analyze_sciatac \
--demux_output "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1007_Bonora-et-al_mm10-CAST-129-Nmasked" \
--genome "mm10_Cast_129_Nmasked" \
--samples \
"F123-CASTx129.embryoid-body-day-3" \
"F123-CASTx129.undifferentiated" \
"F121-6-CASTx129.embryoid-body-day-3" \
"F121-6-CASTx129.embryoid-body-day-11" \
"F123-CASTx129.embryoid-body-day-11" \
"F121-6-CASTx129.embryoid-body-day-7" \
"F121-6-CASTx129.undifferentiated" \
"F121-6-CASTx129.npcs" \
"F123-CASTx129.embryoid-body-day-7"


#  Re-running "mm10_Cast_129_Nmasked" from scratch
#  Unsure how to move past the errors reported under the header labeled "failed
#+ job, job 4" in "strange-errors.txt" (#TODO Rename this file and/or copy
#+ it into your lab notebook), therefore re-running everything from scratch

#  First, turn on tmux
cd "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/src/sciatac_pipeline" || exit 1

mkdir -p "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1023_Bonora-et-al_mm10-CAST-129-Nmasked"

./analyze_sciatac \
--demux_output "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1023_Bonora-et-al_mm10-CAST-129-Nmasked" \
--genome "mm10_Cast_129_Nmasked" \
--samples \
"F123-CASTx129.embryoid-body-day-3" \
"F123-CASTx129.undifferentiated" \
"F121-6-CASTx129.embryoid-body-day-3" \
"F121-6-CASTx129.embryoid-body-day-11" \
"F123-CASTx129.embryoid-body-day-11" \
"F121-6-CASTx129.embryoid-body-day-7" \
"F121-6-CASTx129.undifferentiated" \
"F121-6-CASTx129.npcs" \
"F123-CASTx129.embryoid-body-day-7"

#  Jobs failed because conda environment was not inherited from launch
#+ environment; updated analyze_sciatac.py to activate the conda environment
#+ where necessary; renamed the _log directory so that .o* and .e* will not be
#+ cleared; re-running as yesterday. Let's see how for we can get.
./analyze_sciatac \
--demux_output "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1023_Bonora-et-al_mm10-CAST-129-Nmasked" \
--genome "mm10_Cast_129_Nmasked" \
--samples \
"F123-CASTx129.embryoid-body-day-3" \
"F123-CASTx129.undifferentiated" \
"F121-6-CASTx129.embryoid-body-day-3" \
"F121-6-CASTx129.embryoid-body-day-11" \
"F123-CASTx129.embryoid-body-day-11" \
"F121-6-CASTx129.embryoid-body-day-7" \
"F121-6-CASTx129.undifferentiated" \
"F121-6-CASTx129.npcs" \
"F123-CASTx129.embryoid-body-day-7"

#> CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'.
#> To initialize your shell, run
#>
#>     $ conda init <SHELL_NAME>
#>
#> Currently supported shells are:
#>   - bash
#>   - fish
#>   - tcsh
#>   - xonsh
#>   - zsh
#>   - powershell
#>
#> See 'conda init --help' for more information and options.
#>
#> IMPORTANT: You may need to close and restart your shell after running 'conda init'.

#  Getting rid of the conda calls, unhashing the 'module load' calls, and then
#+ re-running
./analyze_sciatac \
--demux_output "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1023_Bonora-et-al_mm10-CAST-129-Nmasked" \
--genome "mm10_Cast_129_Nmasked" \
--samples \
"F123-CASTx129.embryoid-body-day-3" \
"F123-CASTx129.undifferentiated" \
"F121-6-CASTx129.embryoid-body-day-3" \
"F121-6-CASTx129.embryoid-body-day-11" \
"F123-CASTx129.embryoid-body-day-11" \
"F121-6-CASTx129.embryoid-body-day-7" \
"F121-6-CASTx129.undifferentiated" \
"F121-6-CASTx129.npcs" \
"F123-CASTx129.embryoid-body-day-7"


# -----------------------------------------------------------------------------
#  2021-1024, CAST-EiJ
./analyze_sciatac \
--demux_output "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1024_Bonora-et-al_CAST-EiJ" \
--genome "CAST-EiJ" \
--samples \
"F123-CASTx129.embryoid-body-day-3" \
"F123-CASTx129.undifferentiated" \
"F121-6-CASTx129.embryoid-body-day-3" \
"F121-6-CASTx129.embryoid-body-day-11" \
"F123-CASTx129.embryoid-body-day-11" \
"F121-6-CASTx129.embryoid-body-day-7" \
"F121-6-CASTx129.undifferentiated" \
"F121-6-CASTx129.npcs" \
"F123-CASTx129.embryoid-body-day-7"


# -----------------------------------------------------------------------------
#  2021-1024, 129S1-SvImJ
./analyze_sciatac \
--demux_output "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1024_Bonora-et-al_129S1-SvImJ" \
--genome "129S1-SvImJ" \
--samples \
"F123-CASTx129.embryoid-body-day-3" \
"F123-CASTx129.undifferentiated" \
"F121-6-CASTx129.embryoid-body-day-3" \
"F121-6-CASTx129.embryoid-body-day-11" \
"F123-CASTx129.embryoid-body-day-11" \
"F121-6-CASTx129.embryoid-body-day-7" \
"F121-6-CASTx129.undifferentiated" \
"F121-6-CASTx129.npcs" \
"F123-CASTx129.embryoid-body-day-7"

#  Both the CAST-EiJ and 129S1-SvImJ calls result in empty .bam files: What's
#+ going on? I suspect it might have to do with the bowtie2 indices, so I am
#+ re-generating them now


# -----------------------------------------------------------------------------
#  2021-1024, 129S1-SvImJ
./analyze_sciatac \
--demux_output "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1024_Bonora-et-al_129S1-SvImJ" \
--genome "129S1-SvImJ" \
--samples \
"F123-CASTx129.embryoid-body-day-3" \
"F123-CASTx129.undifferentiated" \
"F121-6-CASTx129.embryoid-body-day-3" \
"F121-6-CASTx129.embryoid-body-day-11" \
"F123-CASTx129.embryoid-body-day-11" \
"F121-6-CASTx129.embryoid-body-day-7" \
"F121-6-CASTx129.undifferentiated" \
"F121-6-CASTx129.npcs" \
"F123-CASTx129.embryoid-body-day-7"

#  What happens with another reference such as 'mm10'?


# -----------------------------------------------------------------------------
#  2021-1024, mm10
./analyze_sciatac \
--demux_output "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1024_Bonora-et-al_mm10" \
--genome "mm10" \
--samples \
"F123-CASTx129.embryoid-body-day-3" \
"F123-CASTx129.undifferentiated" \
"F121-6-CASTx129.embryoid-body-day-3" \
"F121-6-CASTx129.embryoid-body-day-11" \
"F123-CASTx129.embryoid-body-day-11" \
"F121-6-CASTx129.embryoid-body-day-7" \
"F121-6-CASTx129.undifferentiated" \
"F121-6-CASTx129.npcs" \
"F123-CASTx129.embryoid-body-day-7"


# -----------------------------------------------------------------------------
#  2021-1025, 129S1-SvImJ
./analyze_sciatac \
--demux_output "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1025_Bonora-et-al_129S1-SvImJ" \
--genome "129S1-SvImJ" \
--samples \
"F123-CASTx129.embryoid-body-day-3" \
"F123-CASTx129.undifferentiated" \
"F121-6-CASTx129.embryoid-body-day-3" \
"F121-6-CASTx129.embryoid-body-day-11" \
"F123-CASTx129.embryoid-body-day-11" \
"F121-6-CASTx129.embryoid-body-day-7" \
"F121-6-CASTx129.undifferentiated" \
"F121-6-CASTx129.npcs" \
"F123-CASTx129.embryoid-body-day-7"

#  Even after regenerating Bowtie2 indices, still getting empty .bam files...


# -----------------------------------------------------------------------------
#  2021-1025, 129S1-SvImJ (after renaming fasta headers and regenerating Bowtie2 indices)
./analyze_sciatac \
--demux_output "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1026_Bonora-et-al_129S1-SvImJ_test" \
--genome "129S1-SvImJ" \
--samples "F121-6-CASTx129.undifferentiated"


# -----------------------------------------------------------------------------
#  2021-1027, 129S1-SvImJ (after renaming fasta headers and regenerating Bowtie2 indices)
./analyze_sciatac \
--demux_output "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1027_Bonora-et-al_CAST-EiJ_test" \
--genome "CAST-EiJ" \
--samples "F121-6-CASTx129.undifferentiated"

#  Just want to see what these .bams look like...
#  ( Y ) CAST-EiJ
./analyze_sciatac \
--demux_output \
"/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run1" \
"/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run2" \
"/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run3" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1027_Disteche_CAST-EiJ_test" \
--genome "CAST-EiJ" \
--samples "Disteche_sample_1"

# -------------------------------------
#  ( Y ) mm10
./analyze_sciatac \
--demux_output \
"/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run1" \
"/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run2" \
"/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run3" \
--outdir "/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/results/2021-1028_Disteche_mm10_test" \
--genome "mm10" \
--samples "Disteche_sample_1"

# -------------------------------------
#  ( Y ) mm10_Cast_Nmasked
./analyze_sciatac \
--demux_output \
"/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run1" \
"/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run2" \
"/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run3" \
--outdir "/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results/2021-1031_Disteche_mm10_Cast_Nmasked_test" \
--genome "mm10_Cast_Nmasked" \
--samples "Disteche_sample_1"

# -------------------------------------
#  (...) C57BL-6NJ
./analyze_sciatac \
--demux_output \
"/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run1" \
"/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run2" \
"/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run3" \
--outdir "/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results/2021-1101_Disteche_C57BL-6NJ_test" \
--genome "C57BL-6NJ" \
--samples "Disteche_sample_1"


# -----------------------------------------------------------------------------
#  Hacking to get the pipeline to accept the .fastq files, 2021-1007
path="/net/noble/vol7/kga0/2021_kga0_4dn-mouse-cross/data/Bonora-et-al_Genome-Res_2021"
cd "${path}" || 
    {
        echo "Failure to cd to ${path}. Exiting."
        exit 1
    }

typeset -A array_4DN=(
    ["4DNFIVZJ5RRW"]="R2_001.F123-CASTx129.embryoid-body-day-3"
    ["4DNFIUILMYPA"]="R1_001.F123-CASTx129.embryoid-body-day-3"
    ["4DNFICX2I184"]="R2_001.F123-CASTx129.undifferentiated"
    ["4DNFI3JJJT42"]="R1_001.F123-CASTx129.undifferentiated"
    ["4DNFIKZTPHR3"]="R1_001.F121-6-CASTx129.embryoid-body-day-3"
    ["4DNFI3XCABWA"]="R2_001.F121-6-CASTx129.embryoid-body-day-3"
    ["4DNFIVPR1OSK"]="R2_001.F121-6-CASTx129.embryoid-body-day-11"
    ["4DNFIRZCIP6H"]="R1_001.F121-6-CASTx129.embryoid-body-day-11"
    ["4DNFIW8T1BRD"]="R2_001.F123-CASTx129.embryoid-body-day-11"
    ["4DNFIU2DDZCT"]="R1_001.F123-CASTx129.embryoid-body-day-11"
    ["4DNFIBUMQ47F"]="R1_001.F121-6-CASTx129.embryoid-body-day-7"
    ["4DNFIASEMP76"]="R2_001.F121-6-CASTx129.embryoid-body-day-7"
    ["4DNFIKQMIC2N"]="R2_001.F121-6-CASTx129.undifferentiated"
    ["4DNFIHPWGPDU"]="R1_001.F121-6-CASTx129.undifferentiated"
    ["4DNFIS87JUVK"]="R1_001.F121-6-CASTx129.npcs"
    ["4DNFI87X4FRN"]="R2_001.F121-6-CASTx129.npcs"
    ["4DNFIDWXGMHE"]="R1_001.F123-CASTx129.embryoid-body-day-7"
    ["4DNFIDS8GTWG"]="R2_001.F123-CASTx129.embryoid-body-day-7"
)

mkdir -p final_outs
for i in "${!array_4DN[@]}"; do
    echo "  key: ${i}"
    echo "value: ${array_4DN[${i}]}"
    echo "start copy"
    cp "./${i}.fastq.gz" "./final_outs/Undetermined_S0_L001_${array_4DN[${i}]}.trimmed.fastq.gz"
    echo "copy completed"
    echo ""
done
