#################################################
##################################################
##################################################
#
# On Oct 13, 2018, at 9:01 PM, cdistech@u.washington.edu wrote:
#
# See below.
# C
#
# ---------- Forwarded message ----------
# Date: Sat, 13 Oct 2018 16:03:51 -0700
# From: Vijay Ramani <vij.ramani@gmail.com>
# To: Christine M Disteche <cdistech@u.washington.edu>
# Subject: Re: Hi and single cell data
#
# Hi Christine,sorry for the delay! Doing very well--I just returned from Japan (grant was submitted on time!) and am attempting to stave off the jetlag at the moment. My new
# e-mail is vijay.ramani@ucsf.edu, but my gmail address (vij.ramani@gmail.com) will work as well!
#
# Summarizing the sci-ATAC data we've sequenced already is one of my main priorities for this week; it's located at
# '/net/shendure/vol8/projects/single.cell.combinatorics/nobackup/sciATAC_EB_MKII' but I haven't yet converted it to a format that can be easily analyzed! I will work on having
# some figures / statistics computed by Tuesday's meeting (I will have to Skype in from SF).
# -vijay
#
# On Fri, Oct 12, 2018 at 3:54 PM <cdistech@u.washington.edu> wrote:
#      Hi Vijay,
#
#      How are you? I hope that you are getting settled and ready to submit your grant. Please send your new e-mail.
#
#      I am wondering about the ATAC-seq data on single cells that you generated from the ES cell differentiation time points. We would like to start integrating the data
#      with the Hi-C data. However, it seems that no one knows where your sciATAC-seq data is located. Have you already done analysis of the data? Please let us know where
#      to find the data and what analysis was done.
#
#      I apologize if you have already given this information to Giancarlo.
#
#      Thank you very much,
#
#      Best,
#      Christine
#
##################################################
#
# On Wed, Dec 12, 2018 at 5:02 PM Giancarlo <gbonora@uw.edu> wrote:
# Hi Vijay,
#
# Re. the sciATAC-seq, I’m happy to have a go at analyzing the data.
# However, I’m going to be out the office from next Wed (12/19) until 1/7,
# so will likely only make progress with this in the new year.
# Let me know if you do have a chance to work it before then so that we don’t duplicate work.
# I was going to use Hannah’s pipeline. Do you know if Andrew’s pipeline is very different?
# Or does it simply build on Hannah’s?
#
# Giancarlo
#
##################################################
#
# On Dec 28, 2018, at 3:59 PM, Vijay Ramani <vij.ramani@gmail.com> wrote:
#
# Hi all (cc'ing Bill here as well),
# I've finally been able to take a first pass at the sci-ATAC data using Hannah's pipeline and it seems to make sense.
# Here's a (2D) UMAP embedding of the data (a little over 3,000 single cells).
# This is all without separating by allele, and I still haven't looked at any Cicero contacts yet.
# Plots are faceted in order of differentiation time point (F = female, M = male; experiment was performed to process more female cells than male cells).
#
# <image.png>
#
# If someone could forward me summary slides of the cisTopic + XCI findings from analyzing the analogous sci-Hi-C dataset
# that would be great as we can begin to think about integrating the Hi-C contacts with e.g. Cicero-based pairs of co-accessible sites.
#
# Given the Cell paper from Edith Heard's group earlier this week I think it's probably worth accelerating in terms of analyzing +
# coming up with a story for these data, lest we get scooped.
#
# There is also this preprint from Joe Ecker's lab--has anyone applied this to the new datasets yet?
# https://www.biorxiv.org/content/early/2018/12/27/506717
#
# happy holidays / happy new year to all!
# -v
#
##################################################
#
# On Dec 29, 2018, at 9:33 AM, Vijay Ramani <vij.ramani@gmail.com> wrote:
#
# These topic modeling results look amazing—how has interpreting them been?
# And what are the results from applying this to our mESC and patski data?
#
# The UMAP I’m showing is purely for sci-ATAC data, but I think it nicely illustrates
# the formation of the three germ layers as differentiation happens.
# I need to do some TF activity analyses but am pretty sure we’ll see our germ layer markers
# enriched in each separate cluster.
#
##################################################
#
# On Fri, Dec 28, 2018 at 7:07 PM William Stafford Noble <wnoble@uw.edu> wrote:
#
# Hi all (cc'ing Bill here as well),
# I've finally been able to take a first pass at the sci-ATAC data using Hannah's pipeline and it seems to make sense.
#
# > Can you say a bit more about what "making sense" would mean in this context? I'm not sure how to evaluate whether the plots below are good or bad.
#
# Here's a (2D) UMAP embedding of the data (a little over 3,000 single cells).
#
# > What did you apply UMAP to? Just the vectorized Hi-C data?
#
# This is all without separating by allele, and I still haven't looked at any Cicero contacts yet. Plots are faceted in order of differentiation time point (F = female, M = male; experiment was performed to process more female cells than male cells).
#
# > If someone could forward me summary slides of the cisTopic + XCI findings from analyzing the analogous sci-Hi-C dataset that would be great as we can begin to think about integrating the Hi-C contacts with e.g. Cicero-based pairs of co-accessible sites.
# I am attaching Hyeon-Jin's slides from his research reports (he is cc'd).
# These are an early draft -- Hyeon-Jin, can you send your final version?
# I think Hyeon-Jin is working on writing this up, though the Ecker paper will make that more complicated.
#
# > Given the Cell paper from Edith Heard's group earlier this week I think it's probably worth accelerating in terms of analyzing + coming up with a story for these data, lest we get scooped.
# I agree. I will be out of town this week, but it would be great if we could come up with a concrete plan for how to make sense of this data.
#
# > There is also this preprint from Joe Ecker's lab--has anyone applied this to the new datasets yet?  https://www.biorxiv.org/content/early/2018/12/27/506717
# This paper presents a straightforward and sensible way to visualize scHi-C data.
# The steps are as follows:
# (1) replace each count in the matrix with a sum of the counts in a small square around it,
# (2) row-normalize the matrix and perform random walk smoothing on it,
# (3) binarize the matrix so that it is 80% sparse,
# (4) vectorize each nxn matrix and concatenate matrices from m cells to produce an n^2 x m matrix
# (5) project to lower dimension using PCA.
# Thus, though it's called HiCluster, it's not a clustering method per se but a dimensionality reduction method.
# The results do suggest that it works better than HiCRep+MDS.
# We should definitely try it out!
#
# Bill
#
##################################################
#
# On Mon, Jan 28, 2019 at 3:49 PM Giancarlo <gbonora@uw.edu> wrote:
# I see. I'll pick it up from where Vijay left off then.
#
# Vijay, presumably that’s the case and it won’t be necessary to remap/reprocess the read data and that I can work from the data files in this folder:
# /net/shendure/vol8/projects/single.cell.combinatorics/nobackup/sciATAC_EB_MKII/sci-ATAC_out
# Is that correct?
#
# > On Dec 28, 2018, at 3:59 PM, Vijay Ramani <vij.ramani@gmail.com> wrote:
# > I've finally been able to take a first pass at the sci-ATAC data using Hannah's pipeline and it seems to make sense. Here's a (2D) UMAP embedding of the data…
# Would you mind sharing the code that you used to perform that analysis so that I can see how the data was processed and which file(s) are the ones that I should work from.
#
# Also, could you point me to the pipeline that you used to process the data. Last year Hannah granted me access to the following link https://github.com/cole-trapnell-lab/ATAC-Pipeline-Dev, but it no longer appears to work.
#
# Thanks,
# Giancarlo
#
##################################################
#
# On Jan 28, 2019, at 4:09 PM, Vijay Ramani <vij.ramani@gmail.com> wrote:
#
# Okay going backwards:
#
# > Also, could you point me to the pipeline that you used to process the data.
# > Last year Hannah granted me access to the following link https://github.com/cole-trapnell-lab/ATAC-Pipeline-Dev,
# > but it no longer appears to work.
# Link still works for me...I would make sure that the GitHub account you're using is the same one Hannah granted permissions for.
#
# > Would you mind sharing the code that you used to perform that analysis so that I can see how the data was processed and which file(s) are the ones that I should work from.
# Happy to share. I'll comment and try to zip everything and try to send by e.o.d. tomorrow.
#
# > Vijay, presumably that’s the case and it won’t be necessary to remap/reprocess the read data and that I can work from the data files in this folder:
# > Is that correct?
# Yes.
#
# From Bill:
# > I don’t think Vijay has time to do much analysis anymore. His lab is up and running.
# I'm very much interested in seeing this work through to submission given the amount of effort
# I've already put into the project, but realistically only have the bandwidth at this point to
# suggest experiments & visualizations to try, not code and explore the data myself.
#
# The topic analysis results look very interesting;
# I think it would be great to think of ways to view these data in light of our hypotheses about XCI compaction (multi-state vs. biphasic).
# I'll be able to call in for the first 30 min. of Wednesday's meeting.
#
# -v
#
##################################################

source ~/.bashrc
source ~/.bash_profile

PROJDIR="2019_sciATAC_analysis"
# mkdir -p "${PROJDIR}"
WORKDIR="20190208_sciATAC_mouseDiff"
READDIR="data_20190208_sciATAC_mouseDiff"


##################################################
# Working folder
mkdir -p "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"
cd  "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"










#################################################
#################################################
#################################################
# Set up data folder

VOL8HOME=/net/noble/vol8/gbonora/sciOmics
mkdir -p "${VOL8HOME}"
cd "${VOL8HOME}"
pwd -P

# Create a folder for data
mkdir -p "${VOL8HOME}"/"${PROJDIR}"/data/"${READDIR}"

# Link to PROJDIR
cd "${HOME}"/proj/"${PROJDIR}"
ln -sf "${VOL8HOME}"/"${PROJDIR}"/data .
ls -ltrh data/










#################################################
#################################################
#################################################
# Step 0. Link to mouse raw data

# # Make data folder
# ls -ltrh "${HOME}"/proj/"${PROJDIR}"/data
# mkdir -p "${HOME}"/proj/"${PROJDIR}"/data/"${READDIR}"/rawReads
# cd "${HOME}"/proj/"${PROJDIR}"/data/"${READDIR}"/rawReads

# url=http://krishna.gs.washington.edu/content/members/cholilee/sequencing_files/180512_NS500488_0616_AHLYGNBGX5_demultiplexed/
# suffix="_001.fastq.gz"
# for sample in Undetermined_S0; do
#     echo ${sample}
#     for read in I1 R1 R2; do
#       job=wget-${sample}_${read}.job
#       echo wget $url/${sample}_${read}${suffix} > $job
#       qsub -j y -l h_rt=7:59:59 -cwd $job
#     done
# done
# qstat -u gbonora

# SRCDIR=/net/shendure/vol9/seq/NEXTSEQ/180507_NS500488_0611_AH2KHCAFXY
SRCDIR=/net/shendure/vol9/seq/NEXTSEQ/180822_NS500488_0676_AH5FCGAFXY









#################################################
#################################################
#################################################
# Get Hannah pipeline setup
#
# See '20190208_cole-trapnell-lab_ATAC-Pipeline-Dev_githubRepo.sh'










#################################################
#################################################
#################################################
# Vijay work:

ln -sf /net/shendure/vol8/projects/single.cell.combinatorics/nobackup/sciATAC_EB_MKII/ vijayWork

ls -ltrh /net/shendure/vol8/projects/single.cell.combinatorics/nobackup/sciATAC_EB_MKII/*

# Download
WDIR="2019_sciATAC_analysis/results/gbonora/20190208_sciATAC_mouseDiff/vijayWork/"
rsync -auvPKLn --exclude="*.counts.txt" --include "*.R" --include="*.py" --include="*.txt" --include="*.by.sample" --include="*.pdf" --include="*.csv" --include='*/knee-plots/*' --include='*/final-output/*' --prune-empty-dirs --include='*/' --exclude='*' gs:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/
rsync -auvPKLn --include="*.counts.txt.chr12" --prune-empty-dirs --include='*/' --exclude='*' gs:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/










#################################################
#################################################
#################################################
# Set up data description files

cd "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"

# Setup DATADIR
cd "${HOME}"/proj/"${PROJDIR}"/data/"${READDIR}"
pwd -P > "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/DATADIR
cat "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/DATADIR
# /net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20190208_sciATAC_mouseDiff

cd "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"
ls -1 "${HOME}"/proj/"${PROJDIR}"/results/gbonora/"${WORKDIR}"/vijayWork/sciATAC_EB_MKII.split.q10.sort.bam
echo "sciATAC_EB_MKII" > libIDs
cat libIDs
# sciATAC_EB_MKII










#################################################
#################################################
#################################################
# Alellic segregation
# see '20190208_sciATAC_mouseDiff_allelicSegregation_workflow.sh'











##################################################
##################################################
##################################################
# Analyze data
# 20190208
#
# See 'analyzeProcessedData'

#################################################
# Setup environment

module load git/2.19.1
which git
# /net/gs/vol3/software/modules-sw/git/2.19.1/Linux/RHEL6/x86_64/bin/git

module unload R
module load R/3.5.0

# module unload python
# module load python/3.3.3




















##################################################
##################################################
##################################################
# Get results

WDIR="2019_sciATAC_analysis/results/gbonora/20190208_sciATAC_mouseDiff"

# HTML and txt
rsync -auvPKLn --exclude='allelicWork*' --include='*.tsv' --include='*.txt' --include='*.html' --include='*/' --prune-empty-dirs  --exclude='*' gs:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/

# Misc
rsync -auvPKLn --include='20180516_sciATAC_WTandDel1Patski_pilot_segregationResults.txt' --include='*/' --prune-empty-dirs  --exclude='*' gs:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/

# Images
rsync -auvPKLn --include='*.png' --include='*.pdf' --include='*/' --prune-empty-dirs  --exclude='*' gs:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/

# Upload images
rsync -auvPKLn --exclude '*/analysisOutput.cisTopic.testing*/*' --include='*.png' --include='*.pdf' --include='*/' --prune-empty-dirs  --exclude='*' ~/NobleLab/proj/"${WDIR}"/ gs:proj/"${WDIR}"/

# upload code, inc. R analysis
rsync -avuPKLn --include "*.txt" --include "*.py" --include "*.R" --include "*.sh" --include='*/' --prune-empty-dirs  --exclude='*' ~/NobleLab/proj/"${WDIR}"/ gs:proj/"${WDIR}"/
