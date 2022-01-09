##################################################################################################
##################################################################################################
##################################################################################################
# 20191108 
# N-masked Cast/129 assemblies


#################################################
#################################################
#################################################
# 20191108
# N-masked mm10 assemblies using Cast and 129 SNPs (not validated)

mkdir ~/refdata/mm10.Nmasked.Cast129
cd ~/refdata/mm10.Nmasked.Cast129

~/bin/SNPsplit_v0.3.2/SNPsplit_genome_preparation --help | more
# Dual strain mode:
#    1) The VCF file is read and filtered for high-confidence SNPs in the strain specified with --strain <name>
#    2) The reference genome (given with --reference_genome <genome>) is read into memory, and the filtered high-
#       confidence SNP positions are incorporated as full sequence and optionally as N-masking
#    3) The VCF file is read one more time and filtered for high-confidence SNPs in strain 2 specified with --strain2 <name>
#    4) The filtered high-confidence SNP positions of strain 2 are incorporated as full sequence and optionally as N-masking
#    5) The SNP information of strain and strain 2 relative to the reference genome build are compared, and a new Ref/SNP
#       annotation is constructed whereby the new Ref/SNP information will be Strain/Strain2 (and no longer the standard
#       reference genome strain Black6 (C57BL/6J))
#    6) The full genome sequence given with --strain <name> is read into memory, and the high-confidence SNP positions between
#       Strain and Strain2 are incorporated as full sequence and optionally as N-masking
#
# ...
# 
# --dual_hybrid                 Optional. The resulting genome will no longer relate to the original reference specified with '--reference_genome'.
#                               Instead the new Reference (Ref) is defined by '--strain <strain_name>' and the new SNP genome
#                               is defined by '--strain2 <strain_name>'. '--dual_hybrid' automatically sets '--full_sequence'.
# 
#                               This will invoke a multi-step process:
#                                  1) Read/filter SNPs for first strain (specified with '--strain <strain_name>')
#                                  2) Write full SNP incorporated (and optionally N-masked) genome sequence for first strain
#                                  3) Read/filter SNPs for second strain (specified with '--strain2 <strain_name>')
#                                  4) Write full SNP incorporated (and optionally N-masked) genome sequence for second strain
#                                  5) Generate new Ref/Alt SNP annotations for Strain1/Strain2
#                                  6) Set first strain as new reference genome and construct full SNP incorporated (and optionally
#                                     N-masked) genome sequences for Strain1/Strain2

zcat $HOME/refdata/mm10pseudoSpretus/mgp.v5.merged.snps_all.dbSNP142.vcf.gz | more
# Note: uses Ensemble notation i.e. no 'chr' in front of chromosome #.
# 129S1_SvImJ
# CAST_EiJ

# So use Ensemble assembly to generate N-masked assembly then add 'chr' to chrom names.
# See above for explanation.
refFile=/net/noble/vol2/home/gbonora/refdata/iGenomes/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta
altSnpFile=$HOME/refdata/mm10pseudoSpretus/mgp.v5.merged.snps_all.dbSNP142.vcf.gz


jobfile=SNPsplit_genome_preparation.Cast129.job
echo ${jobfile}

cat << EOF > "${jobfile}"
#!/bin/bash -x
#\$ -l mfree=64G
#\$ -cwd
#\$ -j y

source $HOME/.bashrc
source $HOME/.bash_profile

hostname
printf "\n\nstart: %s\n\n" "\$( date )"

# NOTE: --dual_hybrid mode 
#       and no --skip_filtering
~/bin/SNPsplit_v0.3.2/SNPsplit_genome_preparation \
--dual_hybrid \
--strain CAST_EiJ \
--strain2 129S1_SvImJ \
--nmasking \
--vcf_file ${altSnpFile} \
--reference_genome ${refFile} \
--genome_build mm10.Nmasked.Cast129 

printf "\n\nend: %s\n\n" "\$( date )"
EOF

chmod u+x *.job

qsub ${jobfile}

ls -ltrh ~/refdata/mm10.Nmasked.Cast129
ls -ltrh ls -ltrh ~/refdata/mm10.Nmasked.Cast129/CAST_EiJ_129S1_SvImJ_dual_hybrid.based_on_mm10.Nmasked.Cast129_N-masked
# total 2.6G
# -rw-rw-r-- 1 gbonora noblelab 118M Nov  8 18:23 chr11.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 141M Nov  8 18:23 chr7.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  89M Nov  8 18:23 chrY.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 176M Nov  8 18:23 chr2.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  92M Nov  8 18:23 chr17.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 189M Nov  8 18:23 chr1.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  88M Nov  8 18:23 chr18.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 116M Nov  8 18:23 chr13.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  95M Nov  8 18:23 chr16.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 145M Nov  8 18:23 chr6.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 155M Nov  8 18:24 chr3.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 165M Nov  8 18:24 chrX.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  17K Nov  8 18:24 chrMT.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 121M Nov  8 18:24 chr9.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 116M Nov  8 18:24 chr12.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 121M Nov  8 18:24 chr14.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 101M Nov  8 18:24 chr15.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 125M Nov  8 18:24 chr8.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 151M Nov  8 18:24 chr4.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 126M Nov  8 18:24 chr10.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab  60M Nov  8 18:24 chr19.N-masked.fa
# -rw-rw-r-- 1 gbonora noblelab 147M Nov  8 18:24 chr5.N-masked.fa