# Example code to run demux.

# 1.bcl2fastq
# 2.barcodecorrect and split by samples
# 3.trimmomatic

### demux 
./demux_sciatac \ 
    --rundir /net/shendure/vol9/seq/NOVASEQ/200827_Riza_CX_NWGC/net/grc/vol6/data/GRC45NS/raw/200826_GRC45NS_0256_A_HKGH5DSXY \ 
    --outdir /net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/demux_run3 \ 
    --samplesheet /net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/sample_sheet.txt \ 
    --chemistry 3level_384
    
