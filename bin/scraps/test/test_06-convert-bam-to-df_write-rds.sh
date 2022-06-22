#!/bin/bash

#  test_06-lift-strain-to-mm10.sh
#  KA


#  CAST-EiJ test
Rscript bin/workflow/06-convert-bam-to-df_write-rds.R \
-i ./data/2022-0404_prepare_test-datasets_mm10-CAST/Disteche_sample_6.dedup.mm10.prepro.chr19.repair.bam \
-b ./data/2022-0404_prepare_test-datasets_mm10-CAST/Disteche_sample_6.dedup.mm10.prepro.chr19.bam.bai \
-o ./data/2022-0404_prepare_test-datasets_mm10-CAST \
-s mm10 \
-m TRUE \
-u FALSE \
-a FALSE

#  mm10 test
Rscript bin/workflow/06-convert-bam-to-df_write-rds.R \
-i ./data/2022-0404_prepare_test-datasets_mm10-CAST/Disteche_sample_6.dedup.CAST.prepro.chr19.repair.bam \
-b ./data/2022-0404_prepare_test-datasets_mm10-CAST/Disteche_sample_6.dedup.CAST.prepro.chr19.bam.bai \
-o ./data/2022-0404_prepare_test-datasets_mm10-CAST \
-s CAST \
-m TRUE \
-u FALSE \
-a FALSE

# -i, --bam        bam infile, including path <chr>
# -b, --bai        bam index, including path <chr>
# -s, --strain     strain name to be appended to rds outfile columns;
#                  current options are 'mm10', 'CAST', 'SPRET', and
#                  'CAROLI' <chr>
# -o, --outdir     directory for saving rds outfile, including path
#                  <chr>
# -m, --mated      write rds file for mated reads <logical> [default:
#                  TRUE]
# -u, --unmated    write rds file for unmated reads <logical> [default:
#                  FALSE]
# -a, --ambiguous  write rds file for ambiguous reads <logical>
#                  [default: FALSE]
