#!/bin/bash

cd /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross

./bin/sciatac_pipeline/analyze_sciatac \
--demux_output data/demux_run1 data/demux_run2 data/demux_run3 \
--outdir results/2022-02-23/ \
--no_secondary \
--genome CAST-EiJ \
--samples Disteche_sample_1 Disteche_sample_2 Disteche_sample_3 Disteche_sample_4 Disteche_sample_5 Disteche_sample_6 Disteche_sample_7 Disteche_sample_8 Disteche_sample_9 Disteche_sample_10 Disteche_sample_11 Disteche_sample_12 Disteche_sample_13 Disteche_sample_14 Disteche_sample_15 Disteche_sample_16 Disteche_sample_17 Disteche_sample_18 Disteche_sample_19 Disteche_sample_20 Disteche_sample_21 Disteche_sample_22 > src/2022-02-23-CAST/output_2022-02-23-CAST-2 2> src/2022-02-23-CAST/err_2022-02-23-CAST-2

./bin/sciatac_pipeline/analyze_sciatac \
--demux_output data/demux_run1 data/demux_run2 data/demux_run3 \
--outdir results/2022-02-21/ \
--no_secondary \
--genome mm10 \
--samples Disteche_sample_1 Disteche_sample_2 Disteche_sample_3 Disteche_sample_4 Disteche_sample_5 Disteche_sample_6 Disteche_sample_7 Disteche_sample_8 Disteche_sample_9 Disteche_sample_10 Disteche_sample_11 Disteche_sample_12 Disteche_sample_13 Disteche_sample_14 Disteche_sample_15 Disteche_sample_16 Disteche_sample_17 Disteche_sample_18 Disteche_sample_19 Disteche_sample_20 Disteche_sample_21 Disteche_sample_22 > src/2022-02-21-sciatac-mm10/output_2022-02-21_mm10 2> src/2022-02-21-sciatac-mm10/err_2022-02-21_mm10

