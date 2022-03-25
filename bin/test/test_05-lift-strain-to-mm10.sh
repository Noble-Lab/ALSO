#!/bin/bash

#  test_05-lift-strain-to-mm10.sh
#  KA


#  CAST-EiJ test --------------------------------------------------------------
#  POS bed, chr19
start="$(date +%s)"

safe_mode="FALSE"
infile="./data/2022-0324_test_04_chr19/test.CAST-EiJ.300000.chr19.pos.bed"
outpath="./data/2022-0324_test_05_chr19"
strain="CAST-EiJ"
chain="./data/files_chain/CAST-EiJ-to-mm10.over.chain.gz"

bash bin/workflow/05-lift-strain-to-mm10.sh \
-u "${safe_mode}" \
-i "${infile}" \
-o "${outpath}" \
-s "${strain}" \
-c "${chain}"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."  # run time: 3 seconds
echo ""

#  MPOS bed, chr19
start="$(date +%s)"

safe_mode="FALSE"
infile="./data/2022-0324_test_04_chr19/test.CAST-EiJ.300000.chr19.mpos.bed"
outpath="./data/2022-0324_test_05_chr19"
strain="CAST-EiJ"
chain="./data/files_chain/CAST-EiJ-to-mm10.over.chain.gz"

bash bin/workflow/05-lift-strain-to-mm10.sh \
-u "${safe_mode}" \
-i "${infile}" \
-o "${outpath}" \
-s "${strain}" \
-c "${chain}"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."  # run time: 4 seconds
echo ""

#  POS and MPOS bed, all
start="$(date +%s)"

safe_mode="FALSE"
infile="$(find "./data/2022-0324_test_04_all" -name "*.*os.bed" | sort -n)"
outpath="./data/2022-0324_test_05_all"
strain="CAST-EiJ"
chain="./data/files_chain/CAST-EiJ-to-mm10.over.chain.gz"

parallel --header : -k -j 4 \
"bash ./bin/workflow/05-lift-strain-to-mm10.sh \
-u {safe_mode} \
-i {infile} \
-o {outpath} \
-s {strain} \
-c {chain}" \
::: safe_mode "${safe_mode}" \
::: infile "${infile[@]}" \
::: outpath "${outpath}" \
::: strain "${strain}" \
::: chain "${chain}"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."  # run time: 123 seconds
echo ""
