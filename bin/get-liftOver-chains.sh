#!/bin/bash

#  get-liftOver-chains.sh
#  KA

#  Start recording time -------------------------------------------------------
start="$(date +%s)"


#  Functions ------------------------------------------------------------------
# loopEcho() { for i in "${@:-*}"; do echo "${i}"; done; }


#  Handle arguments, assign variables -----------------------------------------
#TODO 1/2 Convert from undocumented positional parameters to documented user-
#TODO 2/2 specified, documented keyword parameters

#  In the associative array below, the keys are made up of two "subterms"
#+ separated by a space: '%\ *' matches the first term (to the left of the
#+ space), '#*\ ' matches the second term (to the right of the space)
unset array_chain
typeset -A array_chain=(
    ["https://hgdownload.soe.ucsc.edu/hubs/GCA/001/624/445/GCA_001624445.1/liftOver/GCA_001624445.1ToMm10.over.chain.gz GCA_001624445.1ToMm10.over.chain.gz"]="CAST-EiJ-to-mm10.over.chain.gz"
    ["https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToGCA_001624445.1.over.chain.gz mm10ToGCA_001624445.1.over.chain.gz"]="mm10-to-CAST-EiJ.over.chain.gz"
    ["https://hgdownload.soe.ucsc.edu/hubs/GCA/001/624/185/GCA_001624185.1/liftOver/GCA_001624185.1ToMm10.over.chain.gz GCA_001624185.1ToMm10.over.chain.gz"]="129S1-SvImJ-to-mm10.over.chain.gz"
    ["https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToGCA_001624185.1.over.chain.gz mm10ToGCA_001624185.1.over.chain.gz"]="mm10-to-129S1-SvImJ.over.chain.gz"
    ["https://hgdownload.soe.ucsc.edu/goldenPath/GCA_001624865.1_SPRET_EiJ_v1/liftOver/GCA_001624865.1_SPRET_EiJ_v1ToMm10.over.chain.gz GCA_001624865.1_SPRET_EiJ_v1ToMm10.over.chain.gz"]="SPRET-EiJ-to-mm10.over.chain.gz"
    ["https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToGCA_001624865.1_SPRET_EiJ_v1.over.chain.gz mm10ToGCA_001624865.1_SPRET_EiJ_v1.over.chain.gz"]="mm10-to-SPRET-EiJ.over.chain.gz"
    ["https://hgdownload.soe.ucsc.edu/goldenPath/GCF_900094665.1/liftOver/GCF_900094665.1ToMm10.over.chain.gz GCF_900094665.1ToMm10.over.chain.gz"]="CAROLI-EiJ-to-mm10.over.chain.gz"
    ["https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToGCF_900094665.1.over.chain.gz mm10ToGCF_900094665.1.over.chain.gz"]="mm10-to-CAROLI-EiJ.over.chain.gz"
)
path_data="${1:-"./data/files_chain"}"
rename="${2:-"TRUE"}"


#  Download files -------------------------------------------------------------
#TODO Remove 'for loop' and, instead, download in parallel
for i in "${!array_chain[@]}"; do
    echo "Working with..."
    echo "    - key 1: ${i%\ *}"
    echo "    - key 2: ${i#*\ }"
    echo "    - value: ${array_chain[${i}]}"
    echo ""
    echo "Started: Processing ${array_chain[${i}]}... "

    curl "${i%\ *}" > "${path_data}/${i#*\ }"
    
    if [[ -f "${path_data}/${i#*\ }" ]]; then
        echo "Successfully downloaded ${i#*\ }"
    else
        echo "Exiting: Download of ${i#*\ } was unsuccessful."
    fi

    case "$(echo "${rename}" | tr '[:upper:]' '[:lower:]')" in
        true | t) \
            echo -e "Renaming ${i#*\ } to ${array_chain[${i}]}."
            mv "${path_data}/${i#*\ }" "${path_data}/${array_chain[${i}]}"
            ;;
        false | f) \
            :
            ;;
        *) \
            echo -e "Exiting: \$rename must be \"TRUE\" or \"FALSE\".\n"
            exit 1
            ;;
    esac
    
    echo "Completed: Processing ${array_chain[${i}]}... "
    echo ""
done


#  End recording time ---------------------------------------------------------
end="$(date +%s)"

#  Return run time
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Completed: ${0}"
echo "${0} run time: ${run_time} seconds."
echo ""

exit 0
