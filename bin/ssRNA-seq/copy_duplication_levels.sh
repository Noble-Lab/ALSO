#!/bin/bash

#  copy_duplication_levels.sh
#  KA

cd "./data" &&
    {
        mkdir -p "./fastqc/_duplicates"

        unset infiles
        typeset -a infiles
        while IFS=" " read -r -d $'\0'; do
            infiles+=( "${REPLY}" )
        done < <(
            for i in */*/*/duplication_levels.png; do printf "%s\0" "${i}"; done
        )

        for i in "${infiles[@]}"; do
            echo "${i}" | tr "/" "."
            cp "${i}" "./fastqc/_duplicates/$(echo "${i}" | tr "/" ".")"
        done
    }
