#!/bin/bash

#  KA

#  Description:
#+ Auxiliary functions for running driver scripts.
#+
#+ Usage:
#+
#+ Notes:
#+


#  Declare variables
# typeset array_resolution
# typeset contact
typeset deposit
# typeset job
# typeset job_prefix
typeset list_prefix
typeset list_suffix
typeset script_make_list
# typeset replicate
typeset stderr
typeset stdout
typeset user
typeset where
typeset when
typeset where_when_what
typeset withdraw


#  Read in functions needed by subsequent functions
echoErrOut() {
    #  Send the same message(s) to stderr and then subsequently to stdout
    >&2 printf "%s\n" "$*"
    >&1 printf "%s\n" "$*"
}

evaluateGreaterEqual() {
    #TODO Documentation
    printf '%s\n%s\n' "${2}" "${1}" | sort -V -C
}


#  Read in functions
checkDependencyBowtie2() {
    command -v bowtie2 &>/dev/null ||
        {
            echoErrOut "Exiting: bowtie2 not found. Install Bowtie2."
            exit 1
        }
}

checkDependencyZsh() {
    #TODO Documentation
    command -v "${ZSH_NAME}" &>/dev/null ||
        {
            echoErrOut "Exiting: You must run this script with zsh."
            exit 1
        }

    version="$("${ZSH_NAME}" --version | awk -F' ' '{ print $2 }')"
    evaluateGreaterEqual "${version}" "5.6" ||
        {
            echoErrOut "Exiting: Install a version of Zsh >= 5.6."
            exit 1
        }
}

checkDependencyParallel() {
    #TODO Documentation
    command -v parallel &>/dev/null ||
        {
            echoErrOut "Exiting: parallel not found. Install GNU Parallel."
            exit 1
        }

    version="$(parallel --version | head -1 | cut -d" " -f3)"
    evaluateGreaterEqual "${version}" "20200101" ||
        {
            echoErrOut "Exiting: Install a version of GNU Parallel from 2020 or" \
            "later."
            exit 1
        }
}

checkJobs() {
    qstat -u "$(whoami)" -xml | grep JB_name | sed 's#</*JB_name>##g'
}

checkPresenceErr() {
    #  If present, remove "${stderr}" with same name
    [[ ! -e "${where_when_what}/${stderr}" ]] || \
        rm "${where_when_what}/${stderr}"
}

checkPresenceOut() {
    #  If present, remove "${stdout}" with same name
    [[ ! -e "${where_when_what}/${stdout}" ]] || \
        rm "${where_when_what}/${stdout}"
}

checkDirectoryPresenceDeposit() {
    #  Check for the presence of the "${deposit}" directory; if it does not
    #+ exist, then prompt the user to create the "${deposit}" directory or,
    #+ alternatively, assign the current working directory as "${deposit}"
    [[ -d "${deposit}" ]] ||
    {
        echo '[[ -d "${deposit}" ]] returned FALSE.'
        echo "[Y]es or [N]o: Make \"deposit\" directory \"${deposit}\"?"
        read -r reply
        # case "$(echo "${reply}" | tr 'A-Z' 'a-z')" in
        case "$(echo "${reply}" | tr '[:upper:]' '[:lower:]')" in
            n | no) \
                echo "No: Assigning \"\$(pwd)\" to \"\${deposit}\""
                echo "deposit=\"$(pwd)\""
                echo ""
                deposit="$(pwd)"
                ;;
            *) \
                echo "Yes: Making directory assigned to \"\${deposit}\""
                echo "deposit=\"${deposit}\""
                echo ""
                mkdir -p "${deposit}"
                ;;    
        esac
    }
}

checkDirectoryPresenceWhereWhenWhat() {
    #  Create "${where_when_what}" directory if it does not exist
    [[ -d "${where_when_what}" ]] || mkdir -p "${where_when_what}"
}

checkDirectoryPresenceWithdrawl() {
    #  Exit with warning if "${withdraw}" directory does not exist
    [[ -d "${withdraw}" ]] || 
    {
        echoErrOut "Exiting because the following \"withdrawl\" directory" \
        "does not exist:"
        echoErrOut "${withdraw}"
        exit 1
    }
}

createOneLineLists() {
    #TODO Documentation
    typeset -a lists
    while IFS=" " read -r -d $'\0'; do lists+=( "${REPLY}" ); done < <(
        find \
        "./${where:-"log/${user}"}" \
        -maxdepth 2 \
        -name "${list_prefix}*${list_suffix}" \
        -type f \
        -printf "%P\0"
    )  #TODO Need to specify depth here?

    # set -Eeuxo pipefail
    for list in "${lists[@]}"; do
        echo ""
        echo "# -------------------------------------------------------"
        echo "Working with ${list}"
        echo ""
        
        list_name="$(basename "${list}")"

        list_path_withdraw="${where}/${when}_${script_make_list%.sh}"
        list_file_withdraw="${list_path_withdraw}/${list_name}"
        line_path_deposit="${where}/${when}_${script_make_list%.sh}_line-files"
     
        #  Create "${line_path_deposit}" directory if it does not exist
        [[ -d "${line_path_deposit}" ]] || mkdir -p "${line_path_deposit}"

        i=0
        sed 1d "${list_file_withdraw}" | while read -r line; do
            #  Increment with each line
            i=$(( i + 1 ))

            #  File for job submission
            line_file_deposit="${line_path_deposit}/${list_name%.txt}.${i}.txt"

            #  If present, remove infile with header and single-line body
            [[ ! -e "${line_file_deposit}" ]] || rm "${line_file_deposit}"

            #  Generate infile with header and single-line body
            # echo "$(head -n 1 ${list_file_withdraw})" >> "${line_file_deposit}"
            head -n 1 "${list_file_withdraw}" >> "${line_file_deposit}"
            echo "${line}" >> "${line_file_deposit}"

            echo "Created file: ${line_file_deposit}"
        done

        echo ""
    done
}

findMultiLineLists() {
    find \
    "./${where}" \
    -name "${list_prefix}*${list_suffix}" \
    -type f | \
        sort -V
}

findOneLineLists() {
    find \
    "./${where}" \
    -name "${list_prefix}*${list_suffix%.txt}.*.txt" \
    -type f | \
        sort -V
}

# findOneLineListsIndividualPooled() {
#     find "./${where}" -maxdepth 2 -name "*_*.${list_prefix}*${list_suffix%.txt}.*.txt" -type f | sort -V
# }

printSaveErr() {
    #  Print stderr to a path and file: "${where_when_what}/${stderr}"
    #+ 
    #+ If the process-substitution one-liner is commented out (and the other
    #+ one-liner is uncommented), then stderr is also printed to the screen
    exec 2> "${where_when_what}/${stderr}"
    # exec 2> >(tee -a "${where_when_what}/${stderr}")
}

printSaveOut() {
    #  Print stdout to a path and file: "${where_when_what}/${stdout}"
    #+ 
    #+ If the process-substitution one-liner is commented out (and the other
    #+ one-liner is uncommented), then stdout is also printed to the screen
    # exec 1> "${where_when_what}/${stdout}"
    exec 1> >(tee -a "${where_when_what}/${stdout}")
}
