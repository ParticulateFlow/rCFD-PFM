#!/bin/bash

EXPECTED_ARGS=1
if [ $# -lt ${EXPECTED_ARGS} ]; then
    echo "Usage: `basename $0 $1` {target directory}"
    exit
fi

# name of target directory
TARGETDIR=$1

# check if target directory exists
if [ ! -d "${TARGETDIR}" ]; then
    echo "Directory ${TARGETDIR} does not exist, exiting"
    exit
fi

DOWNLOAD_IN_PLACE=false
if [[ "${TARGETDIR}" -ef ./ ]]; then
    DOWNLOAD_IN_PLACE=true
fi

if [ "${DOWNLOAD_IN_PLACE}" = false ]; then
    # switch to target dir
    echo "Switching to target directory ..."
    pushd ${TARGETDIR}/tutorials >/dev/null
else
    # switch to tutorials dir
    echo "Switching to tutorials directory ..."
    pushd tutorials >/dev/null
fi

# remove old log file
LOGFILE="download.log"
rm ${LOGFILE} 2> /dev/null

BRED='\033[1;31m'
BGREEN='\033[1;32m'
BYELLOW='\033[1;33m'
NC='\033[0m' # No Color

shopt -s expand_aliases
alias decolorize='sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"'

NSKIPPED=0
NOK=0
NFAILED=0

# download all tutorial files
if [[ $HOST == node* ]]; then
    echo "Unable to download case files on nodes ..."
else
    echo "Downloading case files ..."
    for d in */ ; do
        pushd ${d} >/dev/null
        if [ -f "prep_batch.scm" ]; then
            if [ -f "ansys_fluent/README.md" ]; then
                pushd "ansys_fluent" >/dev/null

                CASFILE=`sed -n 's/^[ \t]*//;/cas.h5$/ p' README.md`
                DATFILE=`sed -n 's/^[ \t]*//;/dat.h5$/ p' README.md`

                echo "Downloading Ansys Fluent cas & dat files for case ${d%/} ..."
                wget -nv -N ${CASFILE}

                if [ $? -ne 0 ]; then
                    echo -e "Download of cas file for case ${d%/} ${BRED}FAILED${NC}" | tee >(decolorize >> ../../${LOGFILE})
                    ((NFAILED++))
                else
                    echo -e "Download of cas file for case ${d%/} ${BGREEN}OK${NC}" | tee >(decolorize >> ../../${LOGFILE})
                    ((NOK++))
                fi

                wget -nv -N ${DATFILE}
                if [ $? -ne 0 ]; then
                    echo -e "Download of dat file for case ${d%/} ${BRED}FAILED${NC}" | tee >(decolorize >> ../../${LOGFILE})
                    ((NFAILED++))
                else
                    echo -e "Download of dat file for case ${d%/} ${BGREEN}OK${NC}" | tee >(decolorize >> ../../${LOGFILE})
                    ((NOK++))
                fi

                popd >/dev/null

                if [ -f "data/csv/README.md" ]; then
                    pushd "data/csv" >/dev/null

                    CSVFILE=`sed -n 's/^[ \t]*//;/.zip$/ p' README.md`

                    echo "Downloading CSV files for case ${d%/} ..."
                    wget -nv -N ${CSVFILE}

                    if [ $? -ne 0 ]; then
                        echo -e "Download of csv zip file for case ${d%/} ${BRED}FAILED${NC}" | tee >(decolorize >> ../../${LOGFILE})
                        ((NFAILED++))
                    else
                        echo -e "Download of csv zip file for case ${d%/} ${BGREEN}OK${NC}" | tee >(decolorize >> ../../${LOGFILE})
                        unzip -o -q '*.zip'
                        ((NOK++))
                    fi

                    popd >/dev/null

                fi
            else
                echo -e "Download of files for case ${d%/} ${BYELLOW}SKIPPED${NC}" | tee >(decolorize >> ../${LOGFILE})
                ((NSKIPPED++))
            fi
        else
            echo -e "Download of files for case ${d%/} ${BYELLOW}SKIPPED${NC}" | tee >(decolorize >> ../${LOGFILE})
            ((NSKIPPED++))
        fi
        popd >/dev/null
    done

    echo "----------------"
    echo "DOWNLOAD SUMMARY"
    echo "----------------"

    PLURALS=$([ $NSKIPPED -eq 1 ] && echo " " || echo "s")
    echo -e "$NSKIPPED file download$PLURALS ${BYELLOW}SKIPPED${NC}"
    PLURALS=$([ $NOK -eq 1 ] && echo " " || echo "s")
    echo -e "$NOK file download$PLURALS ${BGREEN}OK${NC}"
    PLURALS=$([ $NFAILED -eq 1 ] && echo " " || echo "s")
    echo -e "$NFAILED file download$PLURALS ${BRED}FAILED${NC}"
fi

popd >/dev/null


