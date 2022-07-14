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

# copy src folder to target dir
echo "Copying src to target directory ..."
cp -r ./src ${TARGETDIR}

# switch to target dir
echo "Switching to target directory ..."
pushd ${TARGETDIR}/tutorials >/dev/null

# remove old log file
LOGFILE="prep.log"
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
echo "Downloading case files ..."
for d in */ ; do
    pushd ${d} >/dev/null
    if [ -f "prep_batch.scm" ]; then
        if [ -f "ansys_fluent/README.md" ]; then
            pushd "ansys_fluent" >/dev/null

            CASFILE=`sed -n 's/^[ \t]*//;/cas.h5$/ p' README.md`
            DATFILE=`sed -n 's/^[ \t]*//;/dat.h5$/ p' README.md`

            echo "Downloading Ansys Fluent cas & dat files for case ${d%/} ..."
            wget -nv -N --show-progress ${CASFILE}

            if [ $? -ne 0 ]; then
                echo -e "Download of cas file for case ${d%/} ${BRED}FAILED${NC}" | tee >(decolorize >> ../../${LOGFILE})
                ((NFAILED++))
            else
                echo -e "Download of cas file for case ${d%/} ${BGREEN}OK${NC}" | tee >(decolorize >> ../../${LOGFILE})
                ((NOK++))
            fi

            wget -nv -N --show-progress ${DATFILE}
            if [ $? -ne 0 ]; then
                echo -e "Download of dat file for case ${d%/} ${BRED}FAILED${NC}" | tee >(decolorize >> ../../${LOGFILE})
                ((NFAILED++))
            else
                echo -e "Download of dat file for case ${d%/} ${BGREEN}OK${NC}" | tee >(decolorize >> ../../${LOGFILE})
                ((NOK++))
            fi

            popd >/dev/null
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
echo

NSKIPPED=0
NOK=0
NFAILED=0

# prep all tutorials
for d in */ ; do
    pushd ${d} >/dev/null
    if [ -f "prep_batch.scm" ]; then
        echo "Pre-processing case ${d%/} ..."
        fluent 3ddp -t2 -g < prep_batch.scm | tee prep_batch.trn

        # check fluent return code
        if [ $? -ne 0 ]; then
            echo -e "Pre-processing of case ${d%/} ${BRED}FAILED${NC}" | tee >(decolorize >> ../${LOGFILE})
            ((NFAILED++))
        else
            echo -e "Pre-processing of case ${d%/} ${BGREEN}OK${NC}" | tee >(decolorize >> ../${LOGFILE})
            ((NOK++))
        fi
    else
        echo -e "Pre-processing of case ${d%/} ${BYELLOW}SKIPPED${NC}" | tee >(decolorize >> ../${LOGFILE})
        ((NSKIPPED++))
    fi
    popd >/dev/null
done

echo "----------------------"
echo "PRE-PROCESSING SUMMARY"
echo "----------------------"
PLURALS=$([ $NSKIPPED -eq 1 ] && echo " " || echo "s")
echo -e "$NSKIPPED case$PLURALS ${BYELLOW}SKIPPED${NC}"
PLURALS=$([ $NOK -eq 1 ] && echo " " || echo "s")
echo -e "$NOK case$PLURALS ${BGREEN}OK${NC}"
PLURALS=$([ $NFAILED -eq 1 ] && echo " " || echo "s")
echo -e "$NFAILED case$PLURALS ${BRED}FAILED${NC}"

popd >/dev/null


