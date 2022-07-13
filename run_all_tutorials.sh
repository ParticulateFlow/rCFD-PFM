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

# switch to tutorial folder
pushd tutorials >/dev/null

# copy balance_check.m to target dir
if [[ "${TARGETDIR}" = /* ]]; then
    # absolute path
    cp -f "balance_check.m" "${TARGETDIR}/tutorials/balance_check.m" 2>/dev/null
else
    # relative path
    cp -f "balance_check.m" "../${TARGETDIR}/tutorials/balance_check.m" 2>/dev/null
fi

# loop over all case folders
echo "Copying run scripts to target directory ..."
for d in */ ; do
    if [ -f "${d}run_batch.scm" ]; then
        if [[ "${TARGETDIR}" = /* ]]; then
            # absolute path
            cp -f "${d}run_batch.scm" "${TARGETDIR}/tutorials/${d}run_batch.scm" 2>/dev/null
        else
            # relative path
            cp -f "${d}run_batch.scm" "../${TARGETDIR}/tutorials/${d}run_batch.scm" 2>/dev/null
        fi
    fi
done

popd >/dev/null

# switch to target dir and execute tutorials
echo "Switching to target directory ..."
pushd ${TARGETDIR}/tutorials >/dev/null

BRED='\033[1;31m'
BGREEN='\033[1;32m'
BYELLOW='\033[1;33m'
NC='\033[0m' # No Color
shopt -s expand_aliases
alias decolorize='sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"'

# remove old log file
rm balance_check.log 2> /dev/null

NSKIPPED=0
NOK=0
NINCONSISTENT=0
NFAILED=0

# run all tutorials
for d in */ ; do
    pushd ${d} >/dev/null
    if [ -f "run_batch.scm" ]; then
        echo "Executing case ${d%/} ..."
        fluent 3ddp -t2 -g < run_batch.scm | tee run_batch.trn

        # post-process balance
        if [ -f "post/balance_monitor.out" ]; then
            CONFIDENCE=`octave ../balance_check.m`

            if [ $CONFIDENCE -lt 90 ]; then
                echo -e "Case ${d%/} ${BRED}INCONSISTENT${NC}" | tee >(decolorize >> ../balance_check.log)
                ((NINCONSISTENT++))
            else
                echo -e "Case ${d%/} ${BGREEN}OK${NC}" | tee >(decolorize >> ../balance_check.log)
                ((NOK++))
            fi
        else
            echo -e "Case ${d%/} ${BRED}FAILED${NC}" | tee >(decolorize >> ../balance_check.log)
            ((NFAILED++))
        fi
    else
        echo -e "Case ${d%/} ${BYELLOW}SKIPPED${NC}" | tee >(decolorize >> ../balance_check.log)
        ((NSKIPPED++))
    fi
    popd >/dev/null
done

echo "-------"
echo "SUMMARY"
echo "-------"
PLURALS=$([ $NSKIPPED -eq 1 ] && echo " " || echo "s")
echo -e "$NSKIPPED case$PLURALS ${BYELLOW}SKIPPED${NC}"
PLURALS=$([ $NOK -eq 1 ] && echo " " || echo "s")
echo -e "$NOK case$PLURALS ${BGREEN}OK${NC}"
PLURALS=$([ $NINCONSISTENT -eq 1 ] && echo " " || echo "s")
echo -e "$NINCONSISTENT case$PLURALS ${BRED}INCONSISTENT${NC}"
PLURALS=$([ $NFAILED -eq 1 ] && echo " " || echo "s")
echo -e "$NFAILED case$PLURALS ${BRED}FAILED${NC}"

popd >/dev/null


