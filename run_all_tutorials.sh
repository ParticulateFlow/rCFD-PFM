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

RUN_IN_PLACE=false
if [[ "${TARGETDIR}" -ef ./ ]]; then
    RUN_IN_PLACE=true
fi

if [ "${RUN_IN_PLACE}" = false ]; then
    # copy src folder to target dir
    echo "Copying src to target directory ..."
    cp -r ./src ${TARGETDIR}

    # copy balance_check.m to target dir
    cp -f "tutorials/balance_check.m" "${TARGETDIR}/tutorials/balance_check.m" 2> /dev/null

    # switch to tutorial folder
    pushd tutorials >/dev/null

    # loop over all case folders
    echo "Copying user_src and run scripts to target directory ..."
    for d in */ ; do
        if [ -f "${d}run_batch.scm" ]; then
            if [[ "${TARGETDIR}" = /* ]]; then
                # absolute path
                cp -r "${d}user_src/." "${TARGETDIR}/tutorials/${d}user_src/"
                cp -f "${d}run_batch.scm" "${TARGETDIR}/tutorials/${d}run_batch.scm" 2>/dev/null
            else
                # relative path
                cp -r "${d}user_src/." "../${TARGETDIR}/tutorials/${d}user_src/"
                cp -f "${d}run_batch.scm" "../${TARGETDIR}/tutorials/${d}run_batch.scm" 2>/dev/null
            fi
        fi
    done

    popd >/dev/null

    # switch to target dir and execute tutorials
    echo "Switching to target directory ..."
    pushd ${TARGETDIR}/tutorials >/dev/null
else
    # switch to tutorials dir and execute tutorials
    echo "Switching to tutorials directory ..."
    pushd tutorials >/dev/null
fi

BRED='\033[1;31m'
BGREEN='\033[1;32m'
BYELLOW='\033[1;33m'
NC='\033[0m' # No Color
shopt -s expand_aliases
alias decolorize='sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"'

# remove old log file
LOGFILE="balance_check.log"
rm ${LOGFILE} 2> /dev/null

NSKIPPED=0
NOK=0
NINCONSISTENT=0
NFAILED=0

# run all tutorials
for d in */ ; do
    pushd ${d} >/dev/null
    if [ -f "run_batch.scm" ]; then
        if [ -f "./data/c2c/c2c_0_0_0" ]; then
            command -v fluent >/dev/null
            FLUENT_CHECK=$?

            if [ ${FLUENT_CHECK} -eq 0 ]; then
                # get number of processors from number of files for state0, phase0
                pushd ./data/c2c/ >/dev/null
                NPROCESSORS=`ls -1q c2c_0_0_* | wc -l`
                popd >/dev/null
                echo "Executing case ${d%/} on ${NPROCESSORS} processors ..."
                fluent 3ddp -t${NPROCESSORS} -g < run_batch.scm | tee run_batch.trn
            else
                echo "Failed to find fluent installation ..."
            fi

            # post-process balance
            if [ -f "post/balance_monitor.out" ]; then

                CONFIDENCE=0

                command -v matlab >/dev/null
                MATLAB_CHECK=$?
                if [ ${MATLAB_CHECK} -eq 0 ]; then
                    cp ../balance_check.m .
                    matlab -noFigureWindows -batch "balance_check"
                    CONFIDENCE=$?
                else
                    command -v octave >/dev/null
                    OCTAVE_CHECK=$?
                    if [ ${OCTAVE_CHECK} -eq 0 ]; then
                        octave --silent ../balance_check.m
                        CONFIDENCE=$?
                    else
                        echo "Failed to find matlab/octave for post-processing ..."
                    fi
                fi

                echo "Confidence for case ${d%/}: ${CONFIDENCE}%" | tee -a ../${LOGFILE}

                if [ $CONFIDENCE -lt 90 ]; then
                    echo -e "Case ${d%/} ${BRED}INCONSISTENT${NC}" | tee >(decolorize >> ../${LOGFILE})
                    ((NINCONSISTENT++))
                else
                    echo -e "Case ${d%/} ${BGREEN}OK${NC}" | tee >(decolorize >> ../${LOGFILE})
                    ((NOK++))
                fi
            else
                echo -e "Case ${d%/} ${BRED}FAILED${NC}" | tee >(decolorize >> ../${LOGFILE})
                ((NFAILED++))
            fi
        else
            echo -e "Case ${d%/} ${BYELLOW}SKIPPED${NC} (missing data)" | tee >(decolorize >> ../${LOGFILE})
            ((NSKIPPED++))
        fi
    else
        echo -e "Case ${d%/} ${BYELLOW}SKIPPED${NC}" | tee >(decolorize >> ../${LOGFILE})
        ((NSKIPPED++))
    fi
    popd >/dev/null
done

echo "-----------"
echo "RUN SUMMARY"
echo "-----------"

PLURALS=$([ $NSKIPPED -eq 1 ] && echo " " || echo "s")
echo -e "$NSKIPPED case$PLURALS ${BYELLOW}SKIPPED${NC}"
PLURALS=$([ $NOK -eq 1 ] && echo " " || echo "s")
echo -e "$NOK case$PLURALS ${BGREEN}OK${NC}"
PLURALS=$([ $NINCONSISTENT -eq 1 ] && echo " " || echo "s")
echo -e "$NINCONSISTENT case$PLURALS ${BRED}INCONSISTENT${NC}"
PLURALS=$([ $NFAILED -eq 1 ] && echo " " || echo "s")
echo -e "$NFAILED case$PLURALS ${BRED}FAILED${NC}"

popd >/dev/null


