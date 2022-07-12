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
NC='\033[0m' # No Color

# run all tutorials
for d in */ ; do
    pushd ${d} >/dev/null
    if [ -f "run_batch.scm" ]; then
        echo "Executing tutorial ${d%/}"
        fluent 3ddp -t2 -g < run_batch.scm | tee run_batch.trn

        # post-process balance
        if [ -f "post/balance_monitor.out" ]; then
            CONFIDENCE=`octave ../balance_check.m`

            if [ $CONFIDENCE -lt 90 ]; then
                echo -e "Case ${d%/} ${BRED}INCONSISTENT${NC}"
            else
                echo -e "Case ${d%/} ${BGREEN}OK${NC}"
            fi
        else
            echo -e "Case ${d%/} ${BRED}FAILED${NC}"
        fi
    fi
    popd >/dev/null
done

popd >/dev/null


