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

# run all tutorials
for d in */ ; do
    pushd ${d} >/dev/null
    if [ -f "run_batch.scm" ]; then
        echo "Executing tutorial ${d}"
        fluent 3ddp -t2 -g < run_batch.scm >& run_batch.trn
    fi
    popd >/dev/null
done

popd >/dev/null


