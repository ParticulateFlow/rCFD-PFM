#!/bin/csh -f

set EXPECTED_ARGS=1
if ( $#argv < ${EXPECTED_ARGS} ) then
    echo "Usage: `basename $0 $1` {target directory}"
    exit
endif

# name of target directory
set TARGETDIR=$1

# check if target directory exists
if ( ! -d "${TARGETDIR}" ) then
    echo "Directory ${TARGETDIR} does not exist, exiting"
    exit
endif

# copy src folder to target dir
echo "Copying src to target directory ..."
cp -r ./src ${TARGETDIR}

# copy balance_check.m to target dir
cp -f "tutorials/balance_check.m" "${TARGETDIR}/tutorials/balance_check.m" >& /dev/null

# switch to tutorial folder
pushd tutorials >/dev/null

# loop over all case folders
echo "Copying run scripts to target directory ..."
foreach d ( */ )
    if ( -f "${d}run_batch.scm" ) then
        if ( "${TARGETDIR}" =~ /* ) then
            # absolute path
            cp -f "${d}run_batch.scm" "${TARGETDIR}/tutorials/${d}run_batch.scm" >& /dev/null
        else
            # relative path
            cp -f "${d}run_batch.scm" "../${TARGETDIR}/tutorials/${d}run_batch.scm" >& /dev/null
        endif
    endif
end

popd >/dev/null

# switch to target dir and execute tutorials
echo "Switching to target directory ..."
pushd ${TARGETDIR}/tutorials >/dev/null

# remove old log file
set LOGFILE="balance_check.log"
rm ${LOGFILE} >& /dev/null

set NSKIPPED=0
set NOK=0
set NINCONSISTENT=0
set NFAILED=0

# run all tutorials
foreach d ( */ )
    pushd ${d} >/dev/null
    if ( -f "run_batch.scm" ) then
        echo "Executing case ${d} ..."

        command -v fluent
        set FLUENT_CHECK=$status

        if ( ${FLUENT_CHECK} == 0) then
            fluent 3ddp -t2 -g < run_batch.scm | tee run_batch.trn
        else
            echo "Failed to find fluent installation ..."
        endif

        # post-process balance
        if ( -f "post/balance_monitor.out" ) then

            set CONFIDENCE=0

            command -v matlab
            set MATLAB_CHECK=$status
            if ( ${MATLAB_CHECK} == 0 ) then
                cp ../balance_check.m .
                matlab -noFigureWindows -batch "balance_check"
                set CONFIDENCE=$status
            else
                command -v octave
                set OCTAVE_CHECK=$status
                if ( ${OCTAVE_CHECK} == 0 ) then
                    octave --silent ../balance_check.m
                    set CONFIDENCE=$status
                else
                    echo "Failed to find matlab/octave for post-processing ..."
                endif
            endif

            echo "Confidence for case ${d}: ${CONFIDENCE}%" | tee -a ../${LOGFILE}

            if ( ${CONFIDENCE} < 90 ) then
                echo "Case ${d} INCONSISTENT" | tee -a ../${LOGFILE}
                @ NINCONSISTENT++
            else
                echo "Case ${d} OK" | tee -a ../${LOGFILE}
                @ NOK++
            endif
        else
            echo "Case ${d} FAILED" | tee -a ../${LOGFILE}
            @ NFAILED++
        endif
    else
        echo "Case ${d} SKIPPED" | tee -a ../${LOGFILE}
        @ NSKIPPED++
    endif
    popd >/dev/null
end

echo "-----------"
echo "RUN SUMMARY"
echo "-----------"

set PLURALS="s"
if (${NSKIPPED} == 1) set PLURALS=" "
echo "$NSKIPPED case$PLURALS SKIPPED"

set PLURALS="s"
if (${NOK} == 1) set PLURALS=" "
echo "$NOK case$PLURALS OK"

set PLURALS="s"
if (${NINCONSISTENT} == 1) set PLURALS=" "
echo "$NINCONSISTENT case$PLURALS INCONSISTENT"

set PLURALS="s"
if (${NFAILED} == 1) set PLURALS=" "
echo "$NFAILED case$PLURALS FAILED"

popd >/dev/null


