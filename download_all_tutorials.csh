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

set DOWNLOAD_IN_PLACE=false
set FULL_TARGETDIR=`realpath "${TARGETDIR}"`
set FULL_CURRENTDIR=`pwd`
if ( "${FULL_TARGETDIR}" == "${FULL_CURRENTDIR}" ) then
    set DOWNLOAD_IN_PLACE=true
endif

if ( ${DOWNLOAD_IN_PLACE} == false ) then
    # switch to target dir
    echo "Switching to target directory ..."
    pushd ${TARGETDIR}/tutorials >/dev/null
else
    # switch to tutorials dir
    echo "Switching to tutorials directory ..."
    pushd tutorials >/dev/null
endif

# remove old log file
set LOGFILE="download.log"
rm ${LOGFILE} >& /dev/null

set NSKIPPED=0
set NOK=0
set NFAILED=0

# download all tutorial files
if ( "${HOST}" =~ node* ) then
    echo "Unable to download case files on nodes ..."
else
    echo "Downloading case files ..."
    foreach d ( */ )
        pushd ${d} >/dev/null
        if ( -f "prep_batch.scm" ) then
            if ( -f "ansys_fluent/README.md" ) then
                pushd "ansys_fluent" >/dev/null

                set CASFILE=`sed -n 's/^[ \t]*//;/cas.h5$/ p' README.md`
                set DATFILE=`sed -n 's/^[ \t]*//;/dat.h5$/ p' README.md`

                echo "Downloading Ansys Fluent cas & dat files for case ${d} ..."
                wget -nv -N ${CASFILE}

                if ( $status != 0 ) then
                    echo "Download of cas file for case ${d} FAILED" | tee -a ../../${LOGFILE}
                    @ NFAILED++
                else
                    echo "Download of cas file for case ${d} OK" | tee -a ../../${LOGFILE}
                    @ NOK++
                endif

                wget -nv -N ${DATFILE}
                if ( $status != 0 ) then
                    echo "Download of dat file for case ${d} FAILED" | tee -a ../../${LOGFILE}
                    @ NFAILED++
                else
                    echo "Download of dat file for case ${d} OK" | tee -a ../../${LOGFILE}
                    @ NOK++
                endif

                popd >/dev/null

                if ( -f "data/csv/README.md" ) then
                    pushd "data/csv" >/dev/null

                    set CSVFILE=`sed -n 's/^[ \t]*//;/.zip$/ p' README.md`

                    echo "Downloading CSV files for case ${d} ..."
                    wget -nv -N ${CSVFILE}

                    if ( $status != 0 ) then
                        echo "Download of csv zip file for case ${d} FAILED" | tee -a ../../${LOGFILE}
                        @ NFAILED++
                    else
                        echo "Download of csv zip file for case ${d} OK" | tee -a ../../${LOGFILE}
                        unzip -o -q '*.zip'
                        @ NOK++
                    endif

                    popd >/dev/null

                endif
            else
                echo "Download of files for case ${d} SKIPPED" | tee -a ../${LOGFILE}
                @ NSKIPPED++
            endif
        else
            echo "Download of files for case ${d} SKIPPED" | tee -a ../${LOGFILE}
            @ NSKIPPED++
        endif
        popd >/dev/null
    end

    echo "----------------"
    echo "DOWNLOAD SUMMARY"
    echo "----------------"

    set PLURALS="s"
    if (${NSKIPPED} == 1) set PLURALS=" "
    echo "$NSKIPPED file download$PLURALS SKIPPED"

    set PLURALS="s"
    if (${NOK} == 1) set PLURALS=" "
    echo "$NOK file download$PLURALS OK"

    set PLURALS="s"
    if (${NFAILED} == 1) set PLURALS=" "
    echo "$NFAILED file download$PLURALS FAILED"

endif

popd >/dev/null


