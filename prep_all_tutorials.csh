#!/bin/csh

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

# switch to tutorial folder
pushd tutorials >/dev/null

# loop over all case folders
echo "Copying pre-processing scripts to target directory ..."
foreach d ( */ )
    if ( -f "${d}prep_batch.scm" ) then
        if ( "${TARGETDIR}" =~ /* ) then
            # absolute path
            cp -f "${d}prep_batch.scm" "${TARGETDIR}/tutorials/${d}prep_batch.scm" >& /dev/null
        else
            # relative path
            cp -f "${d}prep_batch.scm" "../${TARGETDIR}/tutorials/${d}prep_batch.scm" >& /dev/null
        endif
    endif
end

popd >/dev/null

# switch to target dir
echo "Switching to target directory ..."
pushd ${TARGETDIR}/tutorials >/dev/null

# remove old log file
set LOGFILE="prep.log"
rm ${LOGFILE} >& /dev/null

set BRED='\033[1;31m'
set BGREEN='\033[1;32m'
set BYELLOW='\033[1;33m'
set NC='\033[0m' # No Color

#shopt -s expand_aliases
alias decolorize 'sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"'

set NSKIPPED=0
set NOK=0
set NFAILED=0

# download all tutorial files
echo "Downloading case files ..."
foreach d ( */ )
    pushd ${d} >/dev/null
    if ( -f "prep_batch.scm" ) then
        if ( -f "ansys_fluent/README.md" ) then
            pushd "ansys_fluent" >/dev/null

            set CASFILE=`sed -n 's/^[ \t]*//;/cas.h5$/ p' README.md`
            set DATFILE=`sed -n 's/^[ \t]*//;/dat.h5$/ p' README.md`

            echo "Downloading Ansys Fluent cas & dat files for case ${d} ..."
            wget -nv -N --show-progress ${CASFILE}

            if ( $status != 0 ) then
                echo "Download of cas file for case ${d} FAILED" | tee ../../${LOGFILE}
                @ NFAILED++
            else
                echo "Download of cas file for case ${d} OK" | tee ../../${LOGFILE}
                #printf "Download of cas file for case ${d} ${BGREEN}OK${NC}\n" | tee -a ../../${LOGFILE}
                @ NOK++
            endif

            wget -nv -N --show-progress ${DATFILE}
            if ( $status != 0 ) then
                echo "Download of dat file for case ${d} FAILED" | tee ../../${LOGFILE}
                @ NFAILED++
            else
                echo "Download of dat file for case ${d} OK" | tee ../../${LOGFILE}
                @ NOK++
            endif

            popd >/dev/null
        else
            echo "Download of files for case ${d} SKIPPED" | tee ../${LOGFILE}
            @ NSKIPPED++
        endif
    else
        echo "Download of files for case ${d} SKIPPED" | tee ../${LOGFILE}
        @ NSKIPPED++
    endif
    popd >/dev/null
end

echo "----------------"
echo "DOWNLOAD SUMMARY"
echo "----------------"
set PLURALS=' ' #`($NSKIPPED == 1) && echo ' ' || echo 's'`
echo "$NSKIPPED file download$PLURALS SKIPPED"
#set PLURALS=$([ $NOK == 1 ] && echo " " || echo "s")
echo "$NOK file download$PLURALS OK"
#set PLURALS=$([ $NFAILED == 1 ] && echo " " || echo "s")
echo "$NFAILED file download$PLURALS FAILED"
echo ""


set NSKIPPED=0
set NOK=0
set NFAILED=0

# prep all tutorials
foreach d ( */ )
    pushd ${d} >/dev/null
    if ( -f "prep_batch.scm" ) then
        echo "Pre-processing case ${d} ..."
        #fluent 3ddp -t2 -g < prep_batch.scm | tee prep_batch.trn

        # check fluent return code
        if ( $status != 0 ) then
            #echo -e "Pre-processing of case ${d%/} ${BRED}FAILED${NC}" | tee >(decolorize >> ../${LOGFILE})
            echo "Pre-processing of case ${d} FAILED" | tee ../${LOGFILE}
            @ NFAILED++
        else
            #echo -e "Pre-processing of case ${d%/} ${BGREEN}OK${NC}" | tee >(decolorize >> ../${LOGFILE})
            echo "Pre-processing of case ${d} OK" | tee ../${LOGFILE}
            @ NOK++
        endif
    else
        #echo -e "Pre-processing of case ${d%/} ${BYELLOW}SKIPPED${NC}" | tee >(decolorize >> ../${LOGFILE})
        echo "Pre-processing of case ${d} SKIPPED" | tee ../${LOGFILE}
        @ NSKIPPED++
    endif
    popd >/dev/null
end

echo "----------------------"
echo "PRE-PROCESSING SUMMARY"
echo "----------------------"
set PLURALS=' ' #$([ $NSKIPPED == 1 ] && echo " " || echo "s")
echo "$NSKIPPED case$PLURALS SKIPPED"
#set PLURALS=$([ $NOK == 1 ] && echo " " || echo "s")
echo "$NOK case$PLURALS OK"
#set PLURALS=$([ $NFAILED == 1 ] && echo " " || echo "s")
echo "$NFAILED case$PLURALS FAILED"

popd >/dev/null

