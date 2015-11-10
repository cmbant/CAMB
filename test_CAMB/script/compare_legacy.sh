#!/bin/bash

# This bash script compares the results obtained by the create_spectra.sh script
# with their legacy version.
# This script relies on numdiff.
#
# Developed by: Marco Raveri (mraveri@sissa.it)

# standard initialization of the script:

# colors:
RED='\033[0;31m'   # Red
GREEN='\033[0;32m' # Green
BLUE='\033[0;34m'  # Blue
NC='\033[0m'       # No Color

# get the path of the script:
SCRIPT_PATH="`dirname \"$0\"`"                  # relative
SCRIPT_PATH="`( cd \"$SCRIPT_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$SCRIPT_PATH" ] ; then
  exit 1
fi

# get the path of TBP:
TEST_DIR=$SCRIPT_PATH/..
CAMB_DIR=$TEST_DIR/..
TEST_PARAMS_DIR=$TEST_DIR/parameters
TEST_LOGS_DIR=$TEST_DIR/logs
TEST_LEGACY_DIR=$TEST_DIR/legacy_results
TEST_RESULTS_DIR=$TEST_DIR/results

# protect against empty legacy:
if [ "$(ls -A $TEST_LEGACY_DIR)" ]; then
    echo &> /dev/null
else
    echo 'No legacy results.'
    exit 1
fi

# start the script:
printf "${GREEN}********************************************\n"
printf "CAMB TBP: comparing legacy spectra\n"
printf "********************************************\n${NC}"

COMPARE_SPECTRA_LOG=$TEST_LOGS_DIR/compare_spectra.log

FAILED_TEST=()
FAILURE_REASON=()
TOTAL_TEST=0

echo > $COMPARE_SPECTRA_LOG

for i in $TEST_LEGACY_DIR/*;
do

  filename=$(basename "$i")
  extension="${filename##*.}"
  filename="${filename%.*}"

  file_1=$i
  file_2=$TEST_RESULTS_DIR/$(basename "$i")

  echo >> $COMPARE_SPECTRA_LOG
  echo "Comparing: "$filename >> $COMPARE_SPECTRA_LOG
  TOTAL_TEST=$((TOTAL_TEST+1))

  if [ -e $file_2 ]
    then

    if ( numdiff $file_1 $file_2 >> $COMPARE_SPECTRA_LOG 2>&1 ) ; then
        printf "${GREEN}OK:${NC} ${filename} passed succesfully.\n"
    else
        printf "${RED}TEST FAILED:${NC} ${filename} diff failure.\n"
        FAILED_TEST+=($filename)
        FAILURE_REASON+=('diff failed')
        continue
    fi
  else
    printf "${RED}TEST FAILED:${NC} ${filename} output file not found.\n"
    FAILED_TEST+=($filename)
    FAILURE_REASON+=('output not found')
    continue
  fi

done;

printf "${GREEN}********************************************\n"
printf "CAMB TBP test report:\n"
printf "********************************************\n${NC}"

printf "\n Passed test: "$(( $TOTAL_TEST-${#FAILED_TEST[@]} ))" / "$TOTAL_TEST"\n"

if [ "${#FAILED_TEST[@]}" -gt "0" ]; then
  printf "\n Failed test:\n\n"
  for ind in `seq 1 ${#FAILED_TEST[@]}`
  do
  echo " " ${FAILED_TEST[$ind-1]} ": " ${FAILURE_REASON[$ind-1]}
  done;
  printf "\n"
else
  printf " All test successfully passed.\n\n"
fi

printf "${GREEN}********************************************\n"
printf "CAMB TBP: done comparing legacy spectra\n"
printf "********************************************\n${NC}"

exit 0
