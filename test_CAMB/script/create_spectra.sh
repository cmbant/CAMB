#!/bin/bash

# This bash script runs CAMB for all the parameters files in the folder
# test_CAMB/parameters.
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

# start the script:
printf "${GREEN}********************************************\n"
printf "CAMB TBP: creating test spectra\n"
printf "********************************************\n${NC}"

CREATE_SPECTRA_LOG=$TEST_LOGS_DIR/create_spectra.log

cd $CAMB_DIR

echo "Compiling CAMB"

make clean &> /dev/null
make camb  > $CREATE_SPECTRA_LOG 2>&1

for i in $TEST_PARAMS_DIR/*;
do

  filename=$(basename "$i")
  extension="${filename##*.}"
  filename="${filename%.*}"

  echo | tee -a $CREATE_SPECTRA_LOG
  echo "Doing: " $filename | tee -a $CREATE_SPECTRA_LOG

  ./camb $i | tee -a $CREATE_SPECTRA_LOG &> /dev/null

done;

# cleaning up:
make clean  &> /dev/null
rm -rf camb &> /dev/null

printf "${GREEN}********************************************\n"
printf "CAMB TBP: done creating test spectra\n"
printf "********************************************\n${NC}"

exit 0
