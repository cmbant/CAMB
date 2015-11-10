#!/bin/bash

# This bash script runs CAMB benchmarker for all the parameters files in the folder
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
  # error; for some reason, the path is not accessible
  # to the script (e.g. permissions re-evaled after suid)
  exit 1  # fail
fi

# get the path of the test suite:
TEST_DIR=$SCRIPT_PATH/..
CAMB_DIR=$TEST_DIR/..
TEST_PARAMS_DIR=$TEST_DIR/parameters
TEST_LOGS_DIR=$TEST_DIR/logs
TEST_LEGACY_DIR=$TEST_DIR/legacy_results
TEST_RESULTS_DIR=$TEST_DIR/results

# start the script:

printf "${GREEN}********************************************\n"
printf "CAMB TBP: running benchmarks\n"
printf "********************************************\n${NC}"

BENCHMARKS_LOG=$TEST_LOGS_DIR/run_benchmarks.log

cd $CAMB_DIR

echo "Compiling CAMB"

make clean &> /dev/null
make camb_benchmark  > $BENCHMARKS_LOG 2>&1

for i in $TEST_PARAMS_DIR/*;
do

  filename=$(basename "$i")
  extension="${filename##*.}"
  filename="${filename%.*}"

  echo | tee -a $BENCHMARKS_LOG
  echo "Doing : " $filename | tee -a $BENCHMARKS_LOG

  ./camb_benchmark $i 2 | tee -a $BENCHMARKS_LOG

done;

# cleaning up:
make clean  &> /dev/null
rm -rf camb_benchmark &> /dev/null

printf "${GREEN}********************************************\n"
printf "CAMB TBP: finished benchmarks\n"
printf "********************************************\n${NC}"

exit 0
