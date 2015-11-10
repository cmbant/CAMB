#!/bin/bash

# This bash script runs CAMB profiler for all the parameters files in the folder
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

# create the profiler log folder:
rm -rf $TEST_LOGS_DIR/profiler_logs
mkdir  $TEST_LOGS_DIR/profiler_logs
TEST_LOGS_DIR=$TEST_LOGS_DIR/profiler_logs

# start the script:

printf "${GREEN}********************************************\n"
printf "CAMB TBP: running profiler\n"
printf "********************************************\n${NC}"

PROFILER_LOG=$TEST_LOGS_DIR/run_profile.log

cd $CAMB_DIR

echo "Compiling CAMB"

make clean &> /dev/null
make camb_benchmark FFLAGS=-pg > $PROFILER_LOG 2>&1

for i in $TEST_PARAMS_DIR/*;
do

  filename=$(basename "$i")
  extension="${filename##*.}"
  filename="${filename%.*}"

  echo | tee -a $PROFILER_LOG
  echo "Doing : " $filename | tee -a $PROFILER_LOG

  ./camb_benchmark $i 2 | tee -a $PROFILER_LOG &> /dev/null
  gprof -b camb_benchmark > $TEST_LOGS_DIR/profile_$filename.log

done;

# cleaning up:
make clean  &> /dev/null
rm -rf camb_benchmark &> /dev/null

printf "${GREEN}********************************************\n"
printf "CAMB TBP: finished profilings\n"
printf "********************************************\n${NC}"

exit 0
