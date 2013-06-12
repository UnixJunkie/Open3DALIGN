#!/usr/bin/env bash

# Open3DALIGN self-test
# run after building and installing Open3DALIGN:
# $ ./test.sh
# The output should be "Test completed successfully"
# In case of failure check test/${OUTPUT_FILE}


temp_ref=`mktemp`
temp_test=`mktemp`
VALIDATION_SCRIPT=validation.sh
TEST_COMPLETED_MSG="Validation succeeded"
DATASETS="ace ache bzr cox2 dhfr gpb therm thr"
COPY_FILES=$DATASETS
INPUT_FILE=
OUTPUT_FILE=validation_single.log
TEST_RESULTS=validation_test_results
REFERENCE_RESULTS=validation_reference_results
OPEN3DTOOL=open3dalign


# Before leaving, clean up temporary files and cd
# where the user originally was
clean_exit()
{
  rm -f $temp_ref $temp_test
  cd $old_cwd
  exit $1
}

# In case the test is killed before completion
abrupt_exit()
{
  cat << eof
Test aborted
eof
  clean_exit 1
}

# Get the numerical values and round them up to three decimals
get_atom_based_alignment_val()
{
  grep RMSD | awk '{printf "%.3f\n", $4}'
}


# Save the folder where the user originally was
old_cwd=$PWD
# Catch all premature death signals
trap abrupt_exit SIGTSTP SIGINT SIGTERM SIGKILL
# cd into the "test" folder
cwd=`dirname $0`
if [ -z $cwd ]; then
  cwd=.
fi
cd $cwd
# Make a test_results folder, copy the necessary
# input files there, then run the test in there
rm -rf $TEST_RESULTS
mkdir $TEST_RESULTS
cp -R $COPY_FILES $TEST_RESULTS
cd $TEST_RESULTS

# Run the self-test
O3A_EXE=`which ${OPEN3DTOOL}` ../$VALIDATION_SCRIPT > ${OUTPUT_FILE}
# If the test did not complete successfully,
# print an error message and quit
if (! grep >&/dev/null "$TEST_COMPLETED_MSG" \
  < ${OUTPUT_FILE}); then
  cat << eof
The test case did not complete successfully.
Please check $cwd/${TEST_RESULTS}/${OUTPUT_FILE}
eof
  clean_exit 1
fi
# Get the reference results
get_atom_based_alignment_val \
  < ../${REFERENCE_RESULTS}/${OUTPUT_FILE} > $temp_ref
# Get the results obtained in the current test
get_atom_based_alignment_val \
  < ${OUTPUT_FILE} > $temp_test
# If they differ, print a warning message and exit
if (! diff >&/dev/null $temp_ref $temp_test); then
  cat << eof
$cwd/${OUTPUT_FILE} presents numerical differences
with respect to $cwd/${REFERENCE_RESULTS}/${OUTPUT_FILE}
Please check $cwd/${TEST_RESULTS}/${OUTPUT_FILE}
eof
  clean_exit 1
fi
# The test was OK, remove the test folder and exit
cat << eof
Test completed successfully
eof
cd ..
rm -rf $TEST_RESULTS
clean_exit 0
