#! /bin/sh
# Test suite script for the Elk Code

for i in test-*
do
  cd $i
  echo
  echo "Running test in directory $i..."
  \rm -f *.OUT
  ../../src/elk > test.log
  NERROR=`grep -c Error test.log`
  if test $NERROR -gt 0
  then
    echo " Failed! See test.log and output files"
  else
    echo " Passed"
    \rm -f test.log
    \rm -f *.OUT
  fi
  cd ..
done

