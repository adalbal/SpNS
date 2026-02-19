#!/bin/bash

total_tests=0
failed_tests=0

for X in `ls -d TEST_*` ; do
  ((total_tests++))

  cd $X
  echo "========================================="
  echo $X
  echo "========================================="

  cat stdout.0 | grep "Final " > 1
  cat REF/stdout.0 | grep "Final " > 2

  if ! diff 1 2 > /dev/null 2>&1; then
    ((failed_tests++))
    echo "Test failed! Differences found:"
    diff 1 2
  fi

  rm -rf 1 2
  cd ..
done

echo
echo "======================="
echo "***     SUMMARY     ***"
echo
echo "Total tests run:      $total_tests"
echo "Tests failed:         $failed_tests"
echo
if [ $failed_tests -eq 0 ]; then
  echo "*** QA TETS PASSED! ***"
else
  echo "*** QA TETS FAILED! ***"
fi
echo "======================="
echo
