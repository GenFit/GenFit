#!/bin/bash

mkdir -p coverage-build
cd coverage-build
cmake -DCMAKE_BUILD_TYPE=DEBUG ..
make -j4 
make tests
cd -

lcov --zerocounters --directory .
lcov --capture --initial --directory . --output-file genfit-coverage

# coverage-build/bin/gtests # Googletest not yet in the master
coverage-build/bin/unitTests

lcov --no-checksum --directory . --capture --output-file genfit-coverage.info
genhtml genfit-coverage.info -o coverage-report
