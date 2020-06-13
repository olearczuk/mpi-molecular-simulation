#!/bin/bash
NODES=${1-4}

make clean
make
mpiexec -n "$NODES" ./body3 tests/test_100.txt particles_out 5 0.5
python3 results_checker.py particles_out_5.txt tests/test_100_5_0.5_out

mpiexec -n "$NODES" ./body3 tests/test_100.txt particles_out 50 0.5
python3 results_checker.py particles_out_50.txt tests/test_100_50_0.5_out

mpiexec -n "$NODES" ./body3 tests/test_100.txt particles_out 100 0.5
python3 results_checker.py particles_out_100.txt tests/test_100_100_0.5_out

rm particles_out*

