#!/bin/bash

make
mpiexec -n 4 ./body3 tests/particles_in.txt particles_out 1 0.5 > particles_out
python3 results_checker.py particles_out tests/particles_out.txt

mpiexec -n 4 ./body3 tests/particles_in_5_steps.txt particles_out 5 1 > particles_out
python3 results_checker.py particles_out tests/particles_out_5_steps.txt

python3 sequential.py tests/particles_in.txt 1 0.5 > python_out_1
mpiexec -n 4 ./body3 tests/particles_in.txt python_out_1 1 0.5 > particles_out
python3 results_checker.py particles_out python_out_1

python3 sequential.py tests/particles_in_5_steps.txt 5 1 > python_out_2
mpiexec -n 4 ./body3 tests/particles_in_5_steps.txt particles_out 5 1 > particles_out
python3 results_checker.py particles_out python_out_2

python3 sequential.py tests/particles_in.txt 1 1 > python_out_3
mpiexec -n 4 ./body3 tests/particles_in.txt particles_out 1 1 > particles_out
python3 results_checker.py particles_out python_out_3

python3 sequential.py tests/particles_in.txt 10 1 > python_out_4
mpiexec -n 4 ./body3 tests/particles_in.txt particles_out 10 1 > particles_out
python3 results_checker.py particles_out python_out_4

python3 sequential.py tests/particles_in.txt 100 1 > python_out_5
mpiexec -n 4 ./body3 tests/particles_in.txt particles_out 100 1 > particles_out
python3 results_checker.py particles_out python_out_5

python3 sequential.py tests/particles_in.txt 1000 1 > python_out_6
mpiexec -n 4 ./body3 tests/particles_in.txt particles_out 1000 1 > particles_out
python3 results_checker.py particles_out python_out_6

python3 sequential.py tests/particles_in_5_steps.txt 1000 0.5 > python_out_7
mpiexec -n 4 ./body3 tests/particles_in_5_steps.txt particles_out 1000 0.5 > particles_out
python3 results_checker.py particles_out python_out_7

rm particles_out python_out*
