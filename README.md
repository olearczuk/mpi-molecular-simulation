# Overview
Goal of this project was to implement molecular simulation with three-body interaction.
Algorithm uses [Axilrod-Teller Potential](https://en.wikipedia.org/wiki/Axilrod%E2%80%93Teller_potential)
and [verlocity Verlet algorithm](https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet). <br/>
Algorithm is based on embedded algorithm from [A Computation- and Communication-Optimal Parallel Direct 3-Body Algorithm](https://www.researchgate.net/publication/282380541_A_Computation-_and_Communication-Optimal_Parallel_Direct_3-Body_Algorithm)
paper.
# Usage
```bash
./body3 particles_in.txt particles_out stepcount deltatime [-v]
```
where
- `particles_in.txt` defines the initial positions and velocities, the format is one particle per line, each line consists 
of 3 doubles specifying the x,y,z coordinates of a particle (single space separated); then 3 doubles specifying 
the Vx, Vy, Vz velocities of a particle. This file does not specify accelerations, thus start with step 2 of the 3-step velocity Verlet, calculating the potential gradient for inital positions.
- `particles_out` is the base name of the output file. The actual result (same format as the input file) must be saved in a file particles_out_stepcount.txt (please don't be too smart and subsitute stepcount with the actual number…).
- `stepcount` is the total number of steps of the Verlet algorithm.
- `deltatime` gives Δt between steps
- [-v], if present, puts the result after each step `i` (counted from 1) in a file `particles_out_i.txt`

# Solution
To write solution I used MPI library. In order to make communication more optimised I used asynchronous `Isend/Irecv` communication, where possible. <br/>
For broadcasting messages to multiple receivers I used `Scatterv` and `Gatherv` for collecting messages from multiple sources.

