import sys
import os
import random

number_of_particles = int(sys.argv[1])
amplitude = float(sys.argv[2])
precision = int(sys.argv[3])
filepath = sys.argv[4]

def random_number():
    return round(random.uniform(-amplitude, amplitude), precision)

positions = set()
with open(filepath, "w") as fp:
    for i in range(number_of_particles):
        while True:
            x = random_number()
            y = random_number()
            z = random_number()
            if not (x, y, z) in positions:
                positions.add((x, y, z))
                break
        vx = random_number()
        vy = random_number()
        vz = random_number()
        fp.write("{} {} {} {} {} {}\n".format(x, y, z, vx, vy, vz))