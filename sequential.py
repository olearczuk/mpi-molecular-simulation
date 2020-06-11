import argparse
import math

eps = 4.69041575982343e-08
min_abs = 10e-10
# m = 1.66053892173e-27
m = 1

class Molecule(object):
    def __init__(self, line):
        vars = [float(n) for n in line.split(' ')]
        self.x = vars[0:3]
        self.v = vars[3:6]
        self.a = [0.] * 3
        self.da = [0.] * 3
        self.dv = [0.] * 3
        self.dx = [0.] * 3
        self.f = [0.] * 3

    def update(self):
        for i in range(len(self.x)):
            self.a[i] = self.da[i]
            self.f[i] = 0.


def calc_potential(i, j, k):
    r_i_j = calc_distance(i, j)
    r_i_k = calc_distance(i, k)
    r_k_j = calc_distance(k, j)

    return ((1 / ((r_i_j * r_i_k * r_k_j) ** 3))
            + ((3
                * (-(r_i_j ** 2) + (r_i_k ** 2) + (r_k_j ** 2))
                * ((r_i_j ** 2) - (r_i_k ** 2) + (r_k_j ** 2))
                * ((r_i_j ** 2) + (r_i_k ** 2) - (r_k_j ** 2)))
               / (8 * ((r_i_j * r_i_k * r_k_j) ** 5))))


def calc_derivate(i, j, k):
    res = [0] * len(i.x)

    for idx, val in enumerate(i.x):
        h = eps * val
        if abs(val) < min_abs:
            h = eps * min_abs

        i.x[idx] = val + h
        fx_p_h = calc_potential(i, j, k)
        i.x[idx] = val - h
        fx_m_h = calc_potential(i, j, k)
        i.x[idx] = val

        res[idx] = (fx_p_h - fx_m_h) / ((val + h) - (val - h))

    return res


def calc_force(i, molecules):
    for j, m_j in enumerate(molecules):
        for k, m_k in enumerate(molecules):
            if i != j and j != k and i != k:
                res = calc_derivate(molecules[i], m_j, m_k)
                for idx, el in enumerate(res):
                    molecules[i].f[idx] += el


def calc_distance(i, j):
    distance = math.sqrt((i.x[0] - j.x[0]) ** 2 + (i.x[1] - j.x[1]) ** 2 + (i.x[2] - j.x[2]) ** 2)
    return max(distance, min_abs)


def vertel_integration(molecules, step_count, dt):
    for _ in range(step_count):
        for i, el in enumerate(molecules):
            for j in range(len(el.x)):
                el.x[j] = el.x[j] + el.v[j] * dt + ((el.a[j] * (dt ** 2)) / 2)

        for i, el in enumerate(molecules):
            calc_force(i, molecules)

            for j in range(len(el.x)):
                el.da[j] = -(el.f[j] / m)
                el.v[j] = el.v[j] + (((el.a[j] + el.da[j]) * dt) / 2)

        for el in molecules:
            el.update()


def calc_start_v(molecules):
    for i, el in enumerate(molecules):
        calc_force(i, molecules)
        for j in range(len(el.x)):
            el.da[j] = -(el.f[j] / m)

    for el in molecules:
        el.update()


def parse_molecules(filename):
    molecules = []
    with open(filename, 'r') as f:
        for line in f:
            molecules.append(Molecule(line))

    return molecules


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='File with molcules')
    parser.add_argument('steps', help="Number of steps")
    parser.add_argument('delta', help='Delta time')
    args = parser.parse_args()
    steps = int(args.steps)
    delta = float(args.delta)
    molecules = parse_molecules(args.file)
    calc_start_v(molecules)

    # for el in molecules:
    #     print(el.a[0], el.a[1], el.a[2])

    vertel_integration(molecules, steps, delta)

    for el in molecules:
        print(el.x[0], el.x[1], el.x[2], el.v[0], el.v[1], el.v[2])


if __name__ == '__main__':
    main()