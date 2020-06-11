#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cstdio>
#include "utils.h"

void updateCoords(particle1D &particle, double delta) {
	particle.coord += particle.v * delta + 0.5 * particle.acc * delta * delta;
}

void updateCoords(particle3D &particle, double delta) {
	updateCoords(particle.x, delta);
	updateCoords(particle.y, delta);
	updateCoords(particle.z, delta);
}

void updateCoords(std::vector<particle3D> &v, double delta) {
	for (particle3D &p : v) {
		updateCoords(p, delta);
	}
}

double getH(particle1D &particle) {
	// TODO - should be here abs(particle.coord) or not (shouldn't make a difference here)
    if (std::abs(particle.coord) >= minDistance)
        return particle.coord * auxE;
    else
        return minDistance * auxE;
}

double computeNorm(particle3D &particle1, particle3D &particle2) {
	double x1 = particle1.x.coord, y1 = particle1.y.coord, z1 = particle1.z.coord,
		x2 = particle2.x.coord, y2 = particle2.y.coord, z2 = particle2.z.coord;
	double norm = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
	return fmax(norm, minDistance);
}

double computePotential(particle3D &i, particle3D &j, particle3D &k) {
	double Rij = computeNorm(i, j), Rik = computeNorm(i, k), Rkj = computeNorm(k, j);
	double potential = std::abs(1 / pow(Rij * Rik * Rkj, 3) + 3. *
		(-Rij*Rij + Rik*Rik + Rkj*Rkj) * (Rij*Rij - Rik*Rik + Rkj*Rkj) * (Rij*Rij + Rik*Rik - Rkj*Rkj) /
		(8. * pow(Rij * Rik * Rkj, 5)));
	return 2 * potential;
}

void updatePotential1D(particle1D &i1D, particle3D &i, particle3D &j, particle3D &k) {
	double h = getH(i1D);
	i1D.coord += h;
	i1D.potential += computePotential(i, j, k);
	i1D.coord -= 2 * h;
	i1D.potential -= computePotential(i, j, k);
	i1D.coord += h;
}

void updatePotential(particle3D &i, particle3D &j, particle3D &k) {
	updatePotential1D(i.x, i, j, k);
	updatePotential1D(i.y, i, j, k);
	updatePotential1D(i.z, i, j, k);


	updatePotential1D(j.x, j, i, k);
	updatePotential1D(j.y, j, i, k);
	updatePotential1D(j.z, j, i, k);

	updatePotential1D(k.x, k, j, i);
	updatePotential1D(k.y, k, j, i);
	updatePotential1D(k.z, k, j, i);
}

void updatePotential(std::vector<particle3D> &v1, std::vector<particle3D> &v2, std::vector<particle3D> &v3) {
	for (auto & i : v1) {
		for (auto & j : v2) {
			if (i.number == j.number)
				continue;
			for (auto & k : v3)
				// TODO v3[k].number > v2[j].number (2 * in computePotential)
				if (i.number != k.number && j.number != k.number) {
//					printf("%f %f %f\n", v1[0].number, v2[0].number, v3[0].number);
					updatePotential(i, j, k);
				}
		}
	}
}
void updateAcceleration(particle1D &particle) {
    double h = getH(particle);
    volatile double hh = particle.coord + h - (particle.coord - h);
    particle.acc = -1 / unitMass * particle.potential / hh;
}

void updateVelocity(particle1D &particle, double oldAcc, double delta) {
	particle.v += 0.5 * (oldAcc + particle.acc) * delta;
}

void updateAccelerationVelocity(particle1D &particle, double delta) {
	double oldAcc = particle.acc;
	updateAcceleration(particle);
	updateVelocity(particle, oldAcc, delta);
}

void updateAccelerationVelocity(particle3D &particle, double delta) {
	updateAccelerationVelocity(particle.x, delta);
	updateAccelerationVelocity(particle.y, delta);
	updateAccelerationVelocity(particle.z, delta);
}

void updateAccelerationVelocity(std::vector<particle3D> &v, double delta) {
	for (particle3D &p : v)
		updateAccelerationVelocity(p, delta);
}


void updateAcceleration(particle3D &particle) {
	updateAcceleration(particle.x);
	updateAcceleration(particle.y);
	updateAcceleration(particle.z);
}

void updateAcceleration(std::vector<particle3D> &v) {
	for (particle3D &p : v)
		updateAcceleration(p);
}

void parseCommandLineArgs(int argc, char *argv[], std::string &inFilename, std::string &outFilename, int &steps,
		double &delta, bool &isVerbose) {
	if (argc < 5) {
		printf("Usage: ./body3 <in file> <out file> <stepcount> <deltatime> [-v]\n");
		exit(1);
	}
	inFilename = argv[1];
	outFilename = argv[2];
	steps = atoi(argv[3]);
	delta = atof(argv[4]);
    isVerbose = argc == 6 && strcmp(argv[5], "-v") == 0;
}

// TODO maybe propagate particles in chunks here
std::vector<particle3D> readFile(const std::string& fileName) {
	particle3D p{};
	p.x.acc = p.y.acc = p.z.acc = 0;
	p.x.potential = p.y.potential = p.z.potential = 0;
	std::vector<particle3D> v;
	std::fstream file;
	file.open(fileName);
	std::string s;
	int number = 0;
	while(std::getline(file, s)) {
		std::istringstream stream(s);
		stream >> p.x.coord >> p.y.coord >> p.z.coord >> p.x.v >> p.y.v >> p.z.v;
		p.number = number;
		number++;
		v.emplace_back(p);
	}
	return v;
}

void particlesToArray(std::vector<particle3D> v, int begin, int end, double *arr) {
	int i = 0;
	for (; begin < end; begin++) {
		particle3D p = v[begin];
		arr[i] = p.number;
		i++;
		arr[i] = p.x.coord;
		arr[i+1] = p.x.v;
		arr[i+2] = p.x.acc;
		arr[i+3] = p.x.potential;
		i += 4;
		arr[i] = p.y.coord;
		arr[i+1] = p.y.v;
		arr[i+2] = p.y.acc;
        arr[i+3] = p.y.potential;
        i += 4;
		arr[i] = p.z.coord;
		arr[i+1] = p.z.v;
		arr[i+2] = p.z.acc;
        arr[i+3] = p.z.potential;
        i += 4;
	}
}

std::vector<particle3D> arrayToParticles(const double *arr, int size) {
	std::vector<particle3D> v;
	int i = 0;
	while (i < size) {
		particle3D p{};
		p.number = arr[i];
		i++;
		p.x.coord = arr[i];
		p.x.v = arr[i+1];
		p.x.acc = arr[i+2];
		p.x.potential = arr[i+3];
		i += 4;

		p.y.coord = arr[i];
		p.y.v = arr[i+1];
		p.y.acc = arr[i+2];
        p.y.potential = arr[i+3];
        i += 4;

		p.z.coord = arr[i];
		p.z.v = arr[i+1];
		p.z.acc = arr[i+2];
        p.z.potential = arr[i+3];
        i += 4;
		v.emplace_back(p);
	}
	return v;
}

// TODO remove
void printParticle(particle3D p) {
//    printf("%f: %f %f %f %f %f | %f %f %f %f %f | %f %f %f %f %f\n", p.number,
//            p.x.coord, p.x.v, p.x.acc, p.x.minusPotential, p.x.plusPotential,
//           p.y.coord, p.y.v, p.x.acc, p.y.minusPotential, p.y.plusPotential,
//           p.z.coord, p.z.v, p.z.acc, p.z.minusPotential, p.z.plusPotential);

//    printf("%.0f: %.15f %.15f %.15f\n", p.number, p.x.acc, p.y.acc, p.z.acc);

    printf("%.0f: %.15f %.15f %.15f %.15f %.15f %.15f\n", p.number, p.x.coord, p.y.coord, p.z.coord,
            p.x.v, p.y.v, p.z.v);

//    printf("%f: %f %f %f | %.15f %.15f %.15f\n", p.number,
//           p.x.coord, p.y.coord, p.z.coord,
//           p.x.potential, p.y.potential, p.z.potential);
}