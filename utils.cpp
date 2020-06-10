#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cstdio>
#include "utils.h"

void updateCoords(particle1D &particle, bool delta) {
	particle.coord += particle.v * delta + 0.5 * particle.acc * delta * delta;
}

void updateCoords(particle3D &particle, bool delta) {
	updateCoords(particle.x, delta);
	updateCoords(particle.y, delta);
	updateCoords(particle.z, delta);
}

double getH(particle1D &particle) {
	// TODO - should be here abs(particle.coord) or not (shouldn't make a difference here)
	return abs(particle.coord) < minDistance ? minDistance * auxE : particle.coord * auxE;
}

double computeNorm(particle3D &particle1, particle3D &particle2) {
	double x1 = particle1.x.coord, y1 = particle1.y.coord, z1 = particle1.z.coord,
		x2 = particle2.x.coord, y2 = particle2.y.coord, z2 = particle2.z.coord;
	double norm = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
	return fmax(norm, minDistance);
}

double computePotential(particle3D &i, particle3D &j, particle3D &k) {
	double Rij = computeNorm(i, j), Rik = computeNorm(i, k), Rkj = computeNorm(k, j);
	return abs(1 / pow(Rij * Rik * Rkj, 3) + 3 *
		(-Rij*Rij + Rik*Rik + Rkj*Rkj) * (Rij*Rij - Rik*Rik + Rkj*Rkj) * (Rij*Rij + Rik*Rik - Rkj*Rkj) /
		8 / pow(Rij * Rik * Rkj, 5));
}

void updatePotential1D(particle1D &i1D, particle3D &i, particle3D &j, particle3D &k) {
	double h = getH(i1D);
	i1D.coord += h;
	i1D.plusPotential += computePotential(i, j, k);
	i1D.coord -= 2 * h;
	i1D.minusPotential += computePotential(i, j, k);
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

void updateAcceleration(particle1D &particle) {
	double h = getH(particle);
	volatile double hh = particle.coord + h - (particle.coord - h);
	particle.acc = -(particle.plusPotential - particle.minusPotential) / hh;
}

void updateVelocity(particle1D &particle, double oldAcc, bool delta) {
	particle.v += 0.5 * (oldAcc + particle.acc) * delta;
}

void updateAccelerationVelocity(particle1D &particle, bool delta) {
	double oldAcc = particle.acc;
	updateAcceleration(particle);
	updateVelocity(particle, oldAcc, delta);
}

void updateAccelerationVelocity(particle3D &particle, bool delta) {
	updateAccelerationVelocity(particle.x, delta);
	updateAccelerationVelocity(particle.y, delta);
	updateAccelerationVelocity(particle.z, delta);
}

void updateAcceleration(particle3D &particle) {
	updateAcceleration(particle.x);
	updateAcceleration(particle.y);
	updateAcceleration(particle.z);
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
	if (argc == 6 && strcmp(argv[5], "-v") == 0)
		isVerbose = true;
	else
		isVerbose = false;
}

// TODO maybe propagate particles in chunks here
std::vector<particle3D> readFile(std::string fileName) {
	particle3D p;
	p.x.acc = p.y.acc = p.z.acc = 0;
	p.x.plusPotential = p.y.plusPotential = p.z.plusPotential = 0;
	p.x.minusPotential = p.y.minusPotential = p.z.minusPotential = 0;
	std::vector<particle3D> v;
	std::fstream file;
	file.open(fileName);
	int V, E;
	std::string s;
	while(std::getline(file, s)) {
		std::istringstream stream(s);
		stream >> p.x.coord >> p.y.coord >> p.z.coord >> p.x.v >> p.y.v >> p.z.v;
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
		// TODO arr[i+3] = p.x.potential; i += 4;
		arr[i+3] = p.x.minusPotential;
		arr[i+4] = p.x.plusPotential;
		i += 5;
		arr[i] = p.y.coord;
		arr[i+1] = p.y.v;
		arr[i+2] = p.y.acc;
		// TODO arr[i+3] = p.y.potential; i += 4;
		arr[i+3] = p.y.minusPotential;
		arr[i+4] = p.y.plusPotential;
		i += 5;
		arr[i] = p.z.coord;
		arr[i+1] = p.z.v;
		arr[i+2] = p.z.acc;
		// TODO arr[i+3] = p.z.potential; i += 4;
		arr[i+3] = p.z.minusPotential;
		arr[i+4] = p.z.plusPotential;
		i += 5;
	}
}

std::vector<particle3D> arrayToParticles(double *arr, int size) {
	std::vector<particle3D> v;
	int i = 0;
	while (i < size) {
		particle3D p;
		p.number = arr[i];
		i++;
		p.x.coord = arr[i];
		p.x.v = arr[i+1];
		p.x.acc = arr[i+2];
		// TODO p.x.potential = arr[i+3]; i += 4;
		p.x.minusPotential = arr[i+3];
		p.x.plusPotential = arr[i+4];
		i += 5;

		p.y.coord = arr[i];
		p.y.v = arr[i+1];
		p.y.acc = arr[i+2];
		// TODO p.y.potential = arr[i+3]; i += 4;
		p.y.minusPotential = arr[i+3];
		p.y.plusPotential = arr[i+4];
		i += 5;

		p.z.coord = arr[i];
		p.z.v = arr[i+1];
		p.z.acc = arr[i+2];
		// TODO p.z.potential = arr[i+3]; i += 4;
		p.z.minusPotential = arr[i+3];
		p.z.plusPotential = arr[i+4];
		i += 5;
		v.emplace_back(p);
	}
	return v;
}