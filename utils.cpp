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
	particle.acc = -1 / unitMass * (particle.plusPotential - particle.minusPotential) / hh;
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