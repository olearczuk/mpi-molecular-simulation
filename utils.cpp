#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <cstdio>
#include "utils.h"

void updateCoords(particle1D &particle, double delta) {
	particle.coord += particle.v * delta + 0.5 * particle.acc * delta * delta;
}

void updateCoords(std::vector<particle3D> &v, double delta) {
	for (particle3D &p : v) {
		updateCoords(p.x, delta);
        updateCoords(p.y, delta);
        updateCoords(p.z, delta);
	}
}

double getH(particle1D &particle) {
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
	return 1 / pow(Rij * Rik * Rkj, 3) + 3. *
		(-Rij*Rij + Rik*Rik + Rkj*Rkj) * (Rij*Rij - Rik*Rik + Rkj*Rkj) * (Rij*Rij + Rik*Rik - Rkj*Rkj) /
		(8. * pow(Rij * Rik * Rkj, 5));
}

void updatePotential1D(particle1D &i1D, particle3D &i, particle3D &j, particle3D &k) {
    double h = getH(i1D), oldCoord = i1D.coord;
    i1D.coord = oldCoord + h;
    double potentialPlus = computePotential(i, j, k);
    i1D.coord = oldCoord - h;
    double potentialMinus = computePotential(i, j, k);
    i1D.potential += 2 * (potentialPlus - potentialMinus);
    i1D.coord = oldCoord;
}

void updatePotential(std::vector<particle3D> &v1, std::vector<particle3D> &v2, std::vector<particle3D> &v3,
        int owner1, int owner2, int owner3) {
	for (auto & i : v1) {
		for (auto & j : v2) {
			for (auto & k : v3)
			    if ((owner1 != owner2 && owner1 != owner3 && owner2 != owner3) ||
			        (owner1 == owner2 && owner3 != owner1 && i.number < j.number && k.number != i.number && k.number != j.number) ||
			        (owner1 == owner3 && owner2 != owner1 && i.number < k.number && j.number != i.number && j.number != k.number) ||
			        (owner2 == owner3 && owner1 != owner2 && j.number < k.number && i.number != j.number && i.number != k.number) ||
			        (owner1 == owner2 && owner1 == owner3 && i.number < j.number && j.number < k.number)) {
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

void updateAccelerationVelocity(std::vector<particle3D> &v, double delta) {
	for (particle3D &p : v) {
        updateAccelerationVelocity(p.x, delta);
        updateAccelerationVelocity(p.y, delta);
        updateAccelerationVelocity(p.z, delta);
	}
}

void updateAcceleration(std::vector<particle3D> &v) {
	for (particle3D &p : v) {
        updateAcceleration(p.x);
        updateAcceleration(p.y);
        updateAcceleration(p.z);
	}
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
		p.number = number++;
		v.emplace_back(p);
	}
	file.close();
	return v;
}

void particlesToArray(const std::vector<particle3D>& v, double *arr) {
	int i = 0;
	for (particle3D p : v) {
		arr[i] = p.number;
		i++;
        arr[i] = p.x.potential;
		arr[i+1] = p.x.coord;
        arr[i+2] = p.x.acc;
		arr[i+3] = p.x.v;
		i += 4;
        arr[i] = p.y.potential;
		arr[i+1] = p.y.coord;
        arr[i+2] = p.y.acc;
		arr[i+3] = p.y.v;
        i += 4;
        arr[i] = p.z.potential;
		arr[i+1] = p.z.coord;
        arr[i+2] = p.z.acc;
		arr[i+3] = p.z.v;
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
		p.x = {.potential = arr[i], .coord = arr[i+1], .acc = arr[i+2], .v = arr[i+3]};
		i += 4;
        p.y = {.potential = arr[i], .coord = arr[i+1], .acc = arr[i+2], .v = arr[i+3]};
        i += 4;
        p.z = {.potential = arr[i], .coord = arr[i+1], .acc = arr[i+2], .v = arr[i+3]};
        i += 4;
		v.emplace_back(p);
	}
	return v;
}

void printParticles(std::vector<particle3D>& v, std::ofstream& file) {
    for (particle3D &p : v)
        file << p.x.coord << " " << p.y.coord << " " << p.z.coord << " " << p.x.v << " " << p.y.v << " " << p.z.v << "\n";

}