#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>

const double auxE = 4.69041575982343e-08;
const double unitMass = 1.66053906660e-27;
const double minDistance = 1e-10;


struct particle1D {
	double plusPotential, minusPotential;
	double coord, acc, v;
};

struct particle3D {
	particle1D x, y, z;
};

void parseCommandLineArgs(int argc, char *argv[], std::string &inFilename, std::string &outFilename, int &steps,
		double &delta, bool &isVerbose);

std::vector<particle3D> readFile(std::string fileName);

void updateCoords(particle3D &particle, bool delta);
void updatePotential(particle3D &i, particle3D &j, particle3D &k);
void updateAccelerationVelocity(particle3D &particle, bool delta);

#endif //UTILS_H
