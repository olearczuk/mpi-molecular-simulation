#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>

const double auxE = 4.69041575982343e-08;
const double minDistance = 1e-10;
// TODO 1 + 3 * 4 (after unifying minusPotential and plusPotential)
const int particleSize = 1 + 3 * 5;


struct particle1D {
	// TODO unify plusPotential and minusPotential
	double plusPotential, minusPotential;
	double coord, acc, v;
};

struct particle3D {
	float number;
	particle1D x, y, z;
};

void parseCommandLineArgs(int argc, char *argv[], std::string &inFilename, std::string &outFilename, int &steps,
		double &delta, bool &isVerbose);

std::vector<particle3D> readFile(std::string fileName);

void updateCoords(particle3D &particle, bool delta);
void updatePotential(particle3D &i, particle3D &j, particle3D &k);
void updateAccelerationVelocity(particle3D &particle, bool delta);
void updateAcceleration(particle3D &particle);

void particlesToArray(std::vector<particle3D> v, int begin, int end, double *arr);
std::vector<particle3D> arrayToParticles(double *arr, int size);

#endif //UTILS_H
