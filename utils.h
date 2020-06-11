#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>

const double auxE = 4.69041575982343e-08;
const double minDistance = 1e-10;
const int particleSize = 1 + 3 * 4;
const double unitMass = 1;


struct particle1D {
	double potential;
	double coord, acc, v;
};

struct particle3D {
	double number;
	particle1D x, y, z;
};

void parseCommandLineArgs(int argc, char *argv[], std::string &inFilename, std::string &outFilename, int &steps,
		double &delta, bool &isVerbose);

std::vector<particle3D> readFile(const std::string& fileName);

void updateCoords(std::vector<particle3D> &v, double delta);
void updatePotential(std::vector<particle3D> &v1, std::vector<particle3D> &v2, std::vector<particle3D> &v3,
        int owner1, int owner2, int owner3);
void updateAccelerationVelocity(std::vector<particle3D> &v, double delta);
void updateAcceleration(std::vector<particle3D> &v);

void particlesToArray(std::vector<particle3D> v, int begin, int end, double *arr);
std::vector<particle3D> arrayToParticles(const double *arr, int size);

//TODO remove
void printParticle(particle3D p);

#endif //UTILS_H
