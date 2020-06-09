#include "utils.h"
#include <cstdio>

int main(int argc, char *argv[]) {
	std::string inFilename, outFilename;
	int steps;
	double delta;
	bool isVerbose;
	parseCommandLineArgs(argc, argv, inFilename, outFilename, steps, delta, isVerbose);
	printf("%s %s %d %f %d\n", inFilename.c_str(), outFilename.c_str(), steps, delta, isVerbose);
	auto v = readFile(inFilename);
	printf("--------------\n");
	for (particle3D &p : v) {
		printf("%f %f %f %f %f %f\n", p.x.coord, p.y.coord, p.z.coord,
				p.x.v, p.y.v, p.z.v);
	}
	return 0;
}