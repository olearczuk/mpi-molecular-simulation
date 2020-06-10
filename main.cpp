#include <cstdio>
#include <mpi.h>
#include <algorithm>
#include "utils.h"

int beginIndex(int rank, int p, int n) {
	return rank * (n / p) + std::min(rank, n % p);
}

int particlesNumber(int rank, int p, int n) {
	return beginIndex(rank + 1, p, n) - beginIndex(rank, p, n);
}

void propagateParticles(int rank, int p, int n, std::vector<particle3D> v, double *buf) {
	// TODO - maybe it can be done using ISend
	if (rank == 0) {
		double sendBuf[particleSize * (1 + n/p)];
		for (int i = 1; i < p; i++) {
			int begin = beginIndex(i, p, n), nextBegin = beginIndex(i+1, p, n);
			particlesToArray(v, begin, nextBegin, sendBuf);
			MPI_Send(sendBuf, particleSize * particlesNumber(rank, p, n), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
		int begin = beginIndex(rank, p, n), nextBegin = beginIndex(rank+1, p, n);
		particlesToArray(v, begin, nextBegin, buf);
	} else {
		MPI_Recv(buf, particleSize * particlesNumber(rank, p, n), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
				 MPI_STATUS_IGNORE);
	}
}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	int p, n, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	std::string inFilename, outFilename;
	int steps;
	double delta;
	bool isVerbose;
	std::vector<particle3D> v;
	if (rank == 0) {
		parseCommandLineArgs(argc, argv, inFilename, outFilename, steps, delta, isVerbose);
		v = readFile(inFilename);
		n = v.size();
	}
	MPI_Bcast(&delta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&isVerbose, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	double buf0[particleSize * (1 + n / p)];
	propagateParticles(rank, p, n, v, buf0);
	v = arrayToParticles(buf0, particleSize * particlesNumber(rank, p, n));
	printf("%d %f %f %f %f %f %f\n", rank, v[0].x.coord, v[0].y.coord, v[0].z.coord,
		   v[0].x.v, v[0].y.v, v[0].z.v);
	if (particlesNumber(rank, p, n) > 1)
		printf("%d %f %f %f %f %f %f\n", rank, v[1].x.coord, v[1].y.coord, v[1].z.coord,
			   v[1].x.v, v[1].y.v, v[1].z.v);
	MPI_Finalize();
	return 0;
}