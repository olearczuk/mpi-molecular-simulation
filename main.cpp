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

int previous(int rank, int p) {
	return (rank + p - 1) % p;
}

int next(int rank, int p) {
	return (rank + 1) % p;
}

void propagateParticles(int rank, int p, int n, std::vector<particle3D> v, double *buf) {
	// TODO - maybe it can be done using Isend
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

void shiftRight(int rank, int p, int n, int i, double **buffers, int *owners, std::vector<particle3D> *vectors) {
	int prv = previous(rank, p), nxt = next(rank, p);
	double auxBuffer[particleSize * (1 + n / p)];
	MPI_Request request[2];
	MPI_Status status[2];
    particlesToArray(vectors[i], 0, vectors[i].size(), buffers[i]);
	MPI_Isend(buffers[i], particleSize * particlesNumber(owners[i], p, n), MPI_DOUBLE, nxt, 0, MPI_COMM_WORLD, &request[0]);
	owners[i] = previous(owners[i], p);
    int auxBufferSize = particleSize * particlesNumber(owners[i], p, n);
	MPI_Irecv(auxBuffer, auxBufferSize, MPI_DOUBLE, prv, 0, MPI_COMM_WORLD, &request[1]);
	MPI_Waitall(2, request, status);
	std::copy(auxBuffer, auxBuffer + auxBufferSize, buffers[i]);
	vectors[i] = arrayToParticles(buffers[i], auxBufferSize);
}

std::vector<particle3D> updatePotential(int rank, int p, int n, const std::vector<particle3D>& v, double **buffers) {
	std::vector<particle3D> vectors[3];
	vectors[1] = v;
	particlesToArray(v, 0, v.size(), buffers[1]);
	int prv = previous(rank, p), nxt = next(rank, p);
	int owners[3] = {prv, rank, nxt};
	MPI_Request request[4];
	MPI_Status status[4];
	MPI_Isend(buffers[1], particleSize * particlesNumber(prv, p, n), MPI_DOUBLE, prv, 0, MPI_COMM_WORLD, &request[0]);
	MPI_Isend(buffers[1], particleSize * particlesNumber(nxt, p, n), MPI_DOUBLE, nxt, 0, MPI_COMM_WORLD, &request[1]);
	MPI_Irecv(buffers[0], particleSize * particlesNumber(prv, p, n), MPI_DOUBLE, prv, 0, MPI_COMM_WORLD, &request[2]);
	MPI_Irecv(buffers[2], particleSize * particlesNumber(nxt, p, n), MPI_DOUBLE, nxt, 0, MPI_COMM_WORLD, &request[3]);
	MPI_Waitall(4, request, status);
	vectors[0] = arrayToParticles(buffers[0], particleSize * particlesNumber(prv, p, n));
	vectors[2] = arrayToParticles(buffers[2], particleSize * particlesNumber(nxt, p, n));
	int i = 0;
	for (int s = p - 3; s >= 0; s -= 3) {
		for (int move = 0; move < s; move++) {
			if (move != 0 || s != p - 3)
				shiftRight(rank, p, n, i, buffers, owners, vectors);
			else {
				updatePotential(vectors[1], vectors[1], vectors[1]);
				updatePotential(vectors[1], vectors[1], vectors[2]);
				updatePotential(vectors[0], vectors[0], vectors[2]);
			}
			if (s == p - 3)
				updatePotential(vectors[0], vectors[1], vectors[1]);
			updatePotential(vectors[0], vectors[1], vectors[2]);
		}
		i = (i + 1) % 3;
	}
	if (p % 3 == 0) {
		i = previous(i, 3);
		shiftRight(rank, p, n, i, buffers, owners, vectors);
		if (rank / (p / 3) == 0)
			updatePotential(vectors[0], vectors[1], vectors[2]);
	}
    particlesToArray(vectors[0], 0, vectors[0].size(), buffers[0]);
    particlesToArray(vectors[1], 0, vectors[1].size(), buffers[1]);
    particlesToArray(vectors[1], 0, vectors[1].size(), buffers[2]);
    MPI_Barrier(MPI_COMM_WORLD);
	std::vector<particle3D> vRes;
	double auxBuffer[particleSize * (1 + n / p)];
	for (int j = 0; j < 3; j++) {
		std::vector<particle3D> vectorAux;
		if (owners[j] != rank) {
			MPI_Isend(buffers[j], particleSize * particlesNumber(owners[j], p, n), MPI_DOUBLE, owners[j], 0,
                      MPI_COMM_WORLD, &request[0]);
			MPI_Irecv(auxBuffer, particleSize * particlesNumber(rank, p, n), MPI_DOUBLE, MPI_ANY_SOURCE, 0,
					MPI_COMM_WORLD, &request[1]);
			MPI_Waitall(2, request, status);
			vectorAux = arrayToParticles(auxBuffer, particleSize * particlesNumber(rank, p, n));
		} else
			vectorAux = arrayToParticles(buffers[j], particleSize * particlesNumber(rank, p, n));
		if (vRes.empty())
			vRes = vectorAux;
		else {
			for (size_t k = 0; k < v.size(); k++) {
				vRes[k].x.minusPotential += vectorAux[k].x.minusPotential;
				vRes[k].x.plusPotential += vectorAux[k].x.plusPotential;

				vRes[k].y.minusPotential += vectorAux[k].y.minusPotential;
				vRes[k].y.plusPotential += vectorAux[k].y.plusPotential;

				vRes[k].z.minusPotential += vectorAux[k].z.minusPotential;
				vRes[k].z.plusPotential += vectorAux[k].z.plusPotential;
			}
		}
	}
    MPI_Barrier(MPI_COMM_WORLD);
	return vRes;
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
	parseCommandLineArgs(argc, argv, inFilename, outFilename, steps, delta, isVerbose);
	if (rank == 0) {
		v = readFile(inFilename);
		n = v.size();
	}
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	double *buffers[3];
	buffers[0] = new double[particleSize * (1 + n / p)];
	buffers[1] = new double[particleSize * (1 + n / p)];
	buffers[2] = new double[particleSize * (1 + n / p)];
	propagateParticles(rank, p, n, v, buffers[1]);
	v = arrayToParticles(buffers[1], particleSize * particlesNumber(rank, p, n));
	MPI_Barrier(MPI_COMM_WORLD);
	// TODO - reset potentials
	v = updatePotential(rank, p, n, v, buffers);
	updateAcceleration(v);
//    printf("%d -- ", rank);
//	printParticle(v[0]);

	for (; steps > 0; steps--) {
	    updateCoords(v, delta);
	    v = updatePotential(rank, p, n, v, buffers);
	    updateAccelerationVelocity(v, delta);
	}
	printParticle(v[0]);

	MPI_Finalize();
	delete [] buffers[0];
	delete [] buffers[1];
	delete [] buffers[2];
	return 0;
}