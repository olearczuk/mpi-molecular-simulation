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

void propagateParticles(int rank, int p, int n, const std::vector<particle3D>& v, double *buf) {
	// TODO - it should be done using MPI_Scatterv
	if (rank == 0) {
		double sendBuf[particleSize * (1 + n/p)];
		for (int i = 1; i < p; i++) {
			int begin = beginIndex(i, p, n), nextBegin = beginIndex(i+1, p, n);
			particlesToArray(v, begin, nextBegin, sendBuf);
			MPI_Send(sendBuf, particleSize * particlesNumber(i, p, n), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
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

	MPI_Request request[6];
	MPI_Status status[6];
	// getting neighbour's particles
	MPI_Isend(buffers[1], particleSize * particlesNumber(rank, p, n), MPI_DOUBLE, prv, 0, MPI_COMM_WORLD,
	        &request[0]);
    MPI_Irecv(buffers[2], particleSize * particlesNumber(nxt, p, n), MPI_DOUBLE, nxt, 0, MPI_COMM_WORLD,
            &request[1]);
    MPI_Waitall(2, request, status);
	MPI_Isend(buffers[1], particleSize * particlesNumber(rank, p, n), MPI_DOUBLE, nxt, 0, MPI_COMM_WORLD,
	        &request[0]);
    MPI_Irecv(buffers[0], particleSize * particlesNumber(prv, p, n), MPI_DOUBLE, prv, 0, MPI_COMM_WORLD,
            &request[1]);
	MPI_Waitall(2, request, status);
	vectors[0] = arrayToParticles(buffers[0], particleSize * particlesNumber(prv, p, n));
	vectors[2] = arrayToParticles(buffers[2], particleSize * particlesNumber(nxt, p, n));
	int i = 0;
	for (int s = p - 3; s >= 0; s -= 3) {
		for (int move = 0; move < s; move++) {
			if (move != 0 || s != p - 3)
				shiftRight(rank, p, n, i, buffers, owners, vectors);
			else {
				updatePotential(vectors[1], vectors[1], vectors[1], owners[1], owners[1], owners[1]);
				updatePotential(vectors[1], vectors[1], vectors[2], owners[1], owners[1], owners[2]);
				updatePotential(vectors[0], vectors[0], vectors[2], owners[0], owners[0], owners[2]);
			}
			if (s == p - 3)
				updatePotential(vectors[0], vectors[1], vectors[1], owners[0], owners[1], owners[1]);
			updatePotential(vectors[0], vectors[1], vectors[2], owners[0], owners[1], owners[2]);
		}
		i = (i + 1) % 3;
	}
	if (p % 3 == 0) {
		i = previous(i, 3);
		shiftRight(rank, p, n, i, buffers, owners, vectors);
		if (rank / (p / 3) == 0)
			updatePotential(vectors[0], vectors[1], vectors[2], owners[0], owners[1], owners[2]);
	}
    particlesToArray(vectors[0], 0, vectors[0].size(), buffers[0]);
    particlesToArray(vectors[1], 0, vectors[1].size(), buffers[1]);
    particlesToArray(vectors[2], 0, vectors[2].size(), buffers[2]);
    MPI_Barrier(MPI_COMM_WORLD);

    // sending stored particles to its owners
    double auxBuffers[3][particleSize * (1 + n / p)];
    MPI_Isend(buffers[0], particleSize * particlesNumber(owners[0], p, n), MPI_DOUBLE, owners[0],
            0, MPI_COMM_WORLD, &request[0]);
    MPI_Irecv(auxBuffers[0], particleSize * particlesNumber(rank, p, n), MPI_DOUBLE, MPI_ANY_SOURCE, 0,
              MPI_COMM_WORLD, &request[1]);
    MPI_Isend(buffers[1], particleSize * particlesNumber(owners[1], p, n), MPI_DOUBLE, owners[1],
            0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(auxBuffers[1], particleSize * particlesNumber(rank, p, n), MPI_DOUBLE, MPI_ANY_SOURCE, 0,
              MPI_COMM_WORLD, &request[3]);
    MPI_Isend(buffers[2], particleSize * particlesNumber(owners[2], p, n), MPI_DOUBLE, owners[2],
            0, MPI_COMM_WORLD, &request[4]);
    MPI_Irecv(auxBuffers[2], particleSize * particlesNumber(rank, p, n), MPI_DOUBLE, MPI_ANY_SOURCE, 0,
              MPI_COMM_WORLD, &request[5]);
    MPI_Waitall(6, request, status);
    // summing potentials of all 3 copies of process's particles
	std::vector<particle3D> vRes[3];
	vRes[0] = arrayToParticles(auxBuffers[0], particleSize * particlesNumber(rank, p, n));
    vRes[1] = arrayToParticles(auxBuffers[1], particleSize * particlesNumber(rank, p, n));
    vRes[2] = arrayToParticles(auxBuffers[2], particleSize * particlesNumber(rank, p, n));
    for (size_t k = 0; k < vRes[0].size(); k++) {
        vRes[0][k].x.potential += vRes[1][k].x.potential + vRes[2][k].x.potential;
        vRes[0][k].y.potential += vRes[1][k].y.potential + vRes[2][k].y.potential;
        vRes[0][k].z.potential += vRes[1][k].z.potential + vRes[2][k].z.potential;
    }
    MPI_Barrier(MPI_COMM_WORLD);
	return vRes[0];
}

// TODO - can be done using MPI_Gatherv
void aggregateParticles(int rank, int p, int n, std::vector<particle3D> &v, double *buffer) {
    if (rank == 0) {
        std::vector<particle3D> vAux;
        for (particle3D &particle : v)
            printParticle(particle);
        for (int i = 1; i < p; i++) {
            MPI_Status status;
            int bufferSize = particleSize * particlesNumber(i, p, n);
            MPI_Recv(buffer, bufferSize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
            vAux = arrayToParticles(buffer, bufferSize);
            for (particle3D &particle : vAux)
                printParticle(particle);
        }
    } else {
        particlesToArray(v, 0, v.size(), buffer);
        MPI_Send(buffer, particleSize * particlesNumber(rank, p, n), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
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

	v = updatePotential(rank, p, n, v, buffers);
    updateAcceleration(v);
	for (; steps > 0; steps--) {
	    // resetting potentials
        for (particle3D &particle : v) {
            particle.x.potential = particle.y.potential = particle.z.potential = 0;
        }
	    updateCoords(v, delta);
	    v = updatePotential(rank, p, n, v, buffers);
	    updateAccelerationVelocity(v, delta);
	    if (isVerbose || steps == 1)
	        aggregateParticles(rank, p, n, v, buffers[0]);
	}
	MPI_Finalize();
	delete [] buffers[0];
	delete [] buffers[1];
	delete [] buffers[2];
	return 0;
}