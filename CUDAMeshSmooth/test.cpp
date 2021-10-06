//-----------------------------------------------------------
//  Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//  See LICENSE included with this distribution.

#include "stdafx.h"
#include "test.hpp"
#include "mesh_smooth.cuh"

#include <array>
#include <vector>

void testfunction()
{
	int nSweeps = 40;
	size_t nverts = 1500000;
	size_t nidxs = nverts;
	size_t nadj;
	unsigned int nthreadsPerBlock = 128;

	float *vertsx = new float[nverts];
	float *vertsy = new float[nverts];
	float *vertsz = new float[nverts];
	float* normsx = new float[nverts];
	float* normsy = new float[nverts];
	float* normsz = new float[nverts];
	for (size_t i = 0; i < nidxs; ++i) {
		vertsx[i] = float(rand())/float(RAND_MAX);
		vertsy[i] = float(rand())/float(RAND_MAX);
		vertsz[i] = float(rand())/float(RAND_MAX);
		normsx[i] = 1.0f;
		normsy[i] = 0.0f;
		normsz[i] = 0.0f;
	}

	size_t *vertidxs = new size_t[nidxs];
	for(size_t i=0; i < nidxs; ++i) vertidxs[i] = i;

	size_t *vertadjOffsets = new size_t[nidxs+1];
	size_t offset = 0;
	for (size_t i = 0; i < nidxs; ++i) {
		size_t count = 3 + rand() % 8;
//		size_t count = 3 + rand() % 55;
//		size_t count = (i % 2 == 0 ? 3 : 50);
//		size_t count = (i < nidxs/2 ? 3 : 50);
		vertadjOffsets[i] = offset;
		offset += count;
	}	

	vertadjOffsets[nidxs] = offset;
	nadj = offset;

	size_t *vertadj = new size_t[nadj];
	for (size_t i = 0; i < nidxs; ++i) {
		size_t offset = vertadjOffsets[i];
		size_t offsetNxt = vertadjOffsets[i+1];
		for (size_t j = offset; j < offsetNxt; ++j) {
//			vertadj[j] = rand() % nverts;
			vertadj[j] = (i + 20 - rand() % 40) % nverts;
//			vertadj[j] = (i + 10 - rand() % 20) % nverts;
//			vertadj[j] = (i + rand() % count) % nverts;
//			vertadj[j] = (i + j -1) % nidxs;
//			vertadj[j] = (i + 1) % nidxs;
//			vertadj[j] = 0;
		}	
	}

	for (int i = 0; i < 10; ++i) { //***
		std::cout << vertsx[i] << " ";
		std::cout << vertsy[i] << " ";
		std::cout << vertsz[i] << " ";
	}

	bool success;
	cuda_mesh_smooth(
		nSweeps, nverts, nidxs, nadj, nthreadsPerBlock, true,
		vertidxs, vertadj, vertadjOffsets, 
		vertsx, vertsy, vertsz, normsx, normsy, normsz, success);


	for (int i = 0; i < 10; ++i) { //***
		std::cout << vertsx[i] << " ";
		std::cout << vertsy[i] << " ";
		std::cout << vertsz[i] << " ";
	}

	delete[] vertsx;
	delete[] vertsy;
	delete[] vertsz;
	delete[] normsx;
	delete[] normsy;
	delete[] normsz;
	delete[] vertidxs;
	delete[] vertadj;
	delete[] vertadjOffsets;

}
