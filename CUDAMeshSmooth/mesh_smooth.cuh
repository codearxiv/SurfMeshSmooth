//-----------------------------------------------------------
//  Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//  See LICENSE included with this distribution.

#pragma once
#ifndef MESH_SMOOTH_CUH
#define MESH_SMOOTH_CUH

#ifdef CUDAMESHSMOOTH_EXPORTS
#define CUDAMESHSMOOTH_API __declspec(dllexport)
#else
#define CUDAMESHSMOOTH_API __declspec(dllimport)
#endif

#include "stdafx.h"

extern "C" void CUDAMESHSMOOTH_API cuda_mesh_smooth(
	int nSweeps, size_t nverts, size_t nidxs, size_t nadj,
	unsigned int nthreadsPerBlock, bool localizedVerts,
	const size_t *vertidxs, const size_t *vertadj, 
	const size_t *vertadjOffsets,
	float *vertsx, float *vertsy, float *vertsz, 
	float* normsx, float* normsy, float* normsz, bool& success);

#endif