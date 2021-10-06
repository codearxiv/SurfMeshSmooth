//-----------------------------------------------------------
//  Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//  See LICENSE included with this distribution.

// CUDA laplacian smoothing a surface mesh, adjusted for vertex normals
// to preserve surface curvature.

#include "stdafx.h"
#include "mesh_smooth.cuh"
#include "handle_error.cuh"
#include "devArray.cuh"

#include <array>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>

//---------------------------------------------------------
// Shared memory implementation. Assumes mesh vertices are
// localized into patches, laid out in memory s.t. topologically 
// near vertices tend to be near in memory as well. Otherwise 
// blocking into shared memory won't improve performance. 
// More precisely, assumes that with high probability adjacent 
// vertices v and w satisfy v == verts[i] and w == verts[j] for
// some i and j such that |i-j| < blockDim.x.
//
// ** Note **: Assumes the subset of vertex indices vertidxs
//	to be smoothed is in increasing order and spaced densely, 
//  else again there might not be any performance benefit here.

__global__ void cuda_mesh_smooth_sharedmem(
	size_t nverts, size_t nidxs, 
	const size_t* vertidxs, const size_t *vertadj, 
	const size_t *vertadjOffsets,
	const float *vertsx, const float *vertsy, const float *vertsz,
	const float *normsx, const float *normsy, const float *normsz,
	float* newVertsx, float* newVertsy, float* newVertsz,
	float* newNormsx, float* newNormsy, float* newNormsz)
{
	extern __shared__ char buffer[];
	float3* vertsBlock = (float3*)&buffer[0];
	float3* normsBlock = (float3*)&buffer[blockDim.x*sizeof(float3)];

	const size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	const size_t lower = vertidxs[idx - threadIdx.x];
	const size_t upper = min(nverts, size_t(lower + blockDim.x)) - 1;

	size_t ivert = lower + threadIdx.x;

	if (ivert <= upper) {
		vertsBlock[threadIdx.x].x = vertsx[ivert];
		vertsBlock[threadIdx.x].y = vertsy[ivert];
		vertsBlock[threadIdx.x].z = vertsz[ivert];
		normsBlock[threadIdx.x].x = normsx[ivert];
		normsBlock[threadIdx.x].y = normsy[ivert];
		normsBlock[threadIdx.x].z = normsz[ivert];
	}
	__syncthreads();


	if (idx < nidxs) {
		ivert = vertidxs[idx];

		float3 v0, n0;
		if (ivert < lower || ivert > upper) {
			v0.x = vertsx[ivert];
			v0.y = vertsy[ivert];
			v0.z = vertsz[ivert];
			n0.x = normsx[ivert];
			n0.y = normsy[ivert];
			n0.z = normsz[ivert];
		}
		else {
			v0 = vertsBlock[ivert-lower];
			n0 = normsBlock[ivert-lower];
		}

		float3 v1 = make_float3(0.0f, 0.0f, 0.0f);
		float3 n1 = make_float3(0.0f, 0.0f, 0.0f);

		size_t offset = vertadjOffsets[idx];
		size_t offsetNxt = vertadjOffsets[idx+1];
		size_t size = offsetNxt - offset;
		for (size_t i = offset; i < offsetNxt; ++i) {
			size_t a = vertadj[i];

			float3 v, n;
			if (a < lower || a > upper) {
				v.x = vertsx[a];
				v.y = vertsy[a];
				v.z = vertsz[a];
				n.x = normsx[a];
				n.y = normsy[a];
				n.z = normsz[a];
			}
			else {
				v = vertsBlock[a-lower];
				n = normsBlock[a-lower];
			}

			v.x = v.x - v0.x;
			v.y = v.y - v0.y;
			v.z = v.z - v0.z;
			float vn0 = v.x*n0.x + v.y*n0.y + v.z*n0.z;
			float vn = v.x*n.x + v.y*n.y + v.z*n.z;
			v1.x += 0.5f*v.x + 0.25f*(vn*n.x - vn0*n0.x);
			v1.y += 0.5f*v.y + 0.25f*(vn*n.y - vn0*n0.y);
			v1.z += 0.5f*v.z + 0.25f*(vn*n.z - vn0*n0.z);

			n.x = 0.5f*(n.x + n0.x);
			n.y = 0.5f*(n.y + n0.y);
			n.z = 0.5f*(n.z + n0.z);
			float length = sqrt(n.x*n.x + n.y*n.y + n.z*n.z);	
			n1.x += n.x/length;
			n1.y += n.y/length;
			n1.z += n.z/length;

		}
		newVertsx[ivert] = v1.x/size + v0.x;
		newVertsy[ivert] = v1.y/size + v0.y;
		newVertsz[ivert] = v1.z/size + v0.z;

		float length1 = sqrt(n1.x*n1.x + n1.y*n1.y + n1.z*n1.z);
		newNormsx[ivert] = n1.x/length1;
		newNormsy[ivert] = n1.y/length1;
		newNormsz[ivert] = n1.z/length1;

	}



}


//---------------------------------------------------------
// Naive implementation. Slower in the best case, but without 
// any assumptions attached.

__global__ void cuda_mesh_smooth_naive(
	size_t nidxs, const size_t *vertidxs, 
	const size_t *vertadj, const size_t *vertadjOffsets, 
	const float* vertsx, const float* vertsy, const float* vertsz,
	const float* normsx, const float* normsy, const float* normsz,
	float* newVertsx, float* newVertsy, float* newVertsz,
	float* newNormsx, float* newNormsy, float* newNormsz)
{
	const size_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;

	if (idx < nidxs) {
		size_t ivert = vertidxs[idx];
		float3 v0, n0;
		v0.x = vertsx[ivert];
		v0.y = vertsy[ivert];
		v0.z = vertsz[ivert];
		n0.x = normsx[ivert];
		n0.y = normsy[ivert];
		n0.z = normsz[ivert];

		float3 v1 = make_float3(0.0f, 0.0f, 0.0f);
		float3 n1 = make_float3(0.0f, 0.0f, 0.0f);

		size_t offset = vertadjOffsets[idx];
		size_t offsetNxt = vertadjOffsets[idx+1];
		size_t size = offsetNxt - offset;
		for (size_t i = offset; i < offsetNxt; ++i) {
			size_t a = vertadj[i];
			float3 v, n;
			v.x = vertsx[a];
			v.y = vertsy[a];
			v.z = vertsz[a];
			n.x = normsx[a];
			n.y = normsy[a];
			n.z = normsz[a];

			v.x = v.x - v0.x;
			v.y = v.y - v0.y;
			v.z = v.z - v0.z;
			float vn0 = v.x*n0.x + v.y*n0.y + v.z*n0.z;
			float vn = v.x*n.x + v.y*n.y + v.z*n.z;
			v1.x += 0.5f*v.x + 0.25f*(vn*n.x - vn0*n0.x);
			v1.y += 0.5f*v.y + 0.25f*(vn*n.y - vn0*n0.y);
			v1.z += 0.5f*v.z + 0.25f*(vn*n.z - vn0*n0.z);

			n.x = 0.5f*(n.x + n0.x);
			n.y = 0.5f*(n.y + n0.y);
			n.z = 0.5f*(n.z + n0.z);
			float length = sqrt(n.x*n.x + n.y*n.y + n.z*n.z);
			n1.x += n.x/length;
			n1.y += n.y/length;
			n1.z += n.z/length;
		}
		newVertsx[ivert] = v1.x/size + v0.x;
		newVertsy[ivert] = v1.y/size + v0.y;
		newVertsz[ivert] = v1.z/size + v0.z;

		float length1 = sqrt(n1.x*n1.x + n1.y*n1.y + n1.z*n1.z);
		newNormsx[ivert] = n1.x/length1;
		newNormsy[ivert] = n1.y/length1;
		newNormsz[ivert] = n1.z/length1;

	}

}


//---------------------------------------------------------
// Wrapper calling the kernels

extern "C"  void cuda_mesh_smooth(
	int nSweeps, size_t nverts, size_t nidxs, size_t nadj, 
	unsigned int nthreadsPerBlock, bool localizedVerts,
	const size_t *vertidxs, const size_t *vertadj, 
	const size_t *vertadjOffsets,
	float *vertsx, float *vertsy, float *vertsz, 
	float *normsx, float *normsy, float *normsz, bool& success)
{        
	success = false;
	//cudaEvent_t start, end;
	//float time = 0;

	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	if ( deviceCount == 0 ) return;

	if (nSweeps <= 0 || nverts <= 0 || nidxs <= 0 || nadj <= 0) {
		success = true;
		return;
	}

	DevArray<size_t> d_vertidxs(nidxs);
	DevArray<size_t> d_vertadj(nadj);
	DevArray<size_t> d_vertadjOffsets(nidxs+1);
	DevArray<float> d_vertsx(nverts);
	DevArray<float> d_vertsy(nverts);
	DevArray<float> d_vertsz(nverts);
	DevArray<float> d_normsx(nverts);
	DevArray<float> d_normsy(nverts);
	DevArray<float> d_normsz(nverts);
	DevArray<float> d_newVertsx(nverts);
	DevArray<float> d_newVertsy(nverts);
	DevArray<float> d_newVertsz(nverts);
	DevArray<float> d_newNormsx(nverts);
	DevArray<float> d_newNormsy(nverts);
	DevArray<float> d_newNormsz(nverts);

	d_vertidxs.set(&vertidxs[0], nidxs);
	d_vertadj.set(&vertadj[0], nadj);
	d_vertadjOffsets.set(&vertadjOffsets[0], nidxs+1);
	d_vertsx.set(&vertsx[0], nverts);
	d_vertsy.set(&vertsy[0], nverts);
	d_vertsz.set(&vertsz[0], nverts);
	d_normsx.set(&normsx[0], nverts);
	d_normsy.set(&normsy[0], nverts);
	d_normsz.set(&normsz[0], nverts);

	d_newVertsx.set(&d_vertsx, nverts);
	d_newVertsy.set(&d_vertsy, nverts);
	d_newVertsz.set(&d_vertsz, nverts);
	d_newNormsx.set(&d_normsx, nverts);
	d_newNormsy.set(&d_normsy, nverts);
	d_newNormsz.set(&d_normsz, nverts);


	DevArray<float> *pd_vertsx = &d_vertsx;
	DevArray<float> *pd_vertsy = &d_vertsy;
	DevArray<float> *pd_vertsz = &d_vertsz;
	DevArray<float> *pd_normsx = &d_normsx;
	DevArray<float> *pd_normsy = &d_normsy;
	DevArray<float> *pd_normsz = &d_normsz;
	DevArray<float> *pd_newVertsx = &d_newVertsx;
	DevArray<float> *pd_newVertsy = &d_newVertsy;
	DevArray<float> *pd_newVertsz = &d_newVertsz;
	DevArray<float> *pd_newNormsx = &d_newNormsx;
	DevArray<float> *pd_newNormsy = &d_newNormsy;
	DevArray<float> *pd_newNormsz = &d_newNormsz;

	auto swap_vert_norm_buffers = [&]() ->void {
		std::swap(pd_vertsx, pd_newVertsx);
		std::swap(pd_vertsy, pd_newVertsy);
		std::swap(pd_vertsz, pd_newVertsz);
		std::swap(pd_normsx, pd_newNormsx);
		std::swap(pd_normsy, pd_newNormsy);
		std::swap(pd_normsz, pd_newNormsz);
	};

	size_t nblocks = (nidxs + nthreadsPerBlock - 1)/nthreadsPerBlock;
	unsigned int nshared = 2*nthreadsPerBlock*sizeof(float3);

	cudaFuncSetCacheConfig(cuda_mesh_smooth_naive, cudaFuncCachePreferL1);
//	cudaFuncSetCacheConfig(cuda_mesh_smooth_sharedmem, cudaFuncCachePreferShared);
	cudaFuncSetCacheConfig(cuda_mesh_smooth_sharedmem, cudaFuncCachePreferEqual);

	//cudaEventCreate(&start);
	//cudaEventCreate(&end);
	//cudaEventRecord(start);


	if (localizedVerts) {
		// Vertices ordered into patches on the surface
		for (int i = 0; i < nSweeps; ++i) {
			cuda_mesh_smooth_sharedmem<<<nblocks, nthreadsPerBlock, nshared>>>(
				nverts, nidxs, d_vertidxs.getData(), d_vertadj.getData(),
				d_vertadjOffsets.getData(),
				pd_vertsx->getData(), pd_vertsy->getData(), pd_vertsz->getData(),
				pd_normsx->getData(), pd_normsy->getData(), pd_normsz->getData(),
				pd_newVertsx->getData(), pd_newVertsy->getData(), pd_newVertsz->getData(),
				pd_newNormsx->getData(), pd_newNormsy->getData(), pd_newNormsz->getData());
			//Swap buffer pointers, new coordinates becoming current
			if (i < nSweeps-1) swap_vert_norm_buffers();
		}
	}
	else {
		for (int i = 0; i < nSweeps; ++i) {
			cuda_mesh_smooth_naive<<<nblocks, nthreadsPerBlock>>>(
				nidxs, d_vertidxs.getData(), d_vertadj.getData(),
				d_vertadjOffsets.getData(),
				pd_vertsx->getData(), pd_vertsy->getData(), pd_vertsz->getData(),
				pd_normsx->getData(), pd_normsy->getData(), pd_normsz->getData(),
				pd_newVertsx->getData(), pd_newVertsy->getData(), pd_newVertsz->getData(),
				pd_newNormsx->getData(), pd_newNormsy->getData(), pd_newNormsz->getData());
			//Swap buffer pointers, new coordinates becoming current
			if (i < nSweeps-1) swap_vert_norm_buffers();  
		}
	}

	cudaDeviceSynchronize();
	//cudaEventRecord(end);
	//cudaEventSynchronize(end);
	//cudaEventElapsedTime(&time, start, end);
	//cout << "Execution Time: " << time << endl;

	pd_newVertsx->get(&vertsx[0], nverts);
	pd_newVertsy->get(&vertsy[0], nverts);
	pd_newVertsz->get(&vertsz[0], nverts);
	pd_newNormsx->get(&normsx[0], nverts);
	pd_newNormsy->get(&normsy[0], nverts);
	pd_newNormsz->get(&normsz[0], nverts);

	success = true;

}
