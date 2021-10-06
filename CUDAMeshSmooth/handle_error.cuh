#pragma once
#ifndef HANDLE_ERROR_CUH
#define HANDLE_ERROR_CUH

#include "stdafx.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

static void HandleError(cudaError_t err, const char* file,int line)
{
	if (err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err),file,line);
		exit(EXIT_FAILURE);
	}
}

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

#endif