#pragma once
#ifndef DEVARRAY_CUH
#define DEVARRAY_CUH

#include "stdafx.h"
#include "handle_error.cuh"
#include <cuda.h>

using namespace std;

template <class T>
class DevArray
{
	// public functions
public:
	explicit DevArray()
		: start_(0),
		end_(0)
	{}

	// constructor
	explicit DevArray(size_t size)
	{
		allocate(size);
	}
	// destructor
	~DevArray()
	{
		free();
	}

	// resize the vector
	void resize(size_t size)
	{
		free();
		allocate(size);
	}

	// get the size of the array
	size_t getSize() const
	{
		return end_ - start_;
	}

	// get data
	const T* getData() const
	{
		return start_;
	}

	T* getData()
	{
		return start_;
	}

	// set from host
	void set(const T *src, size_t size)
	{
		size_t min = std::min(size, getSize());
		HANDLE_ERROR(cudaMemcpy(start_, src, min * sizeof(T), cudaMemcpyHostToDevice));
	}
	// get to host
	void get(T *dest, size_t size)
	{
		size_t min = std::min(size, getSize());
		HANDLE_ERROR(cudaMemcpy(dest, start_, min * sizeof(T), cudaMemcpyDeviceToHost));
	}

	// set from device
	void set(const DevArray<T> *src, size_t size)
	{
		size_t min = std::min(size, getSize());
		HANDLE_ERROR(cudaMemcpy(
			start_, src->getData(), min * sizeof(T), cudaMemcpyDeviceToDevice));
	}
	// get to device
	void get(DevArray<T> *dest, size_t size)
	{
		size_t min = std::min(size, getSize());
		HANDLE_ERROR(cudaMemcpy(
			dest->getData(), start_, min * sizeof(T), cudaMemcpyDeviceToDevice));
	}


	// private functions
private:
	// allocate memory on the device
	void allocate(size_t size)
	{
		HANDLE_ERROR(cudaMalloc((void**)&start_, size * sizeof(T)));
		end_ = start_ + size;
	}

	// free memory on the device
	void free()
	{
		if (start_ != 0)
		{
			cudaFree(start_);
			start_ = end_ = 0;
		}
	}

	T* start_;
	T* end_;
};


#endif