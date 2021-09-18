//-----------------------------------------------------------
//  Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//  See LICENSE included with this distribution.

#ifndef SELECTION_SORT_H
#define SELECTION_SORT_H


#include <vector>

namespace std {
template<class T> class function;
}


template <typename T>
void selection_sort(
		std::vector<T>& arr, std::vector<T>& sorted,
		size_t start, size_t end,
		const std::function<bool(T&, T&)> compare)
{
	size_t size = end-start+1;
	sorted.resize(size);

	for(size_t i=0; i<size; ++i) { sorted[i] = arr[start+i]; }

	size_t jmin;
	for(size_t i=0; i<size; ++i) {
		jmin = i;
		for(size_t j=i+1; j<size; ++j) {
			if ( compare(sorted[j],sorted[jmin]) ) jmin = j;
		}
		std::swap(sorted[i], sorted[jmin]);
	}
}



#endif // SELECTION_SORT_H
