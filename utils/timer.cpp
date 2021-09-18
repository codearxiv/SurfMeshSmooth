//-----------------------------------------------------------
//  Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//  See LICENSE included with this distribution.

#include <ctime>
#include <iostream>

//-----------------------------------------------------------

void timer(int startstop, std::clock_t& c_start, std::clock_t& c_end)
{
	if ( startstop == 0 ) {
		c_start = std::clock();
	}
	else {
		c_end = std::clock();
		double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
		std::cout << "\nCPU time used: " << time_elapsed_ms << " ms\n" << std::endl;
	}
}
//-----------------------------------------------------------
