//-----------------------------------------------------------
//     Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//     See LICENSE included with this distribution.

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <limits>

const float float_infinity = std::numeric_limits<float>::infinity();
const float float_tiny = std::numeric_limits<float>::min();
const double double_infinity = std::numeric_limits<double>::infinity();
const float double_tiny = std::numeric_limits<double>::min();
const int int_infinity = std::numeric_limits<int>::max();

#endif // CONSTANTS_H
