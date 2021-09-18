//-----------------------------------------------------------
//  Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//  See LICENSE included with this distribution.


#ifndef VECT_UTILS_H
#define VECT_UTILS_H

#include <Eigen/Dense>


Eigen::Vector3f point_project_triangle(
		const Eigen::Vector3f& p, const Eigen::Vector3f& a,
		const Eigen::Vector3f& b, const Eigen::Vector3f& c,
		bool& hovering, bool& degenerate, Eigen::Vector3f *bary=nullptr);

Eigen::Vector3f point_project_line(
		const Eigen::Vector3f& p, const Eigen::Vector3f& a,
		const Eigen::Vector3f& b, bool segment, float& t);

Eigen::Vector3f point_project_triangle_nondegen(
		const Eigen::Vector3f& p, const Eigen::Vector3f& a,
		const Eigen::Vector3f& b, const Eigen::Vector3f& c,
		bool& hovering, Eigen::Vector3f *bary=nullptr);

Eigen::Vector3f vector_project_line(
		const Eigen::Vector3f& v, const Eigen::Vector3f& w, bool normalized);

Eigen::Vector3f vector_project_plane(
		const Eigen::Vector3f& v, const Eigen::Vector3f& norm, bool normalized);

Eigen::Vector3f triangle_barycoords(
		const Eigen::Vector3f& p, const Eigen::Vector3f& a,
		const Eigen::Vector3f& b, const Eigen::Vector3f& c);

void inverse_2X2(
		const Eigen::Matrix2f& M, bool &degenerate, Eigen::Matrix2f& Minv);

double point_dist_sq(const Eigen::Vector3f& x, const Eigen::Vector3f& y);

#endif // VECT_UTILS_H
