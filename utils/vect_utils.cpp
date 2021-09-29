//-----------------------------------------------------------
//  Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//  See LICENSE included with this distribution.

#include "vect_utils.h"
#include "constants.h"

#include <Eigen/Dense>
#include <vector>


using std::max;
using std::min;
using std::abs;
using Eigen::Index;
using Eigen::Vector2f;
using Eigen::Vector3f;
using Eigen::VectorXf;
using Eigen::Matrix2f;
using Eigen::Matrix3f;
using Eigen::MatrixXf;

Vector3f point_project_triangle(
	const Vector3f& p, const Vector3f& a, const Vector3f& b, const Vector3f& c,
	bool& hovering, bool& degenerate, Vector3f *bary)
{
	Vector3f q;

	Vector3f ab = b - a;
	Vector3f ac = c - a;
	Vector3f bc = c - b;

	double lengthsq[3];
	lengthsq[0] = ab.dot(ab);
	lengthsq[1] = ac.dot(ac);
	lengthsq[2] = bc.dot(bc);

	int imin = (lengthsq[0] <= lengthsq[1] ? 0 : 1);
	int imax = (lengthsq[0] >= lengthsq[1] ? 0 : 1);
	imin = (lengthsq[imin] <= lengthsq[2] ? imin : 2);
	imax = (lengthsq[imax] >= lengthsq[2] ? imax : 2);

	degenerate = (lengthsq[imin] < 1e-15*lengthsq[imax]);

	if ( !degenerate ) {
		q = point_project_triangle_nondegen(p, a, b, c, hovering, bary);
	}
	else {
		hovering = false;
		float t;
		switch(imax) {
		case (0):
			q = point_project_line(p, a, b, true, t);
			if ( bary != nullptr ) *bary = Vector3f(1.0f-t, t, 0.0f);
			break;
		case (1):
			q = point_project_line(p, a, c, true, t);
			if ( bary != nullptr ) *bary = Vector3f(1.0f-t, 0.0f, t);
			break;
		case (2):
			q = point_project_line(p, b, c, true, t);
			if ( bary != nullptr ) *bary = Vector3f(0.0f, 1.0f-t, t);
			break;
		}
	}

	return q;
}

//-----------------------------------------------------------

Vector3f point_project_line(
	const Vector3f& p, const Vector3f& a, const Vector3f& b, bool segment, float& t)
{
	Vector3f q;

	Vector3f ab = b - a;
	Vector3f ap = p - a;

	double abab = ab.dot(ab);

	if (abab > 1e-15) {
		t = ab.dot(ap)/abab;

		if ( segment && t <= 0.0f ) {
			q = a;
			t = 0.0f;
		}
		else if ( segment && t >= 0.0f ) {
			q = b;
			t = 1.0f;
		}
		else {
			q = t*ab + a;
		}
	}
	else {
		 q = a;
		 t = 0.0f;
	}

	return q;
}


//-----------------------------------------------------------


Vector3f point_project_triangle_nondegen(
	const Vector3f& p, const Vector3f& a, const Vector3f& b, const Vector3f& c,
	bool& hovering, Vector3f *bary)
{
	Vector3f q;

	Vector3f ab = b - a;
	Vector3f ac = c - a;
	Vector3f bc = c - b;
	Vector3f ap = p - a;
	Vector3f bp = p - b;

	Vector3f n = ab.cross(ac);
	Vector3f na = n.cross(bc);
	Vector3f nb, nc;

	if (ab.dot(na) <= 0.0f) {
		nb = ac.cross(n);
		nc = n.cross(ab);
	}
	else {
		na = -na;
		nb = n.cross(ac);
		nc = ab.cross(n);
	}

	bool sidea = (na.dot(bp) >= 0.0f);
	bool sideb = (nb.dot(ap) >= 0.0f);
	bool sidec = (nc.dot(ap) >= 0.0f);
	hovering = (sidea && sideb && sidec);

	int type = 0;

	if ( hovering ) {
		q = a + vector_project_plane(ap, n, false);
		type = 7;
		goto findbary;
	}

	float t;

	if ( !sidea ) {
		q = point_project_line(p, b, c, true, t);
		if ( t <= 0.0f ) {
			type = 2;
		}
		else if ( t >= 1.0f ) {
			type = 3;
		}
		else {
			type = 6;
			goto findbary;
		}
	}

	if ( !sideb ) {
		q = point_project_line(p, a, c, true, t);
		if ( t <= 0.0f ) {
			type = 1;
		}
		else if ( t >= 1.0f ) {
			type = 3;
		}
		else {
			type = 5;
			goto findbary;
		}
	}

	if ( !sidec ) {
		q = point_project_line(p, a, b, true, t);
		if ( t <= 0.0f ) {
			type = 1;
		}
		else if ( t >= 1.0f ) {
			type = 2;
		}
		else {
			type = 4;
			goto findbary;
		}
	}

	if ( type == 0 ) {
		q = (a + b + c)/3.0f;
	}

findbary:

	if ( bary != nullptr ) {
		switch (type) {
		case (7):\
			*bary = triangle_barycoords(q, a, b, c);
			break;
		case (6):
			*bary = Vector3f(0.0f, 1.0f-t, t);
			break;
		case (5):
			*bary = Vector3f(1.0f-t, 0.0f, t);
			break;
		case (4):
			*bary = Vector3f(1.0f-t, t, 0.0f);
			break;
		case (3):
			*bary = Vector3f(0.0f, 0.0f, 1.0f);
			break;
		case (2):
			*bary = Vector3f(0.0f, 1.0f, 0.0f);
			break;
		case (1):
			*bary = Vector3f(1.0f, 0.0f, 0.0f);
			break;
		default:
			*bary = Vector3f(1.0f/3.0f, 1.0f/3.0f, 1.0f/3.0f);
			break;
		}
	}


	return q;

}


//-----------------------------------------------------------


Vector3f vector_project_line(
		const Vector3f& v, const Vector3f& w, bool normalized)
{

	float scale;
	if ( normalized ) {
		scale = v.dot(w);
	}
	else {
		float ww = w.dot(w);
		if (ww > float_tiny) {
			scale = v.dot(w)/ww;
		}
		else {
			scale = 0.0f;
		}
	}

	return scale*w;
}


//-----------------------------------------------------------

Vector3f vector_project_plane(
		const Vector3f& v, const Vector3f& norm, bool normalized)
{
	Vector3f normproj = vector_project_line(v, norm, normalized);
	return v - normproj;
}


//-----------------------------------------------------------

Vector3f triangle_barycoords(
		const Vector3f& p, const Vector3f& a, const Vector3f& b, const Vector3f& c)
{
	Vector3f u = b - a;
	Vector3f v = c - a;

	Eigen::Matrix<float, 3, 2> M;
	Matrix2f MTM, MTMinv;
	M.col(0) = u;
	M.col(1) = v;
	MTM.noalias() = M.transpose()*M;
	bool degenerate;
	inverse_2X2(MTM, degenerate, MTMinv);

	Vector2f uv;

	if ( !degenerate ) {
		Vector3f w = p - a;
		Vector2f uv0;
		uv0.noalias() = w.transpose()*M;
		uv.noalias() = uv0.transpose()*MTMinv;
	}
	else {
		float t;
		point_project_line(p, a, a+u, false, t);
		uv(0) = t;
		uv(1) = 0.0f;
	}

	Vector3f bary = Vector3f(1.0f - (uv(0)+uv(1)), uv(0), uv(1));
	return bary;
}

//-----------------------------------------------------------

void inverse_2X2(const Matrix2f& M, bool& degenerate, Matrix2f& Minv)
{
	double det = M(0,0)*M(1,1) - M(0,1)*M(1,0);
	if ( abs(det) < double_tiny ) {
		degenerate = true;
	}
	else {
		degenerate = false;
		double r = 1.0f/det;
		Minv(0,0) = r*M(1,1);
		Minv(1,0) = -r*M(1,0);
		Minv(0,1) = -r*M(0,1);
		Minv(1,1) = r*M(0,0);
	}
}

//-----------------------------------------------------------

double point_dist_sq(const Vector3f& x, const Vector3f& y)
{
	Vector3f diff = x - y;
	return diff.dot(diff);
}

//-----------------------------------------------------------
