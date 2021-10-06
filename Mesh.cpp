//-----------------------------------------------------------
//  Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//  See LICENSE included with this distribution.


#include "Mesh.h"
#include "constants.h"
#include "MessageLogger.h"
#include "rotations.h"
#include "Plane.h"
#include "alignment.h"
#include "PairHash.h"
#include "timer.h"
#include "selection_sort.h"
#include "vect_utils.h"
//#include "mesh_smooth.cuh"


#include <Eigen/Core>
#include <qopengl.h>
#include <QCoreApplication>
#include <QRecursiveMutex>
#include <omp.h>
#include <cstdlib>
#include <numeric>
#include <random>
#include <iterator>
#include <algorithm>
#include <functional>
#include <iostream>
#include <ctime>
#include <unordered_map>
#include <utility>

using std::max;
using std::min;
using std::abs;
using std::floor;
using std::pair;
using std::tuple;
using Eigen::Index;
using Eigen::Vector2f;
using Eigen::Vector3f;
using Eigen::VectorXf;
using Eigen::Matrix3f;
using Eigen::MatrixXf;
using Eigen::MatrixXi;
using MapMtrxf = Eigen::Map<MatrixXf, ALIGNEDX>;
using MapMtrxi = Eigen::Map<MatrixXi, ALIGNEDX>;
template<typename T> using aligned = Eigen::aligned_allocator<T>;
template<typename T> using vector = std::vector<T>;
using vectorfa = std::vector<float, aligned<float>>;
using vectoria = std::vector<int, aligned<int>>;

//---------------------------------------------------------

Mesh::Mesh(MessageLogger* msgLogger, QObject *parent)
	: QObject(parent), m_msgLogger(msgLogger)

{
	m_verts.reserve(2500);
	m_norms.reserve(2500);
	m_vertGL.reserve(2500 * 6);

	m_mesh_smooth_GPU = nullptr;
	m_localizedVerts = false;
	m_staleGPUGeomData = true;
	m_staleGPUTopoData = true;
}

//---------------------------------------------------------

Mesh::~Mesh()
{
	// Don't call any function locked on m_reMutex here,
	// else meshWorker thread may wait on the lock owned
	// by another running process before being killed.

}


//---------------------------------------------------------

void Mesh::clear()
{

	QMutexLocker locker(&m_recMutex);

	m_verts.clear();
	m_norms.clear();
	m_tris.clear();
//	m_tets.clear();
	m_boundVert.clear();
	m_triadj.clear();
	m_vertadj.clear();
	m_vertadjOffsets.clear();
	m_vertTri.clear();

	m_vertsx_GPU.clear();
	m_vertsy_GPU.clear();
	m_vertsz_GPU.clear();
	m_normsx_GPU.clear();
	m_normsy_GPU.clear();
	m_normsz_GPU.clear();
	m_vertidxs_GPU.clear();
	m_vertadj_GPU.clear();
	m_vertadjOffsets_GPU.clear();
	m_localizedVerts = false;
	m_staleGPUGeomData = true;
	m_staleGPUTopoData = true;

	m_debug.clear();
	m_vertGL.clear();
	m_triGL.clear();
	m_normGL.clear();
	m_debugGL.clear();

}
//---------------------------------------------------------


const GLfloat* Mesh::vertGLData()
{
	m_vertGL.resize(6 * m_verts.size());

	for (Index i=0, j=0; i < m_verts.size(); ++i){
		m_vertGL[j] = m_verts[i][0];
		m_vertGL[j+1] = m_verts[i][1];
		m_vertGL[j+2] = m_verts[i][2];
		m_vertGL[j+3] = m_norms[i][0];
		m_vertGL[j+4] = m_norms[i][1];
		m_vertGL[j+5] = m_norms[i][2];
		j += 6;
	}

	return static_cast<const GLfloat*>(m_vertGL.data());
}


//---------------------------------------------------------


const GLuint* Mesh::triGLData()
{
	m_triGL.resize(3 * m_tris.size());

	for (Index i=0, j=0; i < m_tris.size(); ++i){
		m_triGL[j] = m_tris[i][0];
		m_triGL[j+1] = m_tris[i][1];
		m_triGL[j+2] = m_tris[i][2];
		j += 3;
	}

	return static_cast<const GLuint*>(m_triGL.data());
}

//---------------------------------------------------------

const GLfloat* Mesh::normGLData(float scale)
{
	m_normGL.resize(12 * m_verts.size());

	for (Index i=0, j=0; i < m_verts.size(); ++i){
		m_normGL[j] = m_verts[i][0] - scale*m_norms[i][0];
		m_normGL[j+1] = m_verts[i][1] - scale*m_norms[i][1];
		m_normGL[j+2] = m_verts[i][2] - scale*m_norms[i][2];
		m_normGL[j+3] = m_norms[i][0];
		m_normGL[j+4] = m_norms[i][1];
		m_normGL[j+5] = m_norms[i][2];
		m_normGL[j+6] = m_verts[i][0] + scale*m_norms[i][0];
		m_normGL[j+7] = m_verts[i][1] + scale*m_norms[i][1];
		m_normGL[j+8] = m_verts[i][2] + scale*m_norms[i][2];
		m_normGL[j+9] = m_norms[i][0];
		m_normGL[j+10] = m_norms[i][1];
		m_normGL[j+11] = m_norms[i][2];
		j += 12;
	}

	return static_cast<const GLfloat*>(m_normGL.data());
}

//---------------------------------------------------------

const GLfloat* Mesh::debugGLData()
{
	Vector3f norm(0.0f, 0.0f, 1.0f);

	m_debugGL.resize(12 * m_debug.size());

	for (Index i=0, j=0; i < m_debug.size(); ++i){
		Vector3f p = m_debug[i].first;
		Vector3f q = m_debug[i].second;
		m_debugGL[j] = p(0);
		m_debugGL[j+1] = p(1);
		m_debugGL[j+2] = p(2);
		m_debugGL[j+3] = norm(0);
		m_debugGL[j+4] = norm(1);
		m_debugGL[j+5] = norm(2);
		m_debugGL[j+6] = q(0);
		m_debugGL[j+7] = q(1);
		m_debugGL[j+8] = q(2);
		m_debugGL[j+9] = norm(0);
		m_debugGL[j+10] = norm(1);
		m_debugGL[j+11] = norm(2);
		j += 12;
	}

	return static_cast<const GLfloat*>(m_debugGL.data());
}


//---------------------------------------------------------

Index Mesh::addVert(const Vector3f& v, const Vector3f& n, bool threadSafe)
{

	if ( threadSafe ) m_recMutex.lock();

	Index idx = m_verts.size();
	m_verts.push_back(v);
	m_norms.push_back(n);

	if ( threadSafe ) m_recMutex.unlock();

	return idx;
}


//---------------------------------------------------------

Index Mesh::addTri(const Vector3i& iv, bool threadSafe)
{

	if ( threadSafe ) m_recMutex.lock();

	Index idx = m_tris.size();
	m_tris.push_back(iv);

	if ( threadSafe ) m_recMutex.unlock();

	return idx;
}
//---------------------------------------------------------

void Mesh::fromPlaneHeight(
		float normx, float normy, float normz, size_t ndim,
		const std::function<float(float xu, float xv)> heightFun)
{
	if (ndim <= 1) return;

	QMutexLocker locker(&m_recMutex);

	m_msgLogger->logMessage("Generating random mesh...");

	clear();

	Vector3f norm(normx, normy, normz);
	double step = 3.0/ndim;
	Plane plane(Vector3f(0.0f,0.0f,0.0f), norm);
	Vector3f u, v;
	plane.getUVAxes(u,v);

	for (Index i = 0; i < ndim; ++i){
		double ustep = i*step - 1.5;
		Vector3f qu = ustep*u;
		for (Index j = 0; j < ndim; ++j){
			double vstep = j*step - 1.5;
			Vector3f q = qu + vstep*v;
			if ( heightFun != nullptr ){
				float normScale = heightFun(ustep, vstep);
				q = q + normScale*norm;
			}
			addVert(q, norm);
		}
	}

	if (m_msgLogger != nullptr) {
		m_msgLogger->logMessage(QString::number(vertCount()) + " vertices created.");
	}

	for (Index i = 0; i < ndim-1; ++i){
		Vector3i iv, iw;
		Index k = i*ndim;
		for (Index j = 0; j < ndim-1; ++j){
			int l = k+j;
			iv(0) = l;
			iv(1) = l+1;
			iv(2) = l+ndim;
			iw(0) = iv(2);
			iw(1) = iv(1);
			iw(2) = l+ndim+1;
			addTri(iv);
			addTri(iw);
		}
	}

	m_localizedVerts = false;
	m_staleGPUGeomData = true;
	m_staleGPUTopoData = true;

	if (m_msgLogger != nullptr) {
		m_msgLogger->logMessage(QString::number(triCount()) + " triangles created.");
	}


	buildVertAdjacency();
	// Reorder vertices into patches to speed up CUDA smoothing.
	localizeVertIndices(128);
	buildTriAdjacency();
	buildVertToTriMap();
	approxMeshNorms();

//	int number;
//	testfunction(number);
//	m_msgLogger->logMessage("number: " + QString::number(int(number)));
}
//---------------------------------------------------------

//void Mesh::testNearestTri()
//{
//	for (int j=0; j<100; ++j){
//		Vector3f x = Vector3f::Random();

//		Index ntris = m_tris.size();
//		Index itriNear = -1;
//		Vector3f projNear;
//		double dmin = float_infinity;
//		bool hovering, degenerate;

//		for (Index itri=0; itri<ntris; ++itri) {
//			Vector3f va = m_verts[ m_tris[itri][0] ];
//			Vector3f vb = m_verts[ m_tris[itri][1] ];
//			Vector3f vc = m_verts[ m_tris[itri][2] ];
//			Vector3f proj = point_project_triangle(
//						x, va, vb, vc, hovering, degenerate);
//			double distsq = point_dist_sq(proj,x);
//			if ( distsq < dmin ) {
//				dmin = distsq;
//				projNear = proj;
//				itriNear = itri;
//			}
//		}
//		m_debug.push_back(pair<Vector3f,Vector3f>(x, projNear));
//		Vector3f va = m_verts[ m_tris[itriNear][0] ];
//		Vector3f vb = m_verts[ m_tris[itriNear][1] ];
//		Vector3f vc = m_verts[ m_tris[itriNear][2] ];
//		Vector3f bary = triangle_barycoords(projNear, va, vb, vc);
//		Vector3f projNear2 = bary(0)*va + bary(1)*vb + bary(2)*vc;
//		m_debug.push_back(pair<Vector3f,Vector3f>(x, projNear2));
//	}
//}

//---------------------------------------------------------

void Mesh::approxMeshNorms()
{
	QMutexLocker locker(&m_recMutex);

	m_msgLogger->logMessage("Building vertex normals...");

	Index nverts = m_verts.size();
	Index ntris = m_tris.size();

//	size_t threshold = 0, lastPos = 0, progress = 0;
	std::clock_t c_start, c_end;//!***
	timer(0, c_start, c_end);

	const bool useOpenMP = true;
	bool vertwiseNotTriwise = true;//(m_vertTri.size() > 0) && (m_triadj.size() > 0);

	if ( vertwiseNotTriwise ) {
		// No thread contention if we can compute normals per-vertex
		// instead of per-triangle and average.
#pragma omp parallel if (useOpenMP) default(shared)
{
		auto tri_vert_side = [](Vector3i *triVerts, Index ivert) -> int {
			int iside;
			for(iside=0; iside<=2; ++iside) {
				if ((*triVerts)[iside] == ivert) break;
			}
			return iside;
		};

#pragma omp for schedule(static)
		for (Index ivert=0; ivert < nverts; ++ivert) {
			Index itri0 = m_vertTri[ivert];
			int iside0 = tri_vert_side(&m_tris[itri0], ivert);
			Index itrib0 = m_triadj[itri0][iside0!=2 ? iside0+1 : 0];
			Index itric0 = m_triadj[itri0][iside0!=0 ? iside0-1 : 2];
			Vector3f normAvg(0.0f,0.0f,0.0f);
			Vector3f dvert = m_verts[ivert];
			// Walk around ivert computing and averaging adjacent
			// triangle normals as we go.
			for (int idir=0; idir <= 1; ++idir) {
				Index itriPrev = (idir == 0 ? itrib0 : itri0);
				Index itri = (idir == 0 ? itri0 : itric0);

				while (itri >= 0) {
					Vector3i *triVerts = &m_tris[itri];
					int iside = tri_vert_side(triVerts, ivert);
					int isideb = (iside != 2 ? iside+1 : 0);
					int isidec = (iside != 0 ? iside-1 : 2);
					Index b = (*triVerts)[isideb];
					Index c = (*triVerts)[isidec];
					Vector3f vab = (m_verts[b] - dvert);
					Vector3f vac = (m_verts[c] - dvert);
					Vector3f vx = vab.cross(vac);
					float x = vx.dot(vx);
					if ( x <= float_tiny ) continue;
					float ab = vab.dot(vab);
					float ac = vac.dot(vac);
					float wa = 1.0f - vab.dot(vac)/sqrt(ab*ac);
					normAvg += wa*vx;  // weighted average by triangle shape

					Index itrib = m_triadj[itri][iside];
					Index itric = m_triadj[itri][isidec];
					Index itriNext = (itrib != itriPrev ? itrib : itric);
					itriPrev = itri;
					itri = itriNext;
					if ( itri == itri0 ) break;
				}
				// Done if we went all the way around ivert.
				// Else hit the boundary, go opposite direction.
				if ( itri == itri0 ) break;
			}
			m_norms[ivert] = normAvg.normalized();
		}
} //Parallel

	}
	else {
		m_norms.assign(nverts, Vector3f(0.0f,0.0f,0.0f));

		for (Index itri=0; itri < ntris; ++itri) {
			Vector3i *triVerts = &m_tris[itri];
			Index a = (*triVerts)[0];
			Index b = (*triVerts)[1];
			Index c = (*triVerts)[2];
			Vector3f vab = (m_verts[b] - m_verts[a]);
			Vector3f vac = (m_verts[c] - m_verts[a]);
			Vector3f vbc = (m_verts[c] - m_verts[b]);
			Vector3f vx = vab.cross(vac);
			float x = vx.dot(vx);
			if ( x <= float_tiny ) continue;
			float ab = vab.dot(vab);
			float ac = vac.dot(vac);
			float bc = vbc.dot(vbc);
			float wa = 1.0f - vab.dot(vac)/sqrt(ab*ac);
			float wb = 1.0f + vab.dot(vbc)/sqrt(ab*bc); //reflect vab
			float wc = 1.0f - vac.dot(vbc)/sqrt(ac*bc);
			m_norms[a] += wa*vx;
			m_norms[b] += wb*vx;
			m_norms[c] += wc*vx;
			//// Log progress
			//if (m_msgLogger != nullptr) {
			//	++progress;
			//	m_msgLogger->logProgress(
			//				"Building vertex normals",
			//				progress, ntris, 1, threshold, lastPos);
			//}

		}

#pragma omp parallel if (useOpenMP) default(shared)
{
#pragma omp for schedule(static)
		for (Index ivert=0; ivert < nverts; ++ivert) {
			m_norms[ivert].normalize();
		}
} //Parallel

	}


	timer(1, c_start, c_end);


//	m_msgLogger->logMessage("Building vertex normals: 100%...", 0);
	m_msgLogger->logMessage("Done!",1);


}


////---------------------------------------------------------
//// SLOW, use buildTriAdjacency

//void Mesh::buildTriAdjacency2()
//{
//	QMutexLocker locker(&m_recMutex);

//	m_msgLogger->logMessage("Building triangle adjacency...");

//	Index ntris = m_tris.size();
//	Index nverts = m_verts.size();

//	typedef pair<Index, Index> TriEdge;
//	typedef pair<Index, int> TriSide;
//	typedef pair<TriSide, TriSide> TriAdj;

////	std::clock_t c_start, c_end;//!***
////	timer(0, c_start, c_end); //!***

//	// get a map of edges to containing triangles
//	std::unordered_map<TriEdge, TriAdj, PairHash> edgeadj;
//	//edgeadj.max_load_factor(1.0f);

//	for (Index itri=0; itri < ntris; ++itri) {
//		for (int i=0; i <= 2; ++i) {
//			int j = (i == 2 ? 0 : i+1);
//			Index a = m_tris[itri][i];
//			Index b = m_tris[itri][j];
//			TriEdge edge(min(a,b), max(a,b));
//			try {
//				TriAdj *padj = &edgeadj.at(edge);
//				if ( (*padj).first.second == -1 ) {
//					(*padj).first = TriSide(itri,i);
//				}
//				else {
//					(*padj).second = TriSide(itri,i);
//				}

//			}
//			catch(const std::out_of_range &e) {
//				edgeadj[edge] = TriAdj(TriSide(itri,i), TriSide(0,-1));
//			}
//		}
//	}


//	// fill in the adjacency array over each triangle
//	m_triadj.assign(ntris, Vector3i(-1,-1,-1));

//	for ( const auto& key : edgeadj ) {
//		TriAdj adj = key.second;
//		int ia = adj.first.second;
//		int ib = adj.second.second;
//		if ( ia == -1 || ib == -1 ) continue;
//		Index itria = adj.first.first;
//		Index itrib = adj.second.first;
//		m_triadj[itria][ia] = itrib;
//		m_triadj[itrib][ib] = itria;
//////debug
////		Vector3f u = Vector3f(1e-3,1e-3,1e-3);
////		Vector3f centa =
////				m_verts[m_tris[itria][0]] +
////				m_verts[m_tris[itria][1]] +
////				m_verts[m_tris[itria][2]];
////		centa = centa/3.0f;
////		Vector3f centb =
////				m_verts[m_tris[itrib][0]] +
////				m_verts[m_tris[itrib][1]] +
////				m_verts[m_tris[itrib][2]];
////		centb = centb/3.0f;
////		m_debug.push_back(pair<Vector3f,Vector3f>(
////							  centa, centb));
//	}

//	m_boundVert.assign(nverts, false);

//	// flag boundary vertices
//	for (Index itri=0; itri < ntris; ++itri) {
//		for (int i=0; i <= 2; ++i) {
//			if ( m_triadj[itri][i] > -1 ) continue;
//			int j = (i == 2 ? 0 : i+1);
//			Index a = m_tris[itri][i];
//			Index b = m_tris[itri][j];
//			m_boundVert[a] = true;
//			m_boundVert[b] = true;
//////debug
////			Vector3f u = Vector3f(1e-2,1e-2,1e-2);
////			m_debug.push_back(pair<Vector3f,Vector3f>(
////								  m_verts[a], m_verts[a]+u));
////			m_debug.push_back(pair<Vector3f,Vector3f>(
////								  m_verts[b], m_verts[b]+u));

//		}
//	}

////	timer(1, c_start, c_end); //!***

//	m_msgLogger->logMessage("Done!",1);
//}

//---------------------------------------------------------
// Same as buildTriAdjacency2 but avoids hashing slowness
// Build mapping triangles -> 3 adjacent triangles sharing
// an edge.

void Mesh::buildTriAdjacency()
{
	QMutexLocker locker(&m_recMutex);

	m_msgLogger->logMessage("Building triangle adjacency...");

	using std::get;
	Index ntris = m_tris.size();
	Index nverts = m_verts.size();
//	std::clock_t c_start, c_end;//!***
//	timer(0, c_start, c_end); //!***

	// Find the number of vertices adjacent to each vertex.
	// This is the size of each bin to store them.
	vector<Index> binSizes(nverts, 0);

	for (Index itri=0; itri < ntris; ++itri) {
		Vector3i *triVerts = &m_tris[itri];
		Index a = (*triVerts)[0];
		Index b = (*triVerts)[1];
		Index c = (*triVerts)[2];
		binSizes[a]+=2;
		binSizes[b]+=2;
		binSizes[c]+=2;
	}

	// Find total memory to allocate for bins
	vector<Index> binOffsets(nverts);
	Index sizeAccum = 0;
	Index sizeMax = 0;

	for (Index ivert=0; ivert < nverts; ++ivert) {
		binOffsets[ivert] = sizeAccum;
		sizeAccum += binSizes[ivert];
		sizeMax = max(sizeMax, binSizes[ivert]);
	}

	// Populate the bins with vertex adjacency data
	typedef std::tuple<Index, Index, uint> TriVert;
	vector<Index> binCounts(nverts, 0);
	vector<TriVert> bins(sizeAccum);

	for (Index itri=0; itri < ntris; ++itri) {
		for (int i=0; i <= 2; ++i) {
			int j = (i == 2 ? 0 : i+1);
			Index a = m_tris[itri][i];
			Index b = m_tris[itri][j];
			bins[binOffsets[a]+binCounts[a]] = TriVert(b,itri,i);
			bins[binOffsets[b]+binCounts[b]] = TriVert(a,itri,i);
			++binCounts[a];
			++binCounts[b];
		}
	}

	// Sort the vertices in each bin to find duplicates.
	// Duplicates are then adjacent triangles.
	auto compare = [](TriVert v, TriVert w) ->bool {
		return get<0>(v) < get<0>(w);
	};

	m_triadj.assign(ntris, Vector3i(-1,-1,-1));
	vector<TriVert> sorted(sizeMax);

	for (Index ivert=0; ivert < nverts; ++ivert) {
		Index offset = binOffsets[ivert];
		Index size = binSizes[ivert];
		// Sort bin to find duplicates
		if ( size < 100 ) {
			// Quicker to do O(n^2) selection sort for small arrays,
			// also no temporary array allocation under the hood.
			selection_sort<TriVert>(bins, sorted, offset, offset+size-1, compare);
		}
		else{
			sorted.assign(bins.begin()+offset, bins.begin()+offset+size);
			std::sort(sorted.begin(), sorted.end(), compare);
		}
		for (Index i=1; i < size; ++i) {
			TriVert v = sorted[i-1];
			TriVert w = sorted[i];			
			if ( get<0>(v) == get<0>(w) ) {
				Index triv = get<1>(v);
				Index triw = get<1>(w);
				uint sidev = get<2>(v);
				uint sidew = get<2>(w);
				m_triadj[triv][sidev] = triw;
				m_triadj[triw][sidew] = triv;
//debug
//				m_debug.push_back(pair<Vector3f,Vector3f>(
//							m_verts[m_tris[triv][sidev]],
//							m_verts[m_tris[triw][sidew]]));
//				Vector3f centa =
//						m_verts[m_tris[triv][0]] +
//						m_verts[m_tris[triv][1]] +
//						m_verts[m_tris[triv][2]];
//				centa = centa/3.0f;
//				Vector3f centb =
//						m_verts[m_tris[triw][0]] +
//						m_verts[m_tris[triw][1]] +
//						m_verts[m_tris[triw][2]];
//				centb = centb/3.0f;
//				m_debug.push_back(pair<Vector3f,Vector3f>(
//									  centa, centb));

			}
		}
	}

	if ( m_boundVert.size() == 0 ) {
		// flag boundary vertices
		m_boundVert.assign(nverts, false);
		for (Index itri=0; itri < ntris; ++itri) {
			for (int i=0; i <= 2; ++i) {
				if ( m_triadj[itri][i] > -1 ) continue;
				int j = (i == 2 ? 0 : i+1);
				Index a = m_tris[itri][i];
				Index b = m_tris[itri][j];
				m_boundVert[a] = true;
				m_boundVert[b] = true;
				////debug
				//Vector3f u = Vector3f(1e-2,1e-2,1e-2);
				//m_debug.push_back(pair<Vector3f,Vector3f>(
				//					  m_verts[a], m_verts[a]+u));
				//m_debug.push_back(pair<Vector3f,Vector3f>(
				//					  m_verts[b], m_verts[b]+u));

			}
		}
	}

//	timer(1, c_start, c_end); //!***

	m_msgLogger->logMessage("Done!",1);
}



//---------------------------------------------------------
// Build mapping vertices -> adjacent vertices on triangle
// edges.

void Mesh::buildVertAdjacency()
{
	QMutexLocker locker(&m_recMutex);

	m_msgLogger->logMessage("Building vertex adjacency...");

	Index ntris = m_tris.size();
	Index nverts = m_verts.size();

//	std::clock_t c_start, c_end;//!***
//	timer(0, c_start, c_end); //!***

	// Find the number of vertices adjacent to each vertex.
	// This is the size of each bin to store them.
	vector<Index> binSizes(nverts, 0);

	for (Index itri=0; itri < ntris; ++itri) {
		Vector3i *triVerts = &m_tris[itri];
		Index a = (*triVerts)[0];
		Index b = (*triVerts)[1];
		Index c = (*triVerts)[2];
		binSizes[a] += 2;
		binSizes[b] += 2;
		binSizes[c] += 2;
	}

	// Find total memory to allocate for bins
	vector<Index> binOffsets(nverts);
	Index sizeAccum = 0;
	Index sizeMax = 0;

	for (Index ivert=0; ivert < nverts; ++ivert) {
		binOffsets[ivert] = sizeAccum;
		sizeAccum += binSizes[ivert];
		sizeMax = max(sizeMax, binSizes[ivert]);
	}

	// Populate the bins with vertex adjacency data
	vector<Index> binCounts(nverts, 0);
	vector<Index> bins(sizeAccum);

	for (Index itri=0; itri < ntris; ++itri) {
		for (int i=0; i <= 2; ++i) {
			int j = (i == 2 ? 0 : i+1);
			Index a = m_tris[itri][i];
			Index b = m_tris[itri][j];
			bins[binOffsets[a]+binCounts[a]] = b;
			bins[binOffsets[b]+binCounts[b]] = a;
			++binCounts[a];
			++binCounts[b];
		}
	}

	// Sort the vertices in each bin to find duplicates.
	auto compare = [](Index a, Index b) ->bool { return a < b; };

	vector<Index> sorted(sizeMax);
	m_vertadj.reserve(sizeAccum/2);
	m_vertadj.clear();
	m_vertadjOffsets.assign(nverts+1,0);
	m_boundVert.assign(nverts, false);
	m_vertadjMax = 0;

	for (Index ivert=0; ivert < nverts; ++ivert) {
		Index offset = binOffsets[ivert];
		Index size = binSizes[ivert];
		// Sort bin to find duplicates
		if ( size < 100 ) {
			// Quicker to do O(n^2) selection sort for small arrays,
			// also no temporary array allocation under the hood.
			selection_sort<Index>(bins, sorted, offset, offset+size-1, compare);
		}
		else {
			sorted.assign(bins.begin()+offset, bins.begin()+offset+size);
			std::sort(sorted.begin(), sorted.end(), compare);
		}
		m_vertadjOffsets[ivert] = m_vertadj.size();
		Index aprev = sorted[0];
		m_vertadj.push_back(aprev);
		Index nrepeat = 0;
		bool boundVert = false;
		for (Index i=1; i < size; ++i) {
			Index a = sorted[i];
			if ( aprev == a ) {
				++nrepeat;
			}
			else {
				m_vertadj.push_back(a);
				if ( nrepeat == 0 ) boundVert = true;
				nrepeat = 0;
////debug
//				m_debug.push_back(pair<Vector3f,Vector3f>(
//									  m_verts[ivert], m_verts[a]));
			}
			aprev = a;
		}
		Index count = m_vertadj.size() - m_vertadjOffsets[ivert];
		m_vertadjMax = max(m_vertadjMax, count);

		if ( boundVert ) {
			m_boundVert[ivert] = true;
////debug
//			Vector3f u = Vector3f(1e-2,1e-2,1e-2);
//			m_debug.push_back(pair<Vector3f,Vector3f>(
//								  m_verts[ivert], m_verts[ivert]+u));
		}

	}

	m_vertadjOffsets[nverts] = m_vertadj.size();

	m_staleGPUTopoData = true;

//	timer(1, c_start, c_end); //!***

	m_msgLogger->logMessage("Done!",1);

}

//---------------------------------------------------------

void Mesh::buildVertToTriMap()
{
	QMutexLocker locker(&m_recMutex);

	Index nverts = m_verts.size();
	Index ntris = m_tris.size();

	m_vertTri.resize(nverts);

	for (Index itri=0; itri < ntris; ++itri){
		for (int i=0; i<3; ++i){
			Index ivert = m_tris[itri][i];
			m_vertTri[ivert] = itri;
		}
	}

}
//---------------------------------------------------------
// Reorder vertex indices and all data that depends on them
// s.t. they are localized into 'patches' on the mesh - i.e.
// vertices that are topological near tend to have indices
// that are near as well.
// This is to improve cache spatial locality of vertices
// that are topologically near, allowing CUDA blocking into
// shared memory to squeeze out some extra performance
// during smoothing (about 50%-70%).


void Mesh::localizeVertIndices(size_t maxPatchSize)
{
	QMutexLocker locker(&m_recMutex);

	if ( m_vertadj.size() == 0 ) buildVertAdjacency();

	m_msgLogger->logMessage("Localizing vertex indices...");

	Index nverts = m_verts.size();
	Index ntris = m_tris.size();

	//Get the mapping reordering vertices
	vector<Index> map(nverts, -1);
	Index nvertsMapped = 0;
	vector<Index> patch;
	patch.reserve(maxPatchSize);

	for (Index ivertStart=0; ivertStart < nverts; ++ivertStart) {
		if ( map[ivertStart] >= 0 ) continue;
		map[ivertStart] = nvertsMapped;
		++nvertsMapped;
		patch.clear();
		patch.push_back(ivertStart);
		size_t npatch = 1;
		size_t patchIdx = 0;
////debug
//		float scale = 2e-2*(1.0f+float(rand())/float(RAND_MAX));
//		Vector3f dir = scale*Vector3f::Random();
		while ( true ) {
			if ( npatch > maxPatchSize ) break;
			if ( patchIdx >= patch.size() ) break;
			Index ivert = patch[patchIdx];
			Index offset = m_vertadjOffsets[ivert];
			Index offsetNxt = m_vertadjOffsets[ivert+1];
			for (Index i=offset; i < offsetNxt; ++i) {
				if ( npatch > maxPatchSize ) break;
				Index ivertadj = m_vertadj[i];
				if ( map[ivertadj] >= 0 ) continue;
				map[ivertadj] = nvertsMapped;
				++nvertsMapped;
				patch.push_back(ivertadj);
				++npatch;
////debug
//				m_debug.push_back(pair<Vector3f,Vector3f>(
//									  m_verts[ivertadj],
//									  m_verts[ivertadj]+dir));
			}
			++patchIdx;

		}
	}

	// Get inverse map
	vector<Index> mapinv(nverts);
	for (Index ivert=0; ivert < nverts; ++ivert) {
		mapinv[map[ivert]] = ivert;
	}


	// Remap vertex coordinates
	vector<Vector3f> m_vertsOld = m_verts;
	for (Index ivert=0; ivert < nverts; ++ivert) {
		m_verts[ivert] = m_vertsOld[mapinv[ivert]];
	}
	m_vertsOld.clear();
	m_vertsOld.shrink_to_fit(); // clear memory

	// Remap vertex normals
	vector<Vector3f> m_normsOld = m_norms;
	for (Index ivert=0; ivert < nverts; ++ivert) {
		m_norms[ivert] = m_normsOld[mapinv[ivert]];
	}
	m_normsOld.clear();
	m_normsOld.shrink_to_fit(); // clear memory


	// Remap vertex -> containing triangle pointers
	if ( m_vertTri.size() > 0 ) {
		vector<Index> m_vertTriOld = m_vertTri;
		for (Index ivert=0; ivert < nverts; ++ivert) {
			m_vertTri[ivert] = m_vertTriOld[mapinv[ivert]];
		}
//		m_vertTriOld.clear();
//		m_vertTriOld.shrink_to_fit(); // clear memory
	}

	// Remap triangle vertex indices
	for (Index itri=0; itri < ntris; ++itri){
		for (int i=0; i<3; ++i){
			Index ivert = m_tris[itri][i];
			m_tris[itri][i] = map[ivert];
		}
	}

	// Remap boundary vertex flags
	if ( m_boundVert.size() > 0 ) {
		vector<bool> m_boundVertOld = m_boundVert;
		for (Index ivert=0; ivert < nverts; ++ivert) {
			m_boundVert[ivert] = m_boundVertOld[mapinv[ivert]];
		}
		m_boundVertOld.clear();
		m_boundVertOld.shrink_to_fit(); // clear memory
	}

	// Remap vertex adjacency
	vector<Index> m_vertadjOld = m_vertadj;
	vector<Index> m_vertadjOffsetsOld = m_vertadjOffsets;
	Index offset = 0;
	for (Index ivert=0; ivert < nverts; ++ivert) {
		m_vertadjOffsets[ivert] = offset;
		Index ivertmap = mapinv[ivert];
		Index offsetOld = m_vertadjOffsetsOld[ivertmap];
		Index offsetOldNxt = m_vertadjOffsetsOld[ivertmap+1];
		Index size = offsetOldNxt - offsetOld;
		for (Index i=0; i < size; ++i) {
			Index ivertadj = m_vertadjOld[offsetOld+i];
			m_vertadj[offset+i] = map[ivertadj];
		}
		offset += size;
	}

	m_localizedVerts = true;
	m_staleGPUGeomData = true;
	m_staleGPUTopoData = true;

	m_msgLogger->logMessage("Done!",1);

}
////---------------------------------------------------------

//void Mesh::noiseMesh(int nSweeps)
//{
//	QMutexLocker locker(&m_recMutex);

////	Index nverts = m_verts.size();
//	Index ntris = m_tris.size();

//	if ( m_norms.size() == 0 ) approxMeshNorms();
//	if ( m_boundVert.size() == 0 ) buildTriangleAdjacency();

//	m_msgLogger->logMessage("Adding random noise...");

//	vector<Index> tris(ntris);
//	for (Index i=0; i < ntris; ++i){ tris[i] = i; }
//	vector<Index> sides(3);
//	for (Index i=0; i < 3; ++i){ sides[i] = i; }

//	for (int isweep=1; isweep <= nSweeps; ++isweep){
//		std::random_shuffle(tris.begin(), tris.end());

//		for (Index i=0; i < ntris; ++i) {
//			Index itri = tris[i];
//			std::random_shuffle(sides.begin(), sides.end());
//			for (Index j=0; j < 3; ++j) {
//				Index a = m_tris[itri][sides[j]];
//				Index b = m_tris[itri][sides[j != 2 ? j+1 : 0]];
//				Index c = m_tris[itri][sides[j != 0 ? j-1 : 2]];
//				if ( m_boundVert[a] ) continue;
//				Vector3f na = m_norms[a];
//				Vector3f nb = m_norms[b];
//				Vector3f nc = m_norms[c];
//				Vector3f vab = (m_verts[b] - m_verts[a]);
//				Vector3f vac = (m_verts[c] - m_verts[a]);
//				Vector3f bary = Vector3f::Random() + Vector3f(1.1f,1.1f,1.1f);
//				bary = bary/(bary.sum());
//				Vector3f w = bary(1)*vab + bary(2)*vac;
//				Vector3f pa = w - w.dot(na)*na;
//				Vector3f pb = (w - vab) - (w - vab).dot(nb)*nb;
//				Vector3f pc = (w - vac) - (w - vac).dot(nc)*nc;
//				Vector3f w1 = bary(0)*pa + bary(1)*(pb+vab) + bary(2)*(pc+vac);
//				m_verts[a] = w1 + m_verts[a];
//				Vector3f n = bary(0)*na + bary(1)*nb + bary(2)*nc;
//				m_norms[a] = n.normalized();
//			}
//		}
//	}
//
//  staleGPUTopoData = true;
//
//	m_msgLogger->logMessage("Done!",1);

//}



////---------------------------------------------------------

//void Mesh::smoothMesh(int nSweeps)
//{
//	QMutexLocker locker(&m_recMutex);

//	Index nverts = m_verts.size();
//	Index ntris = m_tris.size();
//	const Vector3f baryCent = Vector3f(1.0f,1.0f,1.0f)/3.0f;
//	const Vector3f zeros = Vector3f(0.0f,0.0f,0.0f);

//	if ( m_norms.size() == 0 ) approxMeshNorms();
//	if ( m_boundVert.size() == 0 ) buildTriangleAdjacency2();

//	m_msgLogger->logMessage("Smoothing...");

//	vector<Vector3f> newVerts(nverts, zeros);
//	vector<Vector3f> newNorms(nverts, zeros);
//	vector<Index> counts(nverts, 0);

//	for (Index itri=0; itri < ntris; ++itri) {
//		for (Index j=0; j <= 2 ; ++j) { ++counts[ m_tris[itri][j] ]; }
//	}

////	std::clock_t c_start, c_end;
////	timer(0, c_start, c_end);

//	for (int isweep=1; isweep <= nSweeps; ++isweep){
//		for (Index itri=0; itri < ntris; ++itri) {
//			for (Index j=0; j <= 2; ++j) {
//				Index a = m_tris[itri][j];
//				Index b = m_tris[itri][j != 2 ? j+1 : 0];
//				Index c = m_tris[itri][j != 0 ? j-1 : 2];
//				if ( m_boundVert[a] ) continue;
//				Vector3f na = m_norms[a];
//				Vector3f nb = m_norms[b];
//				Vector3f nc = m_norms[c];
//				Vector3f vab = (m_verts[b] - m_verts[a]);
//				Vector3f vac = (m_verts[c] - m_verts[a]);
//				Vector3f w = baryCent(1)*vab + baryCent(2)*vac;
//				Vector3f pa = w - w.dot(na)*na;
//				Vector3f pb = (w - vab) - (w - vab).dot(nb)*nb;
//				Vector3f pc = (w - vac) - (w - vac).dot(nc)*nc;
//				Vector3f w1 =
//						baryCent(0)*pa + baryCent(1)*(pb+vab) + baryCent(2)*(pc+vac);
//				Vector3f n = baryCent(0)*na + baryCent(1)*nb + baryCent(2)*nc;
//				n.normalize();
//				newVerts[a] += w1;
//				newNorms[a] += n;
//			}
//		}

//		for (Index i=0; i < nverts; ++i) {
//			m_verts[i] = newVerts[i]/counts[i] + m_verts[i];
//			m_norms[i] = newNorms[i].normalized();
//		}
//		newVerts.assign(nverts, zeros);
//		newNorms.assign(nverts, zeros);
//	}

//  staleGPUGeomData = true;


////	timer(1, c_start, c_end);

//	m_msgLogger->logMessage("Done!",1);

//}




//---------------------------------------------------------

void Mesh::noiseMesh(int nSweeps)
{
	QMutexLocker locker(&m_recMutex);

	Index nverts = m_verts.size();
	const Vector3f zeros = Vector3f(0.0f,0.0f,0.0f);

	if ( m_norms.size() == 0 ) approxMeshNorms();
	if ( m_vertadj.size() == 0 ) buildVertAdjacency();

	m_msgLogger->logMessage("Adding random noise...");

	vector<Vector3f> newVerts(nverts, zeros);
	vector<Vector3f> newNorms(nverts, zeros);

//	std::clock_t c_start, c_end;//!***
//	timer(0, c_start, c_end);

	const bool useOpenMP = true;

#pragma omp parallel if (useOpenMP) default(shared)
{

	for (int isweep=1; isweep <= nSweeps; ++isweep){
#pragma omp single
{
		m_msgLogger->logMessage(
					"Sweep " + QString::number(isweep) + "...", isweep > 1 ? 0 : 2);
}
#pragma omp for schedule(static)
		for (Index ivert=0; ivert < nverts; ++ivert) {
			if ( m_boundVert[ivert] ) continue;
			Vector3f v0 = m_verts[ivert];
			Vector3f n0 = m_norms[ivert];

			Index offset = m_vertadjOffsets[ivert];
			Index offsetNxt = m_vertadjOffsets[ivert+1];
			Index size = offsetNxt - offset;
			Index i = rand() % size;
			float t = float(rand())/float(RAND_MAX);
			Index a = m_vertadj[offset+i];
			Vector3f v = m_verts[a];
			Vector3f n = m_norms[a];
			Vector3f w = v - v0;
			Vector3f p0 = w.dot(n0)*n0;
			Vector3f p = w.dot(n)*n;
			Vector3f w1 = t*((1.0f-t)*(p - p0) + w);
			Vector3f n1 = t*n + (1.0f-t)*n0;
			newVerts[ivert] = w1 + v0;
			newNorms[ivert] = n1.normalized();
		}

#pragma omp for schedule(static)
		for (Index ivert=0; ivert < nverts; ++ivert) {
			if ( m_boundVert[ivert] ) continue;
			m_verts[ivert] = newVerts[ivert];
			m_norms[ivert] = newNorms[ivert];
			newVerts[ivert] = zeros;
			newNorms[ivert] = zeros;
		}
	}

} //Parallel

	m_staleGPUGeomData = true;

//	timer(1, c_start, c_end);

	m_msgLogger->logMessage("Done!",1);

}


//---------------------------------------------------------

void Mesh::smoothMesh(int nSweeps)
{
	QMutexLocker locker(&m_recMutex);

	if ( m_mesh_smooth_GPU == nullptr ) {
		smoothMesh_CPU(nSweeps);
	}
	else {
		smoothMesh_GPU(nSweeps);
		if ( m_mesh_smooth_GPU == nullptr ) {
			// GPU smoothing failed
			smoothMesh_CPU(nSweeps);
		}
	}

}

//---------------------------------------------------------

void Mesh::smoothMesh_CPU(int nSweeps)
{
	QMutexLocker locker(&m_recMutex);

	Index nverts = m_verts.size();
	const Vector3f zeros = Vector3f(0.0f,0.0f,0.0f);

	if ( m_norms.size() == 0 ) approxMeshNorms();
	if ( m_vertadj.size() == 0 ) buildVertAdjacency();

	m_msgLogger->logMessage("Smoothing (CPU)...");

	vector<Vector3f> newVerts = m_verts;
	vector<Vector3f> newNorms = m_norms;

	vector<Vector3f> *p_verts = &m_verts;
	vector<Vector3f> *p_norms = &m_norms;
	vector<Vector3f> *p_newVerts = &newVerts;
	vector<Vector3f> *p_newNorms = &newNorms;

//	std::clock_t c_start, c_end;//!***
//	timer(0, c_start, c_end);

	const bool useOpenMP = true;

#pragma omp parallel if (useOpenMP) default(shared)
{

	for (int isweep=1; isweep <= nSweeps; ++isweep){		
#pragma omp single
{
		m_msgLogger->logMessage(
					"Sweep " + QString::number(isweep) + "...", isweep > 1 ? 0 : 2);
}
#pragma omp for schedule(static)
		for (Index ivert=0; ivert < nverts; ++ivert) {
			if ( m_boundVert[ivert] ) continue;
			Vector3f v0 = (*p_verts)[ivert];
			Vector3f n0 = (*p_norms)[ivert];
			Vector3f v1 = zeros;
			Vector3f n1 = zeros;

			Index offset = m_vertadjOffsets[ivert];
			Index offsetNxt = m_vertadjOffsets[ivert+1];
			Index size = offsetNxt - offset;
			for (Index i=offset; i < offsetNxt; ++i) {
				Index a = m_vertadj[i];
				Vector3f v = (*p_verts)[a];
				Vector3f n = (*p_norms)[a];
				v = v - v0;
				float vn0 = v.dot(n0);
				float vn = v.dot(n);
				v1 += 0.25f*(vn*n - vn0*n0) + 0.5f*v;
				n = 0.5f*(n + n0);
				n1 += n.normalized();
			}
			(*p_newVerts)[ivert] = v1/size + v0;
			(*p_newNorms)[ivert] = n1.normalized();
		}

#pragma omp single
{
			std::swap(p_verts, p_newVerts);
			std::swap(p_norms, p_newNorms);
}

	} // sweeps


	if ( nSweeps%2 != 0 ) {
#pragma omp for schedule(static)
		for (Index ivert=0; ivert < nverts; ++ivert) {
			m_verts[ivert] = newVerts[ivert];
			m_norms[ivert] = newNorms[ivert];
		}
	}

} //Parallel

	m_staleGPUGeomData = true;

//	timer(1, c_start, c_end);

	m_msgLogger->logMessage("Done!",1);

}

//---------------------------------------------------------

void Mesh::smoothMesh_GPU(int nSweeps)
{
	QMutexLocker locker(&m_recMutex);

	size_t nverts = m_verts.size();
	const Vector3f zeros = Vector3f(0.0f,0.0f,0.0f);

	if ( m_norms.size() == 0 ) approxMeshNorms();
	if ( m_vertadj.size() == 0 ) buildVertAdjacency();

	m_msgLogger->logMessage("Smoothing (GPU)...");

//	std::clock_t c_start, c_end;//!***
//	timer(0, c_start, c_end);

	if ( m_staleGPUGeomData ) {
		m_vertsx_GPU.resize(nverts);
		m_vertsy_GPU.resize(nverts);
		m_vertsz_GPU.resize(nverts);
		m_normsx_GPU.resize(nverts);
		m_normsy_GPU.resize(nverts);
		m_normsz_GPU.resize(nverts);

		for (size_t ivert = 0; ivert < nverts; ++ivert) {
			m_vertsx_GPU[ivert] = m_verts[ivert][0];
			m_vertsy_GPU[ivert] = m_verts[ivert][1];
			m_vertsz_GPU[ivert] = m_verts[ivert][2];
		}
		for (size_t ivert = 0; ivert < nverts; ++ivert) {
			m_normsx_GPU[ivert] = m_norms[ivert][0];
			m_normsy_GPU[ivert] = m_norms[ivert][1];
			m_normsz_GPU[ivert] = m_norms[ivert][2];
		}
		m_staleGPUGeomData = false;
	}

	if ( m_staleGPUTopoData ) {
		size_t nadj = m_vertadj.size();

		m_vertidxs_GPU.resize(nverts);
		m_vertadj_GPU.resize(nadj);
		m_vertadjOffsets_GPU.resize(nverts+1);

		size_t offset = 0;
		size_t ninterior = 0;
		for(size_t ivert=0; ivert < nverts; ++ivert) {
			if ( m_boundVert[ivert] ) continue;
			m_vertidxs_GPU[ninterior] = ivert;
			m_vertadjOffsets_GPU[ninterior] = offset;
			size_t offsetOrig = m_vertadjOffsets[ivert];
			size_t offsetOrigNxt = m_vertadjOffsets[ivert+1];
			size_t size = offsetOrigNxt - offsetOrig;
			for (size_t i=0; i < size; ++i) {
				m_vertadj_GPU[offset+i] = m_vertadj[offsetOrig+i];
			}
			offset += size;
			++ninterior;
		}
		m_vertadjOffsets_GPU[ninterior] = offset;

		m_vertidxs_GPU.resize(ninterior);
		m_vertadj_GPU.resize(offset);
		m_vertadjOffsets_GPU.resize(ninterior+1);

		m_staleGPUTopoData = false;
	}

	size_t nidxs = m_vertidxs_GPU.size();
	size_t nadjidxs = m_vertadj_GPU.size();

	bool success;
	m_mesh_smooth_GPU(
				nSweeps, nverts, nidxs, nadjidxs, 128,
				m_localizedVerts,
				m_vertidxs_GPU.data(), m_vertadj_GPU.data(),
				m_vertadjOffsets_GPU.data(),
				m_vertsx_GPU.data(), m_vertsy_GPU.data(),
				m_vertsz_GPU.data(),
				m_normsx_GPU.data(), m_normsy_GPU.data(),
				m_normsz_GPU.data(), success);

	if ( success ) {
		for (size_t ivert = 0; ivert < nverts; ++ivert) {
			m_verts[ivert][0] = m_vertsx_GPU[ivert];
			m_verts[ivert][1] = m_vertsy_GPU[ivert];
			m_verts[ivert][2] = m_vertsz_GPU[ivert];
		}

		for (size_t ivert = 0; ivert < nverts; ++ivert) {
			m_norms[ivert][0] = m_normsx_GPU[ivert];
			m_norms[ivert][1] = m_normsy_GPU[ivert];
			m_norms[ivert][2] = m_normsz_GPU[ivert];
		}
		//	timer(1, c_start, c_end);

		m_msgLogger->logMessage("Done!",1);
	}
	else {
		m_msgLogger->logMessage("Failed.",1);
	}
}

//---------------------------------------------------------
// Try to get as close as possible to a point in space by walking
// along the surface of the mesh from a chosen vertex, first doing
// a quick and approximate walk along mesh triangle edges as far as
// possible, then a more thorough walk along interiors of triangles.
//
// As long as the point is relatively close to the starting vertex
// and the mesh is well-behaved around it, the result should be a
// nearest approach of the point and the mesh, arrived at quicker
// than a spatial index can - and without needing to rebuild or
// modify spatial index when the mesh vertex positions are changed
// during smoothing.


void Mesh::meshWalkTo(const Vector3f& pointEnd, Index ivertStart,
				Index& itriEnd, Vector3f& baryEnd, Vector3f *triPointEnd)
{
	if ( m_vertTri.size() == 0 ) buildVertToTriMap();

	Index ivertEnd;
	edgeWalkTo(pointEnd, ivertStart, ivertEnd);
	triWalkTo(pointEnd, m_vertTri[ivertEnd], itriEnd, baryEnd);

	if ( triPointEnd != nullptr ) {
		Vector3i iv = m_tris[itriEnd];
		*triPointEnd =
				baryEnd(0)*m_verts[iv(0)] +
				baryEnd(1)*m_verts[iv(1)] +
				baryEnd(2)*m_verts[iv(2)];
	}

}

//---------------------------------------------------------

void Mesh::edgeWalkTo(const Vector3f& pointEnd, Index ivertStart,
					  Index& ivertEnd)
{
	if ( m_vertadj.size() == 0 ) buildVertAdjacency();

	auto distsq_vert_to_end = [&](Index ivert) ->double {
		Vector3f w = pointEnd - m_verts[ivert];
		return w(0)*w(0) + w(1)*w(1) + w(2)*w(2);
	};

	Index ivertCurr = ivertStart;
	double distsqCurr = distsq_vert_to_end(ivertCurr);
	Index ivertPrev = -1;

	while (true){
		Index ivertCurr0 = ivertCurr;
		Index offset = m_vertadjOffsets[ivertCurr0];
		Index offsetNxt = m_vertadjOffsets[ivertCurr0+1];
		for (Index i=offset; i < offsetNxt; ++i) {
			Index ivert = m_vertadj[i];
			if (ivert == ivertPrev ) continue;
			double distsq = distsq_vert_to_end(ivert);
			if ( distsq < distsqCurr ) {
				ivertCurr = ivert;
				distsqCurr = distsq;
			}
		}

		if ( ivertCurr == ivertCurr0 ) break;

////debug
//		Vector3f v = m_verts[ivertCurr0];
//		Vector3f w = m_verts[ivertCurr];
//		m_debug.push_back(pair<Vector3f,Vector3f>(v,w));

		ivertPrev = ivertCurr0;
	}


	ivertEnd = ivertCurr;

////debug
//	m_debug.push_back(pair<Vector3f,Vector3f>(m_verts[ivertEnd], pointEnd));

}

//---------------------------------------------------------

void Mesh::triWalkTo(const Vector3f& pointEnd, Index itriStart,
					 Index& itriEnd, Vector3f& baryEnd)
{
	if ( m_triadj.size() == 0 ) buildTriAdjacency();

	auto proj_end_to_tri = [&](
			Index itri,
			Vector3f& proj, Vector3f& bary, double& distsq, bool& hovering) ->void {
		bool degenerate;
		Vector3i *triVerts = &m_tris[itri];
		Vector3f va = m_verts[ (*triVerts)[0] ];
		Vector3f vb = m_verts[ (*triVerts)[1] ];
		Vector3f vc = m_verts[ (*triVerts)[2] ];
		proj = point_project_triangle(
					pointEnd, va, vb, vc, hovering, degenerate, &bary);
		 distsq = point_dist_sq(proj,pointEnd);
	};


	Index itriBest = itriStart;
	Vector3f projBest, baryBest;
	double distsqBest;
	bool hovering;
	proj_end_to_tri(itriBest, projBest, baryBest, distsqBest, hovering);
	Index itriCurr = itriBest;
	Index itriPrev = -1;
	Index itriAnchor = itriBest;
	Index itriAdjA0, itriAdjB0, itriAdjC0;
	Index iv0a, iv0b, iv0c;
	double distsqAnchor = distsqBest;
	bool newAnchor = true;

	while (true){
		Index itriNext = -1;
		Vector3f projNext, baryNext;
		double distsqNext = double_infinity;

		if ( newAnchor ) {
			itriAdjA0 = m_triadj[itriAnchor][0];
			itriAdjB0 = m_triadj[itriAnchor][1];
			itriAdjC0 = m_triadj[itriAnchor][2];			
			iv0a = m_tris[itriAnchor][0];
			iv0b = m_tris[itriAnchor][1];
			iv0c = m_tris[itriAnchor][2];
		}

		for (int iside=0; iside <= 2; ++iside){
			Index itri = m_triadj[itriCurr][iside];
			if (itri < 0 || itri == itriPrev ) continue;
			if ( !newAnchor ) {
				if ( itri == itriAnchor || itri == itriAdjA0 ||
					 itri == itriAdjB0  || itri == itriAdjC0 ) continue;
			}

			double distsq;
			Vector3f proj, bary;
			proj_end_to_tri(itri, proj, bary, distsq, hovering);
			if ( distsq < distsqNext ) {
				itriNext = itri;
				projNext = proj;
				baryNext = bary;
				distsqNext = distsq;
			}
		}

		if ( itriNext == -1 ) break;

		if ( distsqNext < distsqBest  ) {
			itriBest = itriNext;
			projBest = projNext;
			baryBest = baryNext;
			distsqBest = distsqNext;
		}

		if ( distsqNext < distsqAnchor ) {
			newAnchor = true;
			itriAnchor = itriNext;
			distsqAnchor = distsqNext;
		}
		else {
			newAnchor = false;
			Index iva = m_tris[itriNext][0];
			Index ivb = m_tris[itriNext][1];
			Index ivc = m_tris[itriNext][2];
			bool anchorNotAdj =
					(iv0a!=iva && iv0a!=ivb && iv0a!=ivc) &&
					(iv0b!=iva && iv0b!=ivb && iv0b!=ivc) &&
					(iv0c!=iva && iv0c!=ivb && iv0c!=ivc);
			if ( anchorNotAdj ) break;
		}
////debug
//		Vector3f va = m_verts[ m_tris[itriCurr][0] ];
//		Vector3f vb = m_verts[ m_tris[itriCurr][1] ];
//		Vector3f vc = m_verts[ m_tris[itriCurr][2] ];
//		Vector3f wa = m_verts[ m_tris[itriNext][0] ];
//		Vector3f wb = m_verts[ m_tris[itriNext][1] ];
//		Vector3f wc = m_verts[ m_tris[itriNext][2] ];
//		m_debug.push_back(pair<Vector3f,Vector3f>(
//							  (va+vb+vc)/3.0f, (wa+wb+wc)/3.0f));

		itriPrev = itriCurr;
		itriCurr = itriNext;

	}


	itriEnd = itriBest;
	baryEnd = baryBest;

////debug
//	m_debug.push_back(pair<Vector3f,Vector3f>(projBest, pointEnd));

}

//---------------------------------------------------------















