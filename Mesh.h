//-----------------------------------------------------------
//  Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//  See LICENSE included with this distribution.

#ifndef MESH_H
#define MESH_H

#include "constants.h"
//#include "MessageLogger.h"

#include <Eigen/Core>
#include <qopengl.h>
#include <QObject>
#include <QRecursiveMutex>
//#include <QMutexLocker>
#include <vector>
#include <functional>

class MessageLogger;
class BoundBox;

class Mesh : public QObject
{
	Q_OBJECT

	template<typename T> using vector = std::vector<T>;
	using Index = Eigen::Index;
	using Vector3f = Eigen::Vector3f;
	using Vector3i = Eigen::Vector3i;
	using Matrix3f = Eigen::Matrix3f;

public:
	Mesh(MessageLogger *msgLogger = nullptr, QObject *parent = nullptr);
	~Mesh();
	void clear();
	const Vector3f vert(Index idx) const { return m_verts[idx]; }
	const GLfloat* vertGLData();
	const GLuint* triGLData();
	const GLfloat* normGLData(float scale);
	const GLfloat* debugGLData();

	Index vertCount() const { return m_verts.size(); }
	Index triCount() const { return m_tris.size(); }
//	Index tetCount() const { return m_tets.size(); }
	Index debugCount() const { return m_debug.size(); }

	Index addVert(const Vector3f& v, const Vector3f& n, bool threadSafe = false);
	Index addTri(const Vector3i& iv, bool threadSafe = false);
	void fromPlaneHeight(
			float normx, float normy, float normz, size_t ndim,
			const std::function<float(float xu, float xv)> heightFun = nullptr);

	void approxMeshNorms();
//	void buildTriAdjacency2();
	void buildTriAdjacency();
	void buildVertAdjacency();
	void buildVertToTriMap();
	void localizeVertIndices(size_t patchSize);
	void localizeTriIndices(size_t patchSize);

	void noiseMesh(int nSweeps);
	void smoothMesh(int nSweeps);
	void meshWalkTo(const Vector3f& pointEnd, Index ivertStart,
					Index& itriEnd, Vector3f& baryEnd,
					Vector3f *triPointEnd = nullptr);
	void edgeWalkTo(const Vector3f& pointEnd, Index ivertStart,
					Index& ivertEnd);
	void triWalkTo(const Vector3f& pointEnd, Index itriStart,
				   Index& itriEnd, Vector3f& baryEnd);
private:

	vector<Vector3f> m_verts;
	vector<Vector3f> m_norms;
	vector<Vector3i> m_tris;
//	vector<Vector3i> m_tets;
	vector<bool> m_boundVert;
	vector<Vector3i> m_triadj;
	vector<Index> m_vertadj;
	vector<Index> m_vertadjOffsets;
	vector<Index> m_vertadjCounts;
	vector<Index> m_vertTri;
	Index m_vertadjMax;
	vector<std::pair<Vector3f,Vector3f>> m_debug;
	vector<GLfloat> m_vertGL;
	vector<GLuint> m_triGL;
	vector<GLfloat> m_normGL;
	vector<GLfloat> m_debugGL;
	MessageLogger *m_msgLogger;
	QRecursiveMutex m_recMutex;

};


#endif // MESH_H
