//-----------------------------------------------------------
//  Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//  See LICENSE included with this distribution.

#include "MeshWorker.h"
#include "Mesh.h"
#include "constants.h"

#include <QVector3D>
//-------------------------------------------------------------------------
MeshWorker::MeshWorker(Mesh& mesh, QObject *parent) :
	QObject(parent)
{
	m_mesh = &mesh;
}

//-------------------------------------------------------------------------
MeshWorker::~MeshWorker()
{
}
//-------------------------------------------------------------------------
void MeshWorker::generateMesh(
		QVector3D norm, size_t ndim,
		const std::function<float(float xu, float xv)> heightFun)
{
	m_mesh->fromPlaneHeight(norm[0], norm[1], norm[2], ndim, heightFun);
	emit finished();
}
//-------------------------------------------------------------------------

void MeshWorker::approxMeshNorms()
{
	if (m_mesh->vertCount() == 0) return;
	m_mesh->approxMeshNorms();
	emit finished();

}

//-------------------------------------------------------------------------
void MeshWorker::noiseMesh(int nSweeps)
{
	if (m_mesh->vertCount() == 0) return;
	m_mesh->noiseMesh(nSweeps);

	emit finished();
}

//-------------------------------------------------------------------------
void MeshWorker::smoothMesh(int nSweeps)
{
	if (m_mesh->vertCount() == 0) return;
	m_mesh->smoothMesh(nSweeps);

	emit finished();
}

//-------------------------------------------------------------------------
