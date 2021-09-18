//-----------------------------------------------------------
//     Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//     See LICENSE included with this distribution.

#ifndef RECONSTRUCTTHREAD_H
#define RECONSTRUCTTHREAD_H

#include "constants.h"

#include <QMutex>
#include <QObject>
#include <QVector3D>

class Mesh;
class BoundBox;


class MeshWorker : public QObject
{
	Q_OBJECT

public:
	explicit MeshWorker(
			Mesh& mesh, QObject *parent = nullptr);
	~MeshWorker();

public slots:
	void generateMesh(
			QVector3D norm, size_t ndim,
			const std::function<float(float xu, float xv)> heightFun);
	void approxMeshNorms();
	void noiseMesh(int nSweeps);
	void smoothMesh(int nSweeps);

signals:
	void finished(int opt=0);

private:
	Mesh *m_mesh;
	QMutex m_mutex;
};


#endif // RECONSTRUCTTHREAD_H
