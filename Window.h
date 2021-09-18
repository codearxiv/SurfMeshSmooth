//-----------------------------------------------------------
//     Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//     See LICENSE included with this distribution.

#ifndef WINDOW_H
#define WINDOW_H

#include "MessageLogger.h"
#include "constants.h"

#include <QWidget>

QT_BEGIN_NAMESPACE
class QSlider;
class QPushButton;
class QMainWindow;
QT_END_NAMESPACE

class GLWidget;


class Window : public QWidget
{
    Q_OBJECT

	typedef void* MeshPtr;

public:
	Window(QMainWindow *mw, MessageLogger *msgLogger = nullptr);

public slots:
	void setMesh(MeshPtr mesh)
	{ emit meshChanged(mesh); }

	void getMesh(MeshPtr& mesh)
	{ emit meshQueried(mesh); }

	void undoMesh() { emit meshUndo(); }

	void viewGLMeshNorms(bool enabled)
	{ emit meshNormsViewGL(enabled); }

	void setRandomMesh(size_t nDim)
	{ emit meshSetRandom(nDim); }

	void approxMeshNorms()
	{ emit meshApproxNorms(); }

	void noiseMesh(int nSweeps)
	{ emit meshNoise(nSweeps); }

	void smoothMesh(int nSweeps)
	{ emit meshSmooth(nSweeps); }

	void setNormScale(float scale)
	{ emit normScaleChanged(scale); }

signals:
	void meshChanged(MeshPtr mesh);
	void meshQueried(MeshPtr& mesh);
	void meshUndo();
	void meshNormsViewGL(bool enabled);
	void meshSetRandom(size_t nDim);
	void meshApproxNorms();
	void meshNoise(int nSweeps);
	void meshSmooth(int nSweeps);
	void normScaleChanged(float scale);

protected:
    void keyPressEvent(QKeyEvent *event) override;

private:
    QSlider *createSlider();

    GLWidget *glWidget;
	QMainWindow *mainWindow;
	MessageLogger* m_msgLogger;

};

#endif
