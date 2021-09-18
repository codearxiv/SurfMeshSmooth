/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** BSD License Usage
** Alternatively, you may use this file under the terms of the BSD license
** as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/
/****************************************************************************
**     Peter Beben: heavily modified this file for the purpose of
**     this project.
****************************************************************************/

#ifndef GLWIDGET_H
#define GLWIDGET_H

#include "Mesh.h"
#include "MeshWorker.h"
//#include "MessageLogger.h"

#include <QOpenGLWidget>
//#include <QOpenGLFunctions>
#include <QOpenGLExtraFunctions>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QMatrix4x4>
#include <QRecursiveMutex>

class MessageLogger;
class MeshWorker;

QT_BEGIN_NAMESPACE
class QThread;
QT_END_NAMESPACE

QT_FORWARD_DECLARE_CLASS(QOpenGLShaderProgram)

class GLWidget : public QOpenGLWidget, protected QOpenGLExtraFunctions
{
    Q_OBJECT

	typedef void* MeshPtr;

public:
	GLWidget(QWidget *parent = nullptr, MessageLogger *msgLogger = nullptr);
	~GLWidget();

	static bool isTransparent() { return m_transparent; }
	static void setTransparent(bool t) { m_transparent = t; }

	QSize minimumSizeHint() const override;
	QSize sizeHint() const override;

public slots:
	void setVectRotation(int angle, QVector3D v);
	void setVectTranslation(QVector3D v);
	void cleanup();
	void setMesh(MeshPtr mesh);
	void setRandomMesh(size_t ndim);
	void getMesh(MeshPtr& mesh);
	void viewGLMeshNorms(bool enabled);
	void approxMeshNorms()
	{ emit meshApproxNorms(); }
	void noiseMesh(int nSweeps)
	{ emit meshNoise(nSweeps); }
	void smoothMesh(int nSweeps)
	{ emit meshSmooth(nSweeps); }
	void setNormScale(float scale) { setGLMeshNorms(scale*m_modelSize); }
	void updateMesh(int opt=0);

signals:
	void vectRotationChanged(int angle, QVector3D v);
	void vectTranslationChanged(QVector3D v);
	void meshGenerate(QVector3D norm, size_t ndim,
					  const std::function<float(float xu, float xv)> heightFun);
	void meshApproxNorms();
	void meshNoise(int nSweeps);
	void meshSmooth(int nSweeps);
	void logMessage(const QString& text);

protected:
	void initializeGL() override;
	void paintGL() override;
	void resizeGL(int width, int height) override;
	void mousePressEvent(QMouseEvent *event) override;
	void mouseMoveEvent(QMouseEvent *event) override;
	void wheelEvent(QWheelEvent *event) override;

private:
	void setupVertexAttribs(QOpenGLBuffer vbo);
	void setGLView();
	void setGLMesh();
	void setGLMeshNorms(float scale);
	void setGLMeshDebug();
	float frobeniusNorm4x4(QMatrix4x4 M);

	int m_vRot;
	QPoint m_lastMousePos;
	//QPoint m_lastWheelPos;
	Mesh m_mesh;
	float m_aspectRatio;
	float m_nearPlane;
	float m_farPlane;
	float m_modelSize;
	static constexpr int m_alignMatCol = 4*sizeof(float);
	static constexpr int m_alignSizMat33 = 3*m_alignMatCol;
	static constexpr int m_alignSizMat44 = 4*m_alignMatCol;
	static constexpr int m_projMatLoc = 0;
	static constexpr int m_projMatSiz = m_alignSizMat44;
	static constexpr int m_movMatLoc = m_projMatLoc + m_projMatSiz;
	static constexpr int m_movMatSiz = m_alignSizMat44;
	static constexpr int m_normMatLoc = m_movMatLoc + m_movMatSiz;
	static constexpr int m_normMatSiz = m_alignSizMat33;
	static constexpr size_t m_matUboSize = 2*m_alignSizMat44 + m_alignSizMat33;
	static constexpr int m_lightPosLoc = 0;
	static constexpr int m_lightPosSiz = 4*sizeof(float);
	static constexpr int m_colorLoc = m_lightPosLoc + m_lightPosSiz;
	static constexpr int m_colorSiz = 4*sizeof(float);
	static constexpr size_t m_lightUboSize = 8*sizeof(float);
	float m_normScale;
	bool m_showMeshNorms;
	QMatrix4x4 m_proj;
	QMatrix4x4 m_camera;
	QMatrix4x4 m_world;
	QVector3D m_rotVect;
	QVector3D m_movVect;
	QOpenGLVertexArrayObject m_meshVao;
	QOpenGLVertexArrayObject m_meshNormsVao;
	QOpenGLVertexArrayObject m_meshDebugVao;
	QOpenGLBuffer m_meshVbo;
	QOpenGLBuffer m_meshEbo;
	unsigned int m_matricesUbo;
	unsigned int m_lightsUbo;
	QOpenGLBuffer m_meshNormsVbo;
	QOpenGLBuffer m_meshDebugVbo;
	QOpenGLShaderProgram *m_shader;
	QOpenGLShaderProgram *m_shaderWireframe;
	static bool m_transparent;
	MessageLogger* m_msgLogger;
	QThread* m_meshThread;
	MeshWorker* m_meshWorker;
	QRecursiveMutex m_recMutex;

};

#endif
