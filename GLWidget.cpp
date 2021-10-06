//-----------------------------------------------------------
//  Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//  See LICENSE included with this distribution.

#include "GLWidget.h"
#include "Mesh.h"
#include "MeshWorker.h"
#include "MessageLogger.h"
#include "constants.h"

#include <Eigen/Core>
#include <QMouseEvent>
#include <QOpenGLShaderProgram>
#include <QOpenGLFunctions>
#include <QCoreApplication>
#include <QDebug>
#include <QThread>
#include <QRecursiveMutex>
#include <math.h>

//#include <QLibrary>
#include <windows.h>

#define SHOW_MESH
#define SHOW_MESH_DBUG


bool GLWidget::m_transparent = false;

#define PROGRAM_VERTEX_ATTRIBUTE 0
#define PROGRAM_NORMAL_ATTRIBUTE 1

//---------------------------------------------------------

GLWidget::GLWidget(QWidget *parent, MessageLogger* msgLogger)
	: QOpenGLWidget(parent),
	  m_msgLogger(msgLogger),
	  m_vRot(0),
	  m_shader(nullptr),
	  m_shaderWireframe(nullptr),
	  m_rotVect(0.0f,0.0f,0.0f),
	  m_movVect(0.0f,0.0f,0.0f),
	  m_showMeshNorms(false),
	  m_aspectRatio(1.0f),
	  m_farPlane(100.0f), m_nearPlane(0.01f),
	  m_modelSize(1.0f),
	  m_meshVbo(QOpenGLBuffer::VertexBuffer),
	  m_meshEbo(QOpenGLBuffer::IndexBuffer),
	  m_meshNormsVbo(QOpenGLBuffer::VertexBuffer),
	  m_meshDebugVbo(QOpenGLBuffer::VertexBuffer),
	  m_mesh(msgLogger, parent)
{
	// --transparent causes the clear color to be transparent. Therefore, on systems that
	// support it, the widget will become transparent apart from the mesh.
	if (m_transparent) {
		QSurfaceFormat fmt = format();
		fmt.setAlphaBufferSize(8);
		setFormat(fmt);
	}


	typedef void (*GPUSmoothPrototype)(int, size_t, size_t, size_t,
									   unsigned int, bool, const size_t*,
									   const size_t*, const size_t*,
									   float*, float*, float*,
									   float*, float*, float*, bool&);
//	if (QLibrary::isLibrary("CUDAMeshSmooth.dll")) {
//		QLibrary lib("CUDAMeshSmooth.dll");
//		lib.load();
//		if (lib.isLoaded()) {
//			auto mesh_smooth_GPU = (GPUSmoothPrototype)lib.resolve("cuda_mesh_smooth");
//			m_mesh.setGPUSmoothing(mesh_smooth_GPU);
//		}
//		else { qDebug() << "Error " << lib.errorString() << "\n"; }
//	}
//	else { qDebug() << "Not a library\n"; }

	HINSTANCE hDLL = LoadLibrary(L"CUDAMeshSmooth.dll");
	if ( hDLL != nullptr ) {
		auto mesh_smooth_GPU =
				(GPUSmoothPrototype)GetProcAddress(hDLL, "cuda_mesh_smooth");
		if ( mesh_smooth_GPU != nullptr ) {
			bool success;
			mesh_smooth_GPU(
						0,0,0,0,0,
						true, nullptr, nullptr, nullptr,
						nullptr, nullptr, nullptr,
						nullptr, nullptr, nullptr, success);
			if ( success ) m_mesh.setGPUSmoothing(mesh_smooth_GPU);
		}
	}

	m_meshThread = new QThread(this);
	m_meshWorker = new MeshWorker(m_mesh);
	m_meshWorker->moveToThread(m_meshThread);

	connect(this, &GLWidget::meshGenerate,
			m_meshWorker, &MeshWorker::generateMesh);

	connect(this, &GLWidget::meshApproxNorms,
			m_meshWorker, &MeshWorker::approxMeshNorms);

	connect(this, &GLWidget::meshNoise,
			m_meshWorker, &MeshWorker::noiseMesh);

	connect(this, &GLWidget::meshSmooth,
			m_meshWorker, &MeshWorker::smoothMesh);

	connect(m_meshWorker, &MeshWorker::finished,
			this, &GLWidget::updateMesh);

	connect(m_meshThread, &QThread::finished,
			m_meshWorker, &QObject::deleteLater);

//	connect(qApp, &QCoreApplication::aboutToQuit,
//			m_meshThread, &QThread::quit);

	m_meshThread->start();

}
//---------------------------------------------------------

GLWidget::~GLWidget()
{
	m_meshThread->quit();
	m_meshThread->wait();
//	delete m_meshThread;
	delete m_meshWorker;
	cleanup();
}
//---------------------------------------------------------

QSize GLWidget::minimumSizeHint() const
{
	return QSize(50, 50);
}
//---------------------------------------------------------

QSize GLWidget::sizeHint() const
{
	return QSize(400, 400);
}
//---------------------------------------------------------

static void qNormalizeAngle(int &angle)
{
	while (angle < 0)
		angle += 360 * 16;
	while (angle > 360 * 16)
		angle -= 360 * 16;
}

//---------------------------------------------------------

void GLWidget::setVectRotation(int angle, QVector3D v)
{
	qNormalizeAngle(angle);
	m_vRot = angle;
	m_rotVect = v;
	m_movVect = QVector3D(0.0f,0.0f,0.0f);
	emit vectRotationChanged(angle, v);
	update();
}
//---------------------------------------------------------

void GLWidget::setVectTranslation(QVector3D v)
{
	m_vRot = 0.0f;
	m_rotVect = QVector3D(0.0f,0.0f,0.0f);
	m_movVect = v;
	emit vectTranslationChanged(v);
	update();
}
//---------------------------------------------------------

void GLWidget::cleanup()
{
	makeCurrent();
	m_meshVbo.destroy();
	m_meshEbo.destroy();
	m_meshNormsVbo.destroy();
	m_meshDebugVbo.destroy();
	if (m_shader != nullptr) {
		delete m_shader;
		m_shader = nullptr;
	}
	if (m_shaderWireframe != nullptr) {
		delete m_shaderWireframe;
		m_shaderWireframe = nullptr;
	}
	doneCurrent();
}
//---------------------------------------------------------

static const char *vertexShaderSource =
		"#version 420 core\n"
		"layout (std140, binding=0) uniform Matrices\n"
		"{\n"
		"   mat4 projMatrix;\n"
		"   mat4 movMatrix;\n"
		"   mat3 normMatrix;\n"
		"};\n"
		"in highp vec4 vertex;\n"
		"in mediump vec3 normal;\n"
		"out highp vec3 vertPos;\n"
		"out mediump vec3 vertNormal;\n"
		"void main() {\n"
		"   vertPos = vertex.xyz;\n"
		"   vertNormal = normMatrix * normal;\n"
		"   gl_Position = projMatrix * movMatrix * vertex;\n"
		"}\n";

static const char *geometryShaderSource =		
		//Render filled mesh wireframe in a single-pass
		//'Single-Pass Wireframe Rendering':
		//http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/4884/pdf/imm4884.pdf
		//https://strattonbrazil.blogspot.com/2011/09/single-pass-wireframe-rendering_10.html
		"#version 420 core\n"
		"layout(triangles) in;\n"
		"layout(triangle_strip, max_vertices = 3) out;\n"
		"in vec3 vertPos[3];\n"
		"in vec3 vertNormal[3];\n"
		"out vec3 normal;\n"
		"out vec3 pos;\n"
		"uniform vec2 WIN_SCALE;\n"
		"noperspective out vec3 dist;\n"
		"void main(void){\n"
		"   vec2 p0 = WIN_SCALE * gl_in[0].gl_Position.xy/gl_in[0].gl_Position.w;\n"
		"   vec2 p1 = WIN_SCALE * gl_in[1].gl_Position.xy/gl_in[1].gl_Position.w;\n"
		"   vec2 p2 = WIN_SCALE * gl_in[2].gl_Position.xy/gl_in[2].gl_Position.w;\n"
		"   vec2 v0 = p2-p1;\n"
		"   vec2 v1 = p2-p0;\n"
		"   vec2 v2 = p1-p0;\n"
		"   float area = abs(v1.x*v2.y - v1.y * v2.x);\n"
		"   dist = vec3(area/length(v0),0,0);\n"
		"   pos = vertPos[0];\n"
		"   normal = vertNormal[0];\n"
		"   gl_Position = gl_in[0].gl_Position;\n"
		"   EmitVertex();\n"
		"   dist = vec3(0,area/length(v1),0);\n"
		"   pos = vertPos[1];\n"
		"   normal = vertNormal[1];\n"
		"   gl_Position = gl_in[1].gl_Position;\n"
		"   EmitVertex();\n"
		"   dist = vec3(0,0,area/length(v2));\n"
		"   pos = vertPos[2];\n"
		"   normal = vertNormal[2];\n"
		"   gl_Position = gl_in[2].gl_Position;\n"
		"   EmitVertex();\n"
		"   EndPrimitive();\n"
		"}\n";

static const char *fragmentShaderSourceWireframe =
		"#version 420 core\n"
		"in highp vec3 pos;\n"
		"in mediump vec3 normal;\n"
		"out highp vec4 fragColor;\n"
		"layout (std140, binding=1) uniform Lights\n"
		"{\n"
		"   vec3 lightPos;\n"
		"   vec3 vertColor;\n"
		"};\n"
		"noperspective in vec3 dist;\n"
		"void main() {\n"
		"   float mediump nearD = min(min(dist[0],dist[1]),dist[2]);\n"
		"   float mediump t = exp2(-5.0*nearD*nearD);\n"
		"   highp vec3 L = normalize(lightPos - pos);\n"
		"   highp float NL = max(abs(dot(normalize(normal), L)), 0.0);\n"
		"   highp vec3 color = vertColor;\n"
		"   highp vec3 col = clamp(color * (0.3 + 0.8*NL), 0.0, 1.0);\n"
		"   fragColor = t*vec4(0.1,0.1,0.1,1.0) + (1.0-t)*vec4(col, 1.0);\n"
//		"   fragColor = vec4(col, 1.0);\n"
		"}\n";

static const char *fragmentShaderSource =
		"#version 420 core\n"
		"in highp vec3 vertPos;\n"
		"in mediump vec3 vertNormal;\n"
		"out highp vec4 fragColor;\n"
		"layout (std140, binding=1) uniform Lights\n"
		"{\n"
		"   vec3 lightPos;\n"
		"   vec3 vertColor;\n"
		"};\n"
		"void main() {\n"
		"   highp vec3 L = normalize(lightPos - vertPos);\n"
		"   highp float NL = max(abs(dot(normalize(vertNormal), L)), 0.0);\n"
		"   highp vec3 color = vertColor;\n"
		"   highp vec3 col = clamp(color * (0.3 + 0.8*NL), 0.0, 1.0);\n"
		"   fragColor = vec4(col, 1.0);\n"
		"}\n";



//---------------------------------------------------------

void GLWidget::initializeGL()
{
	connect(context(), &QOpenGLContext::aboutToBeDestroyed,
			this, &GLWidget::cleanup);

	initializeOpenGLFunctions();
	glClearColor(0, 0, 0, m_transparent ? 0 : 1);

	glGenBuffers(1, &m_matricesUbo);
	glBindBuffer(GL_UNIFORM_BUFFER, m_matricesUbo);
	glBufferData(GL_UNIFORM_BUFFER, m_matUboSize, NULL, GL_STATIC_DRAW);
	glBindBuffer(GL_UNIFORM_BUFFER, 0);
	glBindBufferRange(GL_UNIFORM_BUFFER, 0, m_matricesUbo, 0, m_matUboSize);

	glGenBuffers(1, &m_lightsUbo);
	glBindBuffer(GL_UNIFORM_BUFFER, m_lightsUbo);
	glBufferData(GL_UNIFORM_BUFFER, m_lightUboSize, NULL, GL_STATIC_DRAW);
	glBindBuffer(GL_UNIFORM_BUFFER, 0);
	glBindBufferRange(GL_UNIFORM_BUFFER, 1, m_lightsUbo, 0, m_lightUboSize);

	glBindBuffer(GL_UNIFORM_BUFFER, m_lightsUbo);
	float lightPos[3] = {0.0f, 0.0f, 2.0f*m_modelSize};
	glBufferSubData(
				GL_UNIFORM_BUFFER, m_lightPosLoc, m_lightPosSiz, (void*)(lightPos));
	glBindBuffer(GL_UNIFORM_BUFFER, 0);

	bool success;
	m_shader = new QOpenGLShaderProgram;
	success = m_shader->addShaderFromSourceCode(
				QOpenGLShader::Vertex, vertexShaderSource);
	if ( !success ) qWarning() << m_shader->log();
	success = m_shader->addShaderFromSourceCode(
				QOpenGLShader::Fragment, fragmentShaderSource);
	if ( !success ) qWarning() << m_shader->log();
	m_shader->bindAttributeLocation("vertex", PROGRAM_VERTEX_ATTRIBUTE);
	m_shader->bindAttributeLocation("normal", PROGRAM_NORMAL_ATTRIBUTE);
	success = m_shader->link();
	if ( !success ) qWarning() << "Could not link shader program:" << m_shader->log();
//	m_shader->bind();

	m_shaderWireframe = new QOpenGLShaderProgram;
	success = m_shaderWireframe->addShaderFromSourceCode(
				QOpenGLShader::Vertex, vertexShaderSource);
	if ( !success ) qWarning() << m_shaderWireframe->log();
	success = m_shaderWireframe->addShaderFromSourceCode(
				QOpenGLShader::Geometry, geometryShaderSource);
	if ( !success ) qWarning() << m_shaderWireframe->log();
	success = m_shaderWireframe->addShaderFromSourceCode(
				QOpenGLShader::Fragment, fragmentShaderSourceWireframe);
	if ( !success ) qWarning() << m_shaderWireframe->log();
	m_shaderWireframe->bindAttributeLocation("vertex", PROGRAM_VERTEX_ATTRIBUTE);
	m_shaderWireframe->bindAttributeLocation("normal", PROGRAM_NORMAL_ATTRIBUTE);
	success = m_shaderWireframe->link();
	if ( !success ) qWarning() << "Could not link shader program:" << m_shaderWireframe->log();
	m_shaderWireframe->bind();
	m_shaderWireframe->setUniformValue("WIN_SCALE", QVector2D(this->width(),this->height()));

	m_meshVao.create();
	QOpenGLVertexArrayObject::Binder vaoBinderMesh(&m_meshVao);
	setGLMesh();

	m_meshNormsVao.create();
	QOpenGLVertexArrayObject::Binder vaoBinderMeshNorms(&m_meshNormsVao);
	setGLMeshNorms(1.0f);

#ifdef SHOW_MESH_DBUG
	m_meshDebugVao.create();
	QOpenGLVertexArrayObject::Binder vaoBinderMeshDebug(&m_meshDebugVao);
	setGLMeshDebug();
#endif

	setGLView();

	m_shaderWireframe->release();

}

//---------------------------------------------------------
void GLWidget::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
//	glEnable(GL_CULL_FACE);
//	glEnable(GL_PROGRAM_POINT_SIZE);
//#ifdef SHOW_MESH_DBUG
//	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//	glLineWidth(10.0f);
//#endif

	size_t nverts = m_mesh.vertCount();
	size_t ntris = m_mesh.triCount();

	m_world.rotate(-m_vRot / 8.0f, m_rotVect);
	m_camera.translate(-m_movVect / 2000.0f);
	QMatrix4x4 moveMatrix = m_camera * m_world.transposed();
	QMatrix3x3 normalMatrix = m_world.normalMatrix();

	glBindBuffer(GL_UNIFORM_BUFFER, m_matricesUbo);
	glBufferSubData(
				GL_UNIFORM_BUFFER, m_projMatLoc, m_projMatSiz,
				(void*)(m_proj.constData()));
	glBufferSubData(
				GL_UNIFORM_BUFFER, m_movMatLoc, m_movMatSiz,
				(void*)(moveMatrix.constData()));
	glBufferSubData(
				GL_UNIFORM_BUFFER, m_normMatLoc, m_normMatSiz,
				(void*)(normalMatrix.constData()));
	glBindBuffer(GL_UNIFORM_BUFFER, 0);

	glBindBuffer(GL_UNIFORM_BUFFER, m_lightsUbo);
	float lightPos[3] = {0.0f, 0.0f, 2.0f*m_modelSize};
	glBufferSubData(
				GL_UNIFORM_BUFFER, m_lightPosLoc, m_lightPosSiz,
				(void*)(lightPos));
	float meshColor[3] = {0.4f, 1.0f, 0.0f};
	glBufferSubData(
				GL_UNIFORM_BUFFER, m_colorLoc, m_colorSiz,
				(void*)(meshColor));

#ifdef SHOW_MESH
	m_shaderWireframe->bind();
	QOpenGLVertexArrayObject::Binder vaoBinder(&m_meshVao);
//	glDrawArrays(GL_POINTS, 0, nverts);
//	glDrawArrays(GL_TRIANGLES, 0, nverts);
	glDrawElements(GL_TRIANGLES, 3*ntris, GL_UNSIGNED_INT, 0);
	m_shaderWireframe->release();
#endif

	if ( m_showMeshNorms ){
		m_shader->bind();
		float normColor[3] = {1.0f, 0.0f, 1.0f};
		glBufferSubData(
					GL_UNIFORM_BUFFER, m_colorLoc, m_colorSiz,
					(void*)(normColor));
		QOpenGLVertexArrayObject::Binder vaoBinder3(&m_meshNormsVao);
		glDrawArrays(GL_LINES, 0, 2*nverts);
		m_shader->release();
	}

#ifdef SHOW_MESH_DBUG
	glDisable(GL_DEPTH_TEST);
	m_shader->bind();
	float debugColor[3] = {1.0f, 0.5f, 1.0f};
	glBufferSubData(
				GL_UNIFORM_BUFFER, m_colorLoc, m_colorSiz,
				(void*)(debugColor));
	QOpenGLVertexArrayObject::Binder vaoBinder4(&m_meshDebugVao);
	glDrawArrays(GL_LINES, 0, 2*m_mesh.debugCount());
	m_shader->release();
	glEnable(GL_DEPTH_TEST);
#endif

}
//---------------------------------------------------------

void GLWidget::resizeGL(int w, int h)
{
	m_aspectRatio = float(w)/h;
	m_proj.setToIdentity();
	m_proj.perspective(45.0f, m_aspectRatio, m_nearPlane, m_farPlane);
}

//---------------------------------------------------------

void GLWidget::setGLView()
{
	m_rotVect = QVector3D(0.0f,0.0f,0.0f);
	m_movVect = QVector3D(0.0f,0.0f,0.0f);

	m_nearPlane = 1e-5*m_modelSize;
	m_farPlane = 100.0f*m_modelSize;
	m_proj.setToIdentity();
	m_proj.perspective(45.0f, m_aspectRatio, m_nearPlane, m_farPlane);

	m_world.setToIdentity();
	m_camera.setToIdentity();
	m_camera.translate(0, 0, -m_modelSize);
}
//---------------------------------------------------------

void GLWidget::setGLMesh()
{
	size_t nverts = m_mesh.vertCount();
	size_t ntris = m_mesh.triCount();

	makeCurrent();

	if ( !m_meshVbo.isCreated() ) m_meshVbo.create();
	m_meshVbo.bind();
	if ( nverts > 0 ) {
		m_meshVbo.allocate(
					m_mesh.vertGLData(), 6*nverts*sizeof(GLfloat));
	}

	if ( !m_meshEbo.isCreated() ) m_meshEbo.create();
	m_meshEbo.bind();
	if ( ntris > 0 ) {
		m_meshEbo.allocate(
					m_mesh.triGLData(), 3*ntris*sizeof(GLuint));
	}
	setupVertexAttribs(m_meshVbo);

}

//---------------------------------------------------------

void GLWidget::setGLMeshNorms(float scale)
{
	size_t nverts = m_mesh.vertCount();

	makeCurrent();

	if ( !m_meshNormsVbo.isCreated() ) m_meshNormsVbo.create();
	m_meshNormsVbo.bind();
	if ( nverts > 0 ) {
		m_meshNormsVbo.allocate(
					m_mesh.normGLData(scale),
					12*nverts*sizeof(GLfloat));
	}
	setupVertexAttribs(m_meshNormsVbo);

}

//---------------------------------------------------------

void GLWidget::setGLMeshDebug()
{
	size_t ndebug = m_mesh.debugCount();

	makeCurrent();

	if ( !m_meshDebugVbo.isCreated() ) m_meshDebugVbo.create();
	m_meshDebugVbo.bind();
	if ( ndebug > 0 ) {
		m_meshDebugVbo.allocate(
					m_mesh.debugGLData(),
					12*ndebug*sizeof(GLfloat));
	}
	setupVertexAttribs(m_meshDebugVbo);

}

//---------------------------------------------------------

void GLWidget::setupVertexAttribs(QOpenGLBuffer vbo)
{

	vbo.bind();
	QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();

	f->glEnableVertexAttribArray(PROGRAM_VERTEX_ATTRIBUTE);
	f->glVertexAttribPointer(
				PROGRAM_VERTEX_ATTRIBUTE, 3, GL_FLOAT, GL_FALSE,
				6*sizeof(GLfloat), 0);

	f->glEnableVertexAttribArray(PROGRAM_NORMAL_ATTRIBUTE);
	f->glVertexAttribPointer(
				PROGRAM_NORMAL_ATTRIBUTE, 3, GL_FLOAT, GL_FALSE,
				6*sizeof(GLfloat), reinterpret_cast<void *>(3*sizeof(GLfloat)));
	vbo.release();


}

//---------------------------------------------------------

void GLWidget::mousePressEvent(QMouseEvent *event)
{
	m_lastMousePos = event->pos();
	//***
//	if (event->buttons() & Qt::MidButton){
//		setRandomMesh(35000);
//	}

}
//---------------------------------------------------------

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
	int dx = event->x() - m_lastMousePos.x();
	int dy = event->y() - m_lastMousePos.y();

	bool leftButton = (event->buttons() & Qt::LeftButton);
	bool rightButton = (event->buttons() & Qt::RightButton);

	if ( leftButton && rightButton ) {
		float norm = frobeniusNorm4x4(m_camera);
		setVectTranslation(QVector3D(-norm*dx,norm*dy,0.0f));
	}
	else if ( leftButton ) {
		int angle = abs(dx) + abs(dy);
		setVectRotation(angle, QVector3D(dy,dx,0.0f));
	}
	else if ( rightButton ) {
		int angle = abs(dx) + abs(dy);
		setVectRotation(angle, QVector3D(dy,0.0f,-dx));
	}

	m_lastMousePos = event->pos();
}
//---------------------------------------------------------

void GLWidget::wheelEvent(QWheelEvent *event)
{
	float norm = frobeniusNorm4x4(m_camera) + 1e-3;
	setVectTranslation(QVector3D(0.0f,0.0f,norm*event->delta()));
}
//---------------------------------------------------------

void GLWidget::setMesh(MeshPtr mesh)
{
	QMutexLocker locker(&m_recMutex);

	//m_mesh.fromPCL(mesh);

	updateMesh();
	setGLView();
}


//---------------------------------------------------------

void GLWidget::setRandomMesh(size_t ndim)
{
//	QMutexLocker locker(&m_recMutex);

	Eigen::VectorXf C = Eigen::VectorXf::Random(5);
	Eigen::VectorXf D = Eigen::VectorXf::Random(8);
	Eigen::VectorXf E = Eigen::VectorXf::Random(5);
	Eigen::VectorXf F = Eigen::VectorXf::Random(8);

	auto heightFun = [=](float xu, float xv){
		float pu = 15*xu, pv = 15*xv;
		float height =
				C(0)*cos(D(0)*pu) +
				C(1)*cos(D(1)*pv) +
				C(2)*cos(D(2)*pu)*cos(D(3)*pu) +
				C(3)*cos(D(4)*pu)*cos(D(5)*pv) +
				C(4)*cos(D(6)*pv)*cos(D(7)*pv) +
				E(0)*sin(F(0)*pu) +
				E(1)*sin(F(1)*pv) +
				E(2)*sin(F(2)*pu)*sin(F(3)*pu) +
				E(3)*sin(F(4)*pu)*sin(F(5)*pv) +
				E(4)*sin(F(6)*pv)*sin(F(7)*pv);
		return 0.1f*height;
	};

	QVector3D norm(0.0f,1.0f,0.0f);
	emit meshGenerate(norm, ndim, heightFun);

//	m_mesh.fromPlaneHeight(norm[0], norm[1], norm[2], ndim, heightFun);
//	updateMesh(0);
//	setGLView();

}

//---------------------------------------------------------

void GLWidget::getMesh(MeshPtr& mesh)
{
	QMutexLocker locker(&m_recMutex);

	//m_mesh.toPCL(mesh);
}

//---------------------------------------------------------

void GLWidget::viewGLMeshNorms(bool enabled)
{
	m_showMeshNorms = enabled;
	setGLMeshNorms(2.5e-2*m_modelSize);
	update();
}

//---------------------------------------------------------

void GLWidget::updateMesh(int opt)
{
	QMutexLocker locker(&m_recMutex);

	setGLMesh();
	if ( m_showMeshNorms ) setGLMeshNorms(1e-2*m_modelSize);

#ifdef SHOW_MESH_DBUG
	setGLMeshDebug();
#endif

	update();

	switch (opt) {
	case 1:
		setGLView();
		break;
	}
}

//---------------------------------------------------------

float GLWidget::frobeniusNorm4x4(QMatrix4x4 M)
{
	float normSq =
			M(0,0)*M(0,0) +
			M(1,0)*M(1,0) +
			M(2,0)*M(2,0) +
			M(3,0)*M(3,0) +
			M(0,1)*M(0,1) +
			M(1,1)*M(1,1) +
			M(2,1)*M(2,1) +
			M(3,1)*M(3,1) +
			M(0,2)*M(0,2) +
			M(1,2)*M(1,2) +
			M(2,2)*M(2,2) +
			M(3,2)*M(3,2) +
			M(0,3)*M(0,3) +
			M(1,3)*M(1,3) +
			M(2,3)*M(2,3) +
			M(3,3)*M(3,3);
	return sqrt(normSq);
}
//---------------------------------------------------------
