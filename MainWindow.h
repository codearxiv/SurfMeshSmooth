//     Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//     See LICENSE included with this distribution.

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "constants.h"

#include <QMainWindow>

QT_BEGIN_NAMESPACE
class QAction;
class QMenu;
class QPlainTextEdit;
class QListWidget;
QT_END_NAMESPACE

class GLWidget;
class Window;
class MessageLogger;
class RandomSurfDialog;
class NormalsDialog;
class NoiseDialog;
class NoiseDialog;
class SmoothDialog;
class OptionsDialog;

class MainWindow : public QMainWindow
{
	Q_OBJECT

	typedef void* MeshPtr;

public:
	explicit MainWindow(QWidget *parent = nullptr);
	void badInputMessageBox(const QString& info);

private:
	void open();
	void saveAs();
	void viewGLNorms(bool enabled);
	void setRandom();
	void approxNorms();
	void noise();
	void smooth();
    void options();
    void about();

public slots:
	void appendLogText(const QString& text);
	void insertLogText(const QString& text);
	void replaceLogText(const QString& text);

signals:
	void meshChanged(MeshPtr mesh);
	void meshQueried(MeshPtr& mesh);
	void meshNormsViewGL(bool enabled);
	void meshSetRandom(size_t nDim);
	void meshApproxNorms();
	void meshNoise(int nSweeps);
	void meshSmooth(int nSweeps);
	void normScaleChanged(float scale);

private:
	Window *centralWidget;
	QToolBar *mainToolBar;
	QStatusBar *statusBar;
	QToolBar *toolBar;
	QPlainTextEdit *logText;
    RandomSurfDialog *randomSurfDialog;
	NormalsDialog *normalsDialog;
	NoiseDialog *noiseDialog;
	SmoothDialog *smoothDialog;
    OptionsDialog *optionsDialog;

	MessageLogger *msgLogger;

};

#endif // MAINWINDOW_H
