//     Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//     See LICENSE included with this distribution.


#include "MainWindow.h"
#include "Window.h"
#include "MessageLogger.h"
#include "RandomSurfDialog.h"
#include "NormalsDialog.h"
#include "NoiseDialog.h"
#include "SmoothDialog.h"
#include "OptionsDialog.h"
#include "constants.h"

#include <QApplication>
#include <QCoreApplication>
#include <QToolBar>
#include <QDockWidget>
#include <QPlainTextEdit>
#include <QScrollBar>
#include <QFileDialog>
#include <QMenuBar>
#include <QMenu>
#include <QMessageBox>
#include <QMainWindow>


MainWindow::MainWindow(QWidget *parent) :
	QMainWindow(parent)
{
	if (objectName().isEmpty())
		setObjectName(QString::fromUtf8("MainWindow"));
	resize(1204, 640);

	setWindowTitle(QCoreApplication::translate("MainWindow", "SurfMeshSmooth", nullptr));

	//------
	// Add docks
	QDockWidget *dock = new QDockWidget("Log Window", this);
	dock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
	logText = new QPlainTextEdit;
	logText->setReadOnly(true);
	dock->setWidget(logText);
	addDockWidget(Qt::RightDockWidgetArea, dock);

	//------
	// Add actions
	QMenu *fileMenu = menuBar()->addMenu("&File");
	QMenu *editMenu = menuBar()->addMenu("&Edit");
	QMenu *viewMenu = menuBar()->addMenu("&View");
	QMenu *toolsMenu = menuBar()->addMenu("&Tools");
	QMenu *helpMenu = menuBar()->addMenu("&Help");


	QToolBar *fileToolBar = addToolBar("File");

	const QIcon openIcon =
			QIcon::fromTheme("document-open", QIcon(":/images/open.png"));
	QAction *openAct = new QAction(openIcon, "&Open...", this);
	openAct->setShortcuts(QKeySequence::Open);
	openAct->setStatusTip("Open an existing PCD file");
	connect(openAct, &QAction::triggered, this, &MainWindow::open);
	fileMenu->addAction(openAct);
	fileToolBar->addAction(openAct);

	const QIcon saveAsIcon =
			QIcon::fromTheme("document-save-as", QIcon(":/images/save.png"));
	QAction *saveAsAct = new QAction(saveAsIcon, "Save &As...", this);
	saveAsAct->setShortcuts(QKeySequence::SaveAs);
	saveAsAct->setStatusTip("Save PCD to disk");
	connect(saveAsAct, &QAction::triggered, this, &MainWindow::saveAs);
	fileMenu->addAction(saveAsAct);
	fileToolBar->addAction(saveAsAct);

	const QIcon exitIcon = QIcon::fromTheme("application-exit");
	QAction *exitAct =
			fileMenu->addAction(exitIcon, "E&xit", this, &QWidget::close);
	exitAct->setShortcuts(QKeySequence::Quit);
	exitAct->setStatusTip("Exit SurfMeshSmooth");


	QToolBar *editToolBar = addToolBar("Edit");

	viewMenu->addAction(dock->toggleViewAction());

	QAction *viewNormsAct = new QAction("&View normals", this);
	viewNormsAct->setStatusTip("View mesh vertex normals.");
	viewNormsAct->setCheckable(true);
	connect(viewNormsAct, &QAction::triggered, this, &MainWindow::viewGLNorms);
	viewMenu->addAction(viewNormsAct);

    QAction *randomSurfAct = new QAction("&Random surface", this);
    randomSurfAct->setStatusTip("Create a mesh from a random surface");
	connect(randomSurfAct, &QAction::triggered, this, &MainWindow::setRandom);
    toolsMenu->addAction(randomSurfAct);

	QAction *normalsAct = new QAction("&Approx. Normals", this);
	normalsAct->setStatusTip("Approximate mesh surface normals");
	connect(normalsAct, &QAction::triggered, this, &MainWindow::approxNorms);
	toolsMenu->addAction(normalsAct);

	QAction *noiseAct = new QAction("&Noise", this);
	noiseAct->setStatusTip("Add noise a random subset of the mesh vertices");
	connect(noiseAct, &QAction::triggered, this, &MainWindow::noise);
	toolsMenu->addAction(noiseAct);

	QAction *smoothAct = new QAction("&Smooth", this);
	smoothAct->setStatusTip("Smooth mesh");
	connect(smoothAct, &QAction::triggered, this, &MainWindow::smooth);
	toolsMenu->addAction(smoothAct);

    toolsMenu->addSeparator();

    QAction *optionsAct = new QAction("&Options", this);
    optionsAct->setStatusTip("Change app settings.");
	connect(optionsAct, &QAction::triggered, this, &MainWindow::options);
    toolsMenu->addAction(optionsAct);


	QAction *aboutAct =
			helpMenu->addAction("&About", this, &MainWindow::about);
	aboutAct->setStatusTip("About");

	QAction *aboutQtAct =
			helpMenu->addAction("About &Qt", qApp, &QApplication::aboutQt);
	aboutQtAct->setStatusTip("About Qt");


	//------
	// Central widget

	msgLogger = new MessageLogger(logText);
	centralWidget = new Window(this, msgLogger);
	centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
	setCentralWidget(centralWidget);

	connect(msgLogger, &MessageLogger::logTextAppend,
			this, &MainWindow::appendLogText);

	connect(msgLogger, &MessageLogger::logTextInsert,
			this, &MainWindow::insertLogText);

	connect(msgLogger, &MessageLogger::logTextReplace,
			this, &MainWindow::replaceLogText);

	connect(this, &MainWindow::meshChanged,
			centralWidget, &Window::setMesh);

	connect(this, &MainWindow::meshQueried,
			centralWidget, &Window::getMesh);

	connect(this, &MainWindow::meshNormsViewGL,
			centralWidget, &Window::viewGLMeshNorms);

	connect(this, &MainWindow::meshSetRandom,
			centralWidget, &Window::setRandomMesh);

	connect(this, &MainWindow::meshApproxNorms,
			centralWidget, &Window::approxMeshNorms);

	connect(this, &MainWindow::meshNoise,
			centralWidget, &Window::noiseMesh);

	connect(this, &MainWindow::meshSmooth,
			centralWidget, &Window::smoothMesh);

	connect(this, &MainWindow::normScaleChanged,
			centralWidget, &Window::setNormScale);

	//------
	// Dialogs

    randomSurfDialog = new RandomSurfDialog(this);
	normalsDialog = new NormalsDialog(this);
	noiseDialog = new NoiseDialog(this);
	smoothDialog = new SmoothDialog(this);
    optionsDialog = new OptionsDialog(this);


	//------

	QMetaObject::connectSlotsByName(this);

}

//---------------------------------------------------------

void MainWindow::open()
{

//	emit meshChanged(mesh);

}

//---------------------------------------------------------

void MainWindow::saveAs()
{

}

//---------------------------------------------------------

void MainWindow::viewGLNorms(bool enabled)
{
	emit meshNormsViewGL(enabled);
}

//---------------------------------------------------------

void MainWindow::setRandom()
{
	// Show the dialog as modal
	if (randomSurfDialog->exec() == QDialog::Accepted){
		size_t nDim;
		bool ok = randomSurfDialog->getFields(nDim);
		if (!ok){
			badInputMessageBox("All fields should be integers bigger than zero.");
			return;
		}
		emit meshSetRandom(nDim);
	}

}


//---------------------------------------------------------

void MainWindow::approxNorms()
{
	// Show the dialog as modal
	if (normalsDialog->exec() == QDialog::Accepted){
		emit meshApproxNorms();
	}

}

//---------------------------------------------------------

void MainWindow::noise()
{
	// Show the dialog as modal
	if (noiseDialog->exec() == QDialog::Accepted){
		int nSweeps;
		int ok = noiseDialog->getFields(nSweeps);
		switch ( ok ){
		case -1:
			badInputMessageBox("Number of noising sweeps.");
			break;
		default:
			emit meshNoise(nSweeps);
		}
	}

}


//---------------------------------------------------------

void MainWindow::smooth()
{
	// Show the dialog as modal
	if (smoothDialog->exec() == QDialog::Accepted){
		int nSweeps;

		int ok = smoothDialog->getFields(nSweeps);
		switch ( ok ){
		case -1:
			badInputMessageBox(
						"Number of iterations field must be bigger than zero.");
			break;
		default:
			emit meshSmooth(nSweeps);
		}
	}

}


//---------------------------------------------------------

void MainWindow::options()
{
    // Show the dialog as modal
	if (optionsDialog->exec() == QDialog::Accepted){
		float normScale;
		bool ok = optionsDialog->getFields(normScale);
		if (!ok){
			badInputMessageBox("All field should be greater than zero.");
            return;
        }
		emit normScaleChanged(normScale);
	}

}

//---------------------------------------------------------

void MainWindow::about()
{
   QMessageBox::about(this, "About",
					  QString("SurfMeshSmooth.\n") +
					  QString("Copyright 2021 Piotr (Peter) Beben.\n") +
					  QString("pdbcas2@gmail.com\n"));
}

//---------------------------------------------------------

void MainWindow::badInputMessageBox(const QString& info)
{
	QMessageBox msgBox;
	msgBox.setText("Invalid input.");
	msgBox.setInformativeText(info);
	msgBox.setStandardButtons(QMessageBox::Ok);
	msgBox.setDefaultButton(QMessageBox::Ok);
	msgBox.setIcon(QMessageBox::Information);
	msgBox.exec();
}

//---------------------------------------------------------

void MainWindow::appendLogText(const QString& text)
{
	logText->appendPlainText(text);
	logText->verticalScrollBar()->setValue(
				logText->verticalScrollBar()->maximum());
}

//---------------------------------------------------------

void MainWindow::insertLogText(const QString& text)
{
	logText->insertPlainText(text);
	logText->verticalScrollBar()->setValue(
				logText->verticalScrollBar()->maximum());
}
//---------------------------------------------------------

void MainWindow::replaceLogText(const QString& text)
{
	logText->undo();
	logText->appendPlainText(text);
	logText->verticalScrollBar()->setValue(
				logText->verticalScrollBar()->maximum());
	logText->moveCursor(QTextCursor::End, QTextCursor::MoveAnchor);
}

//---------------------------------------------------------
