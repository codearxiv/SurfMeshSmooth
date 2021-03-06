//     Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//     See LICENSE included with this distribution.

#include "GLWidget.h"
#include "MainWindow.h"
#include "constants.h"

#include <Eigen/Dense>
#include <QApplication>
#include <QDesktopWidget>
#include <QSurfaceFormat>
#include <QCommandLineParser>
#include <QCommandLineOption>
#include <QTextBlock>
#include <QTextCursor>
#include <cstdlib>
#include <iostream>

void out_of_mem_handler();

int main(int argc, char *argv[])
{
	srand(2);
	std::set_new_handler(out_of_mem_handler);
	Eigen::initParallel();

	qRegisterMetaType<size_t>("size_t");
	qRegisterMetaType<QTextBlock>("QTextBlock");
	qRegisterMetaType<QTextCursor>("QTextCursor");
	qRegisterMetaType<std::function<float(float xu, float xv)>>(
				"std::function<float(float xu, float xv)>");

	QCoreApplication::setApplicationName("SurfMeshSmooth");
	QCoreApplication::setOrganizationName("Piotr Beben");
	QCoreApplication::setApplicationVersion(QT_VERSION_STR);
	QCoreApplication::addLibraryPath("./");

	QApplication app(argc, argv);

	QCommandLineParser parser;
	parser.setApplicationDescription(QCoreApplication::applicationName());
	parser.addHelpOption();
	parser.addVersionOption();
	QCommandLineOption multipleSampleOption("multisample", "Multisampling");
	parser.addOption(multipleSampleOption);
	QCommandLineOption coreProfileOption("coreprofile", "Use core profile");
	parser.addOption(coreProfileOption);
	QCommandLineOption transparentOption("transparent", "Transparent window");
	parser.addOption(transparentOption);
	parser.process(app);

	QSurfaceFormat fmt;
	fmt.setDepthBufferSize(24);
	if (parser.isSet(multipleSampleOption))
		fmt.setSamples(4);
	if (parser.isSet(coreProfileOption)) {
		fmt.setVersion(3, 2);
		fmt.setProfile(QSurfaceFormat::CoreProfile);
	}
	QSurfaceFormat::setDefaultFormat(fmt);

	MainWindow mainWindow;

	GLWidget::setTransparent(parser.isSet(transparentOption));
	if (GLWidget::isTransparent()) {
		mainWindow.setAttribute(Qt::WA_TranslucentBackground);
		mainWindow.setAttribute(Qt::WA_NoSystemBackground, false);
	}
	mainWindow.resize(mainWindow.sizeHint());
//	int desktopArea = QApplication::desktop()->width() *
//					 QApplication::desktop()->height();
//	int widgetArea = mainWindow.width() * mainWindow.height();
//	if (((float)widgetArea / (float)desktopArea) < 0.75f)
//		mainWindow.show();
//	else
//		mainWindow.showMaximized();

	mainWindow.showMaximized();

	return app.exec();

}


void out_of_mem_handler()
{
	std::cerr << "Unable to allocate memory.\n";
	std::abort();
}
