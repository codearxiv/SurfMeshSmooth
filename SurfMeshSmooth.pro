#-------------------------------------------------
#
# Project created by QtCreator 2021-09-28T10:50:52
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = SurfMeshSmooth
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS
DEFINES += EIGEN_NO_MALLOC
DEFINES += EIGEN_NO_DEBUG

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11

PRECOMPILED_HEADER  = stable.h

SOURCES += \
	GLWidget.cpp \
	MainWindow.cpp \
	Mesh.cpp \
	MeshWorker.cpp \
	MessageLogger.cpp \
	Window.cpp \
	dialogs/NoiseDialog.cpp \
	dialogs/NormalsDialog.cpp \
	dialogs/OptionsDialog.cpp \
	dialogs/RandomSurfDialog.cpp \
	dialogs/SmoothDialog.cpp \
	dialogs/get_field.cpp \
	main.cpp \
	utils/Plane.cpp \
	utils/rotations.cpp \
	utils/timer.cpp \
	utils/vect_utils.cpp

HEADERS += \
	GLWidget.h \
	MainWindow.h \
	Mesh.h \
	MeshWorker.h \
	MessageLogger.h \
	PairHash.h \
	Window.h \
	alignment.h \
	constants.h \
	dialogs/NoiseDialog.h \
	dialogs/NormalsDialog.h \
	dialogs/OptionsDialog.h \
	dialogs/RandomSurfDialog.h \
	dialogs/SmoothDialog.h \
	dialogs/get_field.h \
	stable.h \
	utils/Plane.h \
	utils/ensure_buffer_size.h \
	utils/rotations.h \
	utils/selection_sort.h \
	utils/timer.h \
	utils/vect_utils.h

FORMS +=

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

INCLUDEPATH += $$PWD/'../../lib/eigen3' \
		   $$PWD/'../../lib/Boost/include/boost-1_68' \
		   $$PWD/'../../lib/OpenNI2/Include' \
		   $$PWD/'../../lib/FLANN/include' \
		   $$PWD/'dialogs' \
		   $$PWD/'utils'

LIBS += -L$$PWD/'../../lib/Boost/lib/'
LIBS += -L$$PWD/'../../lib/OpenNI2/Lib/' -lOpenNI2

RESOURCES += \
	SurfMeshSmooth.qrc

msvc{
    QMAKE_CXXFLAGS += -openmp
}

gcc{
	QMAKE_CXXFLAGS += -fopenmp -Wno-attributes -Wno-unused-parameter
	QMAKE_CXXFLAGS += -Wno-ignored-attributes
    QMAKE_LFLAGS += -fopenmp
    LIBS += -fopenmp
}
