# QT += opengl xml
QT += core gui widgets
# CONFIG += qt

TARGET = renderer
TEMPLATE = app
CONFIG += console

# TODO mac os x only
QMAKE_CXX = /usr/local/Cellar/llvm/4.0.0_1/bin/clang++
# QMAKE_CXX = g++-8
# QMAKE_CXXFLAGS -= -stdlib=libc++
QMAKE_CXXFLAGS += -msse2
QMAKE_CXXFLAGS += -std=c++14
CONFIG += c++14

# allow debug and release modes, and set -DNDEBUG for release mode
CONFIG += debug_and_release
CONFIG(release, debug|release): DEFINES += NDEBUG

Release:DESTDIR = release
Release:OBJECTS_DIR = release/.obj
Release:MOC_DIR = release/.moc
Release:RCC_DIR = release/.rcc
Release:UI_DIR = release/.ui

Debug:DESTDIR = debug
Debug:OBJECTS_DIR = debug/.obj
Debug:MOC_DIR = debug/.moc
Debug:RCC_DIR = debug/.rcc
Debug:UI_DIR = debug/.ui

win32 {
    DEFINES += GLEW_STATIC
    LIBS += -lopengl32 -lglu32
}

linux {
    # fixes an issue with linux and vr libraries
    QMAKE_CXXFLAGS += -D_GLIBCXX_USE_CXX11_ABI=0
}

HEADERS += \
    libs/tiny_obj_loader/tiny_obj_loader.h \
    src/BSDF.h \
    src/camera.h \
    src/common.h \
    src/distance_estimation.h \
    src/geom.h \
    src/kdtree.h \
    src/mesh.h \
    src/raster.h \
    src/rendering.h \
    src/renderthreads.h \
    src/sampling.h \
    src/scene.h \
    src/ui/mainwindow.h \
    src/ui/view.h

SOURCES += \
    libs/tiny_obj_loader/tiny_obj_loader.cpp \
    src/BSDF.cpp \
    src/camera.cpp \
    src/common.cpp \
    src/distance_estimation.cpp \
    src/geom.cpp \
    src/kdtree.cpp \
    src/main.cpp \
    src/mesh.cpp \
    src/raster.cpp \
    src/rendering.cpp \
    src/renderthreads.cpp \
    src/sampling.cpp \
    src/scene.cpp \
    src/ui/mainwindow.cpp \
    src/ui/view.cpp

FORMS += src/ui/mainwindow.ui
INCLUDEPATH += libs src src/ui libs/glm libs/tiny_obj_loader # libs/glew-1.10.0/include
DEPENDPATH += libs src src/ui libs/glm libs/tiny_obj_loader # libs/glew-1.10.0/include

DEFINES += _USE_MATH_DEFINES
# DEFINES += TIXML_USE_STL
DEFINES +=  GLM_FORCE_RADIANS GLM_ENABLE_EXPERIMENTAL
# OTHER_FILES += \
#     shaders/*

# RCC_DIR = shaders/
# RESOURCES += \
#     resources.qrc
