#-------------------------------------------------
#
# Project created by QtCreator 2014-09-11T20:15:29
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = junction_simulator
TEMPLATE = app


SOURCES += main.cpp\
    rootfinder.cpp \
    cubic.cpp \
    junctionModeler.cpp

HEADERS  += \
    givensQR.hpp \
    matrix.hpp \
    rootfinder.h \
    cubic.h \
    persistence1d.hpp \
    junctionModeler.h

FORMS    += mainwindow.ui

INCLUDEPATH += $$PWD/C:/Qt/Tools/mingw482_32/i686-w64-mingw32/include
DEPENDPATH += $$PWD/C:/Qt/Tools/mingw482_32/i686-w64-mingw32/include

CONFIG += qwt
CONFIG += console

QMAKE_CXXFLAGS += -std=c++11

RESOURCES += \
    resources.qrc
