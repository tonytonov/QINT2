######################################################################
# Automatically generated by qmake (3.0) Fri May 30 12:48:26 2014
######################################################################

TEMPLATE = app
TARGET = QINT2
INCLUDEPATH += .

# Input
HEADERS += \
    rqmcintegrator.h \
    qintintegrator.h \
    interceptableintegrand.h \
    utils.h
SOURCES += main.cpp \
    rqmcintegrator.cpp \
    qintintegrator.cpp

QMAKE_CXXFLAGS += -std=c++0x

LIBS += -L/usr/ -lhintlib
