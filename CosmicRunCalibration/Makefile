CC = g++
CFLAGS = -Wall -c -g


INCLUDES = -I. -I$(ROOTSYS)/include 

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib  -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit


all: calibrateFibresWithCosmics FitTools.o

calibrateFibresWithCosmics: calibrateFibresWithCosmics.cpp FitTools.o
	$(CC) -Wall $(INCLUDES) -o calibrateFibresWithCosmics calibrateFibresWithCosmics.cpp $(ROOTFLAG) $(EXTRALIBS)

FitTools.o: ../PositionAnalysis/src/FitTools.cc
	$(CC) $(CFLAGS) -c ../PositionAnalysis/src/FitTools.cc $(ROOTFLAG) $(EXTRALIBS)






clean:
	rm *.o
