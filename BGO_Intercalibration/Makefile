## source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.12/x86_64-slc5-gcc46-opt/root/bin/thisroot.sh

CC = g++
CFLAGS = -Wall -c -g


INCLUDES = -I. -I$(ROOTSYS)/include 

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib  -lHtml -lMathCore -lGenVector -lMinuit -lTMVA -lRooFit -lRooFitCore



all: bgoIntercalibration


bgoIntercalibration: bgoIntercalibration.cc
	$(CC) -Wall $(INCLUDES) -o bgoIntercalibration bgoIntercalibration.cc $(ROOTFLAG) $(EXTRALIBS)

clean:
	rm *.o
	rm bgoIntercalibration

