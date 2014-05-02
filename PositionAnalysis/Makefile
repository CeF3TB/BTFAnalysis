CC = g++
CFLAGS = -Wall -c -g


INCLUDES = -I. -I$(ROOTSYS)/include 

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib  -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit



all: positionAnalysis drawPositionAnalysis hodoscopeAnalysis calibrateCeF3 drawHodoClustersVsEnergy drawShitVsDiagonal


positionAnalysis: HodoCluster.o positionAnalysis.cpp hodo_efficiency.dat
	$(CC) -Wall $(INCLUDES) -o positionAnalysis positionAnalysis.cpp HodoCluster.o $(ROOTFLAG) $(EXTRALIBS)

drawPositionAnalysis: DrawTools.o drawPositionAnalysis.cpp
	$(CC) -Wall $(INCLUDES) -o drawPositionAnalysis drawPositionAnalysis.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

calibrateCeF3: calibrateCeF3.cpp
	$(CC) -Wall $(INCLUDES) -o calibrateCeF3 calibrateCeF3.cpp $(ROOTFLAG) $(EXTRALIBS)

hodoscopeAnalysis: hodoscopeAnalysis.cpp
	$(CC) -Wall $(INCLUDES) -o hodoscopeAnalysis hodoscopeAnalysis.cpp $(ROOTFLAG) $(EXTRALIBS)

drawHodoClustersVsEnergy: DrawTools.o drawHodoClustersVsEnergy.cpp
	$(CC) -Wall $(INCLUDES) -o drawHodoClustersVsEnergy drawHodoClustersVsEnergy.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

drawShitVsDiagonal: drawShitVsDiagonal.cpp
	$(CC) -Wall $(INCLUDES) -o drawShitVsDiagonal drawShitVsDiagonal.cpp $(ROOTFLAG) $(EXTRALIBS)




HodoCluster.o: src/HodoCluster.cc
	$(CC) $(CFLAGS) -c src/HodoCluster.cc $(ROOTFLAG) $(EXTRALIBS)

DrawTools.o: src/DrawTools.cc
	$(CC) $(CFLAGS) -c src/DrawTools.cc $(ROOTFLAG) $(EXTRALIBS)



clean:
	rm positionAnalysis
	rm drawPositionAnalysis 
	rm hodoscopeAnalysis
	rm calibrateCeF3
	rm drawHodoClustersVsEnergy
	rm drawShitVsDiagonal
	rm *.o
