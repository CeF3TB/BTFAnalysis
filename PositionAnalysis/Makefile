CC = g++
CFLAGS = -Wall -c -g


INCLUDES = -I. -I$(ROOTSYS)/include 

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib  -lHtml -lMathCore -lGenVector -lMinuit -lTMVA



all: positionAnalysis drawPositionAnalysis hodoscopeAnalysis calibrateCeF3 drawHodoClustersVsEnergy drawShitVsDiagonal calibrateBGO drawEnergyScan calibrateCaloPositioning drawCaloPositioning calibrateCeF3Positioning makePedestalFile makeAnalysisTree drawPedestalsVsRun


checkRuns: checkRuns.cpp
	$(CC) -Wall $(INCLUDES) -o checkRuns checkRuns.cpp $(ROOTFLAG) $(EXTRALIBS)

makePedestalFile: makePedestalFile.cpp 
	$(CC) -Wall $(INCLUDES) -o makePedestalFile makePedestalFile.cpp $(ROOTFLAG) $(EXTRALIBS)

makeAnalysisTree: makeAnalysisTree.cpp RunHelper.o HodoCluster.o EnergyCalibration.o CalibrationUtility.o 
	$(CC) -Wall $(INCLUDES) -o makeAnalysisTree makeAnalysisTree.cpp RunHelper.o HodoCluster.o EnergyCalibration.o CalibrationUtility.o $(ROOTFLAG) $(EXTRALIBS)

positionAnalysis: HodoCluster.o RunHelper.o positionAnalysis.cpp hodo_efficiency.dat
	$(CC) -Wall $(INCLUDES) -o positionAnalysis positionAnalysis.cpp HodoCluster.o RunHelper.o $(ROOTFLAG) $(EXTRALIBS)

drawPositionAnalysis: DrawTools.o RunHelper.o drawPositionAnalysis.cpp 
	$(CC) -Wall $(INCLUDES) -o drawPositionAnalysis drawPositionAnalysis.cpp DrawTools.o RunHelper.o $(ROOTFLAG) $(EXTRALIBS)

drawPedestalsVsRun: DrawTools.o drawPedestalsVsRun.cpp 
	$(CC) -Wall $(INCLUDES) -o drawPedestalsVsRun drawPedestalsVsRun.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

calibrateCeF3: calibrateCeF3.cpp DrawTools.o
	$(CC) -Wall $(INCLUDES) -o calibrateCeF3 calibrateCeF3.cpp DrawTools.o  $(ROOTFLAG) $(EXTRALIBS)

calibrateBGO: DrawTools.o calibrateBGO.cpp
	$(CC) -Wall $(INCLUDES) -o calibrateBGO calibrateBGO.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

calibrateCaloPositioning: RunHelper.o calibrateCaloPositioning.cpp
	$(CC) -Wall $(INCLUDES) -o calibrateCaloPositioning calibrateCaloPositioning.cpp RunHelper.o $(ROOTFLAG) $(EXTRALIBS)

calibrateCeF3Positioning: RunHelper.o DrawTools.o calibrateCeF3Positioning.cpp
	$(CC) -Wall $(INCLUDES) -o calibrateCeF3Positioning calibrateCeF3Positioning.cpp RunHelper.o DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

hodoscopeAnalysis: hodoscopeAnalysis.cpp
	$(CC) -Wall $(INCLUDES) -o hodoscopeAnalysis hodoscopeAnalysis.cpp $(ROOTFLAG) $(EXTRALIBS)

drawHodoClustersVsEnergy: DrawTools.o drawHodoClustersVsEnergy.cpp
	$(CC) -Wall $(INCLUDES) -o drawHodoClustersVsEnergy drawHodoClustersVsEnergy.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

drawShitVsDiagonal: drawShitVsDiagonal.cpp DrawTools.o
	$(CC) -Wall $(INCLUDES) -o drawShitVsDiagonal drawShitVsDiagonal.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

drawEnergyScan: drawEnergyScan.cpp DrawTools.o
	$(CC) -Wall $(INCLUDES) -o drawEnergyScan drawEnergyScan.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

drawCaloPositioning: drawCaloPositioning.cpp DrawTools.o RunHelper.o
	$(CC) -Wall $(INCLUDES) -o drawCaloPositioning drawCaloPositioning.cpp DrawTools.o RunHelper.o $(ROOTFLAG) $(EXTRALIBS)




HodoCluster.o: src/HodoCluster.cc
	$(CC) $(CFLAGS) -c src/HodoCluster.cc $(ROOTFLAG) $(EXTRALIBS)

DrawTools.o: src/DrawTools.cc
	$(CC) $(CFLAGS) -c src/DrawTools.cc $(ROOTFLAG) $(EXTRALIBS)

RunHelper.o: src/RunHelper.cc
	$(CC) $(CFLAGS) -c src/RunHelper.cc $(ROOTFLAG) $(EXTRALIBS)

CalibrationUtility.o: src/CalibrationUtility.cc
	$(CC) $(CFLAGS) -c src/CalibrationUtility.cc $(ROOTFLAG) $(EXTRALIBS)

EnergyCalibration.o: src/EnergyCalibration.cc
	$(CC) $(CFLAGS) -c src/EnergyCalibration.cc $(ROOTFLAG) $(EXTRALIBS)



clean:
	rm *.o
	rm positionAnalysis
	rm drawCaloPositioning
	rm drawPositionAnalysis 
	rm hodoscopeAnalysis
	rm calibrateCeF3
	rm calibrateBGO
	rm calibrateCaloPositioning
	rm drawHodoClustersVsEnergy
	rm drawShitVsDiagonal
	rm drawEnergyScan
	rm drawPedestalsVsRun
