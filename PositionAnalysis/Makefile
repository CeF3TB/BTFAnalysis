CC = g++
CFLAGS = -Wall -c -g


INCLUDES = -I. -I$(ROOTSYS)/include 

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib  -lHtml -lMathCore -lGenVector -lMinuit -lTMVA




all: positionAnalysis drawPositionAnalysis hodoscopeAnalysis calibrateCeF3 drawHodoClustersVsEnergy drawShitVsDiagonal calibrateBGO drawEnergyScan calibrateCaloPositioning drawCaloPositioning calibrateCeF3Positioning makePedestalFile makeAnalysisTree drawPedestalsVsRun drawResolutionStudies calibrateBGOPositioning checkBGOPositioning calibrateCef3Positioning drawCef3Positioning drawAngularStudies singleElectronSelection drawCalibComparison drawCeF3Intercalibration make2DPlots



checkRuns: checkRuns.cpp
	$(CC) -Wall $(INCLUDES) -o checkRuns checkRuns.cpp $(ROOTFLAG) $(EXTRALIBS)

makePedestalFile: makePedestalFile.cpp 
	$(CC) -Wall $(INCLUDES) -o makePedestalFile makePedestalFile.cpp $(ROOTFLAG) $(EXTRALIBS)

makeAnalysisTree: makeAnalysisTree.cpp RunHelper.o HodoCluster.o EnergyCalibration.o CalibrationUtility.o 
	$(CC) -Wall $(INCLUDES) -o makeAnalysisTree makeAnalysisTree.cpp RunHelper.o HodoCluster.o EnergyCalibration.o CalibrationUtility.o $(ROOTFLAG) $(EXTRALIBS)

positionAnalysis: HodoCluster.o RunHelper.o PositionTools.o positionAnalysis.cpp hodo_efficiency.dat
	$(CC) -Wall $(INCLUDES) -o positionAnalysis positionAnalysis.cpp HodoCluster.o RunHelper.o PositionTools.o $(ROOTFLAG) $(EXTRALIBS)

make2DPlots: HodoCluster.o RunHelper.o PositionTools.o DrawTools.o make2DPlots.cpp hodo_efficiency.dat
	$(CC) -Wall $(INCLUDES) -o make2DPlots make2DPlots.cpp HodoCluster.o RunHelper.o DrawTools.o PositionTools.o $(ROOTFLAG) $(EXTRALIBS)


singleElectronSelection: HodoCluster.o RunHelper.o PositionTools.o singleElectronSelection.cpp hodo_efficiency.dat
	$(CC) -Wall $(INCLUDES) -o singleElectronSelection singleElectronSelection.cpp HodoCluster.o RunHelper.o PositionTools.o $(ROOTFLAG) $(EXTRALIBS)

drawPositionAnalysis: DrawTools.o RunHelper.o PositionTools.o drawPositionAnalysis.cpp 
	$(CC) -Wall $(INCLUDES) -o drawPositionAnalysis drawPositionAnalysis.cpp DrawTools.o RunHelper.o PositionTools.o $(ROOTFLAG) $(EXTRALIBS)

drawPedestalsVsRun: DrawTools.o drawPedestalsVsRun.cpp 
	$(CC) -Wall $(INCLUDES) -o drawPedestalsVsRun drawPedestalsVsRun.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

calibrateCeF3: calibrateCeF3.cpp DrawTools.o
	$(CC) -Wall $(INCLUDES) -o calibrateCeF3 calibrateCeF3.cpp DrawTools.o  $(ROOTFLAG) $(EXTRALIBS)

calibrateBGO: DrawTools.o calibrateBGO.cpp
	$(CC) -Wall $(INCLUDES) -o calibrateBGO calibrateBGO.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

calibrateCaloPositioning: PositionTools.o calibrateCaloPositioning.cpp
	$(CC) -Wall $(INCLUDES) -o calibrateCaloPositioning calibrateCaloPositioning.cpp PositionTools.o $(ROOTFLAG) $(EXTRALIBS)

calibrateBGOPositioning: PositionTools.o DrawTools.o calibrateBGOPositioning.cpp
	$(CC) -Wall $(INCLUDES) -o calibrateBGOPositioning calibrateBGOPositioning.cpp DrawTools.o PositionTools.o $(ROOTFLAG) $(EXTRALIBS)

calibrateCef3Positioning: PositionTools.o DrawTools.o calibrateCef3Positioning.cpp
	$(CC) -Wall $(INCLUDES) -o calibrateCef3Positioning calibrateCef3Positioning.cpp DrawTools.o PositionTools.o $(ROOTFLAG) $(EXTRALIBS)

hodoscopeAnalysis: hodoscopeAnalysis.cpp
	$(CC) -Wall $(INCLUDES) -o hodoscopeAnalysis hodoscopeAnalysis.cpp $(ROOTFLAG) $(EXTRALIBS)

drawHodoClustersVsEnergy: DrawTools.o drawHodoClustersVsEnergy.cpp
	$(CC) -Wall $(INCLUDES) -o drawHodoClustersVsEnergy drawHodoClustersVsEnergy.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

drawShitVsDiagonal: drawShitVsDiagonal.cpp DrawTools.o
	$(CC) -Wall $(INCLUDES) -o drawShitVsDiagonal drawShitVsDiagonal.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

drawEnergyScan: drawEnergyScan.cpp DrawTools.o FitTools.o
	$(CC) -Wall $(INCLUDES) -o drawEnergyScan drawEnergyScan.cpp DrawTools.o FitTools.o $(ROOTFLAG) $(EXTRALIBS)

drawCaloPositioning: drawCaloPositioning.cpp DrawTools.o RunHelper.o
	$(CC) -Wall $(INCLUDES) -o drawCaloPositioning drawCaloPositioning.cpp DrawTools.o RunHelper.o $(ROOTFLAG) $(EXTRALIBS)

drawCef3Positioning: drawCef3Positioning.cpp DrawTools.o RunHelper.o
	$(CC) -Wall $(INCLUDES) -o drawCef3Positioning drawCef3Positioning.cpp DrawTools.o RunHelper.o $(ROOTFLAG) $(EXTRALIBS)

drawResolutionStudies: drawResolutionStudies.cpp DrawTools.o FitTools.o
	$(CC) -Wall $(INCLUDES) -o drawResolutionStudies drawResolutionStudies.cpp DrawTools.o FitTools.o $(ROOTFLAG) $(EXTRALIBS)

checkBGOPositioning: checkBGOPositioning.cpp DrawTools.o
	$(CC) -Wall $(INCLUDES) -o checkBGOPositioning checkBGOPositioning.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

drawCeF3Intercalibration: drawCeF3Intercalibration.cpp DrawTools.o
	$(CC) -Wall $(INCLUDES) -o drawCeF3Intercalibration drawCeF3Intercalibration.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

drawCalibComparison: drawCalibComparison.cpp DrawTools.o
	$(CC) -Wall $(INCLUDES) -o drawCalibComparison drawCalibComparison.cpp DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

drawAngularStudies: drawAngularStudies.cpp DrawTools.o FitTools.o
	$(CC) -Wall $(INCLUDES) -o drawAngularStudies drawAngularStudies.cpp DrawTools.o FitTools.o $(ROOTFLAG) $(EXTRALIBS)




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

FitTools.o: src/FitTools.cc
	$(CC) $(CFLAGS) -c src/FitTools.cc $(ROOTFLAG) $(EXTRALIBS)

PositionTools.o: src/PositionTools.cc
	$(CC) $(CFLAGS) -c src/PositionTools.cc $(ROOTFLAG) $(EXTRALIBS)




clean:
	rm *.o
	rm positionAnalysis
	rm make2DPlots	
	rm singleElectronSelection
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
	rm drawCeF3Intercalibration
