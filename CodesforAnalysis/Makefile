#define compiler
Compiler = g++ -std=c++11
#exe
exedue = Analysis
exetre = QualityLKLHOOD.exe
exesex = QualityLKLHOOD_NaF.exe
exesev = QualityLKLHOOD_Agl.exe
#exequad = Preliminar.exe
exequint = TimeAnalysis.exe
#file to compile
CCdue = Analysis.cpp
CCtre = QualityLKLHOOD.cpp
CCsex = QualityLKLHOOD_NaF.cpp
CCsev = QualityLKLHOOD_Agl.cpp
#CCquad = Preliminar.cpp
CCquint = TimeAnalysis.cpp

#root libs
rootlibs=`root-config --libs --cflags` 
all:
	$(Compiler) -o $(exedue) $(CCdue) $(rootlibs) -lTMVA 
	$(Compiler) -o $(exetre) $(CCtre) $(rootlibs) -lTMVA
	$(Compiler) -o $(exesex) $(CCsex) $(rootlibs) -lTMVA
	$(Compiler) -o $(exesev) $(CCsev) $(rootlibs) -lTMVA
#	$(Compiler) -o $(exequad) $(CCquad) $(rootlibs) -lTMVA
	$(Compiler) -o $(exequint) $(CCquint) $(rootlibs) -lTMVA


