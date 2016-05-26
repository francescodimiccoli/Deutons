#define compiler
Compiler = g++ -std=c++11

#exe
exedue = Analysis
exequint = TimeAnalysis.exe

#file to compile
CCdue = Analysis.cpp
CCquint = TimeAnalysis.cpp

#root libs
rootlibs=`root-config --libs --cflags` 
all:
	$(Compiler) -o $(exedue) $(CCdue) $(rootlibs) -lTMVA 
	$(Compiler) -o $(exequint) $(CCquint) $(rootlibs) -lTMVA


