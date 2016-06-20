#define compiler
Compiler = g++ -std=c++11 -g 

#exe
exedue = Analysis
exequint = TimeAnalysis.exe

#file to compile
CCdue = Analysis.cpp
CCquint = TimeAnalysis.cpp

#root libs
rootlibs=`root-config --libs --cflags` 
all:
	$(Compiler) -o $(exedue) $(CCdue) $(rootlibs) -I../include/
	$(Compiler) -o $(exequint) $(CCquint) $(rootlibs) -I../include/
