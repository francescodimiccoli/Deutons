#define compiler
Compiler = g++ -std=c++11 -g 

#exe
exedue = Analysis
exequint = TimeAnalysis.exe
exe = prova
#file to compile
CCdue = Analysis.cpp
CCquint = TimeAnalysis.cpp
CC = prova.C
#root libs
rootlibs=`root-config --libs --cflags` 
all:
	$(Compiler) -o $(exedue) $(CCdue) $(rootlibs) -I../include/
	$(Compiler) -o $(exequint) $(CCquint) $(rootlibs) -I../include/
	$(Compiler) -o $(exe) $(CC) $(rootlibs) -I../include/

