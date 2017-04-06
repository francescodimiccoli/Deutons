#define compiler
Compiler = g++
FLAGS = -std=c++11 -g

#exe
exe = Resolution_Calibrations
exeplot = PrettyPlotting
#file to compile
CC = Resolution_Calibrations.cpp
CCplot = PrettyPlotting.cpp
#root libs
rootlibs=`root-config --libs --cflags` 
all:
	$(Compiler) $(FLAGS)  -o $(exe) $(CC)  $(rootlibs) -lTMVA
	$(Compiler) $(FLAGS)  -o $(exeplot) $(CCplot)  $(rootlibs) -lTMVA
