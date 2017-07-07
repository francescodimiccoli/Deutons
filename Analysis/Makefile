#define compiler
Compiler = g++
FLAGS = -std=c++11 -g

#exe
exe = CountsExtraction
exe1 = HeliumContamination 
exe2 = MCEfficiency
exe3 = Fluxes
exe4 = EffCorr
#file to compile
CC = sources/CountsExtraction.cpp
CC1 = sources/HeliumContamination.cpp
CC2 = sources/MCEfficiency.cpp
CC3 = sources/Fluxes.cpp
CC4 = sources/EffCorr.cpp
#root libs
rootlibs=`root-config --libs --cflags` 
all:
	$(Compiler) $(FLAGS)  -o $(exe) $(CC)  $(rootlibs) -lTMVA
	$(Compiler) $(FLAGS)  -o $(exe1) $(CC1)  $(rootlibs) -lTMVA
	$(Compiler) $(FLAGS)  -o $(exe2) $(CC2)  $(rootlibs) -lTMVA
	$(Compiler) $(FLAGS)  -o $(exe3) $(CC3)  $(rootlibs) -lTMVA
	$(Compiler) $(FLAGS)  -o $(exe4) $(CC4)  $(rootlibs) -lTMVA
