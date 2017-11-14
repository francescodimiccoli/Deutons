#define compiler
Compiler = g++
FLAGS = -std=c++11 -g

#exe
exe = CountsExtraction_Parallel
exe1 = HeliumContamination_Parallel 
exe2 = MCEfficiency_Parallel
exe3 = Fluxes_Parallel
exe4 = EffCorr_Parallel
#file to compile
CC = sources/CountsExtraction_Parallel.cpp
CC1 = sources/HeliumContamination_Parallel.cpp
CC2 = sources/MCEfficiency_Parallel.cpp
CC3 = sources/Fluxes_Parallel.cpp
CC4 = sources/EffCorr_Parallel.cpp
GCC = include/GCR_data.cpp

#root libs
rootlibs=`root-config --libs --cflags` 
all:
	$(Compiler) $(FLAGS)  -o $(exe)  $(CC)  $(GCC) $(rootlibs) -lTMVA -Iinclude
	$(Compiler) $(FLAGS)  -o $(exe1) $(CC1) $(GCC) $(rootlibs) -lTMVA -Iinclude
	$(Compiler) $(FLAGS)  -o $(exe2) $(CC2) $(GCC) $(rootlibs) -lTMVA -Iinclude
	$(Compiler) $(FLAGS)  -o $(exe3) $(CC3) $(GCC) $(rootlibs) -lTMVA -Iinclude
	$(Compiler) $(FLAGS)  -o $(exe4) $(CC4) $(GCC) $(rootlibs) -lTMVA -Iinclude

