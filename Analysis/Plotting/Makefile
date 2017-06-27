#define compiler
Compiler = g++
FLAGS = -std=c++11 -g

#exe
exe = CountsExtraction
exe1= MCEfficiency
exe2= P_Fluxes
exe3= D_Fluxes


#file to compile
CC = CountsExtraction.cpp
CC1= MCEfficiency.cpp
CC2= P_Fluxes.cpp
CC3= D_Fluxes.cpp
#
#root libs
rootlibs=`root-config --libs --cflags` 
all:
	$(Compiler) $(FLAGS)  -o $(exe) $(CC)  $(rootlibs) -lTMVA
	$(Compiler) $(FLAGS)  -o $(exe1) $(CC1)  $(rootlibs) -lTMVA
	$(Compiler) $(FLAGS)  -o $(exe2) $(CC2)  $(rootlibs) -lTMVA
	$(Compiler) $(FLAGS)  -o $(exe3) $(CC3)  $(rootlibs) -lTMVA

