#define compiler
Compiler = g++
FLAGS = -std=c++11 -g

#exe
exe = CountsExtraction
exe1= MCEfficiency
exe2= P_Fluxes
exe3= D_Fluxes
exe4= RICHEffCorr 
exe5= QualEffCorr 




#file to compile
CC = CountsExtraction.cpp
CC1= MCEfficiency.cpp
CC2= P_Fluxes.cpp
CC3= D_Fluxes.cpp
CC4= RICHEffCorr.cpp
CC5= QualEffCorr.cpp



#
#root libs
rootlibs=`root-config --libs --cflags` 
all:
	$(Compiler) $(FLAGS)  -o $(exe) $(CC)  $(rootlibs) -lTMVA
	$(Compiler) $(FLAGS)  -o $(exe1) $(CC1)  $(rootlibs) -lTMVA
	$(Compiler) $(FLAGS)  -o $(exe2) $(CC2)  $(rootlibs) -lTMVA
	$(Compiler) $(FLAGS)  -o $(exe3) $(CC3)  $(rootlibs) -lTMVA
	$(Compiler) $(FLAGS)  -o $(exe4) $(CC4)  $(rootlibs) -lTMVA
	$(Compiler) $(FLAGS)  -o $(exe5) $(CC5)  $(rootlibs) -lTMVA

