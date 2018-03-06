#define compiler
Compiler = g++
FLAGS = -std=c++11 -g

#exe
exe = CountsExtraction
exe1= HeliumContamination
#exe2= P_Fluxes
exe2= MCEfficiency
exe3= D_Fluxes
exe4= RICHEffCorr 
exe5= QualEffCorr 

cpp = ../src/Variables.cpp
cpp2 = ../src/reweight.cpp
cpp3 = ../src/GCR_data.cpp
cpp4 = ../src/GlobalBinning.cpp
cpp5 = ../src/filesaver.cpp
cpp6 = ../src/Cuts.cpp


#file to compile
CC = CountsExtraction.cpp
CC1= HeliumContamination.cpp
CC2= MCEfficiency.cpp
CC3= D_Fluxes.cpp
CC4= RICHEffCorr.cpp
CC5= QualEffCorr.cpp



#
#root libs
rootlibs=`root-config --libs --cflags` 
all:
	$(Compiler) $(FLAGS)  -o $(exe) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6) $(CC)  $(rootlibs) -lTMVA -I../include
	$(Compiler) $(FLAGS)  -o $(exe1) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6) $(CC1)  $(rootlibs) -lTMVA -I../include
	$(Compiler) $(FLAGS)  -o $(exe2) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6) $(CC2)  $(rootlibs) -lTMVA -I../include

