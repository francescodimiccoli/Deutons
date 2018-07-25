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
exe6= P_Fluxes

cpp = ../src/Variables.cpp
cpp2 = ../src/reweight.cpp
cpp3 = ../src/GCR_data.cpp
cpp4 = ../src/GlobalBinning.cpp
cpp5 = ../src/filesaver.cpp
cpp6 = ../src/Cuts.cpp
cpp7 = ../src/LatReweighter.cpp
cpp8 = ../src/DBarReader.cpp
cpp9 = ../src/Livetime.cpp

#file to compile
CC = CountsExtraction.cpp
CC1= HeliumContamination.cpp
CC2= MCEfficiency.cpp
CC3= D_Fluxes.cpp
CC4= RICHEffCorr.cpp
CC6= P_Fluxes.cpp


DBARPATH=../../data_reduction


rootlibs=`root-config --libs --cflags` 
all:
	$(Compiler) $(FLAGS)  -o $(exe) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6)  $(cpp7)  $(cpp8)  $(cpp9) $(CC)  $(rootlibs) -lTMVA -I../include -L$(DBARPATH)/lib/ -lntp 
	$(Compiler) $(FLAGS)  -o $(exe1) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6)  $(cpp7)  $(cpp8)  $(cpp9) $(CC1)  $(rootlibs) -lTMVA -I../include -L$(DBARPATH)/lib/ -lntp 
	#$(Compiler) $(FLAGS)  -o $(exe2) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6)  $(cpp7)  $(cpp8)  $(cpp9) $(CC2)  $(rootlibs) -lTMVA -I../include  -L$(DBARPATH)/lib/ -lntp
	$(Compiler) $(FLAGS)  -o $(exe3) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6)  $(cpp7)  $(cpp8)  $(cpp9) $(CC3)  $(rootlibs) -lTMVA -I../include  -L$(DBARPATH)/lib/ -lntp
	#$(Compiler) $(FLAGS)  -o $(exe4) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6)  $(cpp7)  $(cpp8)  $(cpp9) $(CC4)  $(rootlibs) -lTMVA -I../include  -L$(DBARPATH)/lib/ -lntp
	$(Compiler) $(FLAGS)  -o $(exe6) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6)  $(cpp7)  $(cpp8)  $(cpp9) $(CC6)  $(rootlibs) -lTMVA -I../include  -L$(DBARPATH)/lib/ -lntp
