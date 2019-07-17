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
exe7 = TimeAnalysis
exe8 = PlotBinnedHisto

cpp = ../src/Variables.cpp
cpp2 = ../src/reweight.cpp
cpp3 = ../src/GCR_data.cpp
cpp4 = ../src/GlobalBinning.cpp
cpp5 = ../src/filesaver.cpp
cpp6 = ../src/Cuts.cpp
cpp7 = ../src/LatReweighter.cpp
cpp8 = ../src/DBarReader.cpp
cpp9 = ../src/Livetime.cpp
cpp10 = ../src/rundb.cpp
cpp11 = ../src/binning.cpp
cpp12 = ../src/particle.cpp
cpp13 = ../src/Efficiency.cpp
cpp14 = ../src/TemplateFITbetasmear.cpp
cpp15 = ../src/effcorr.cpp

#file to compile
CC = CountsExtraction.cpp
CC1= HeliumContamination.cpp
CC2= MCEfficiency.cpp
CC3= D_Fluxes.cpp
CC4= RICHEffCorr.cpp
CC6= P_Fluxes.cpp
CC7= TimeAnalysis.cpp
CC8= PlotBinnedHisto.cpp

DBARPATH=/data1/home/fdimicco/Deutons/data_reduction


rootlibs=`root-config --libs --cflags` 
all:
#	$(Compiler) $(FLAGS) $(DBARPATH)/obj/RichBDT.o  -o  $(exe)  $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6)  $(cpp7)  $(cpp8)  $(cpp9)  $(cpp10) $(cpp11) $(cpp12)  $(cpp13) $(cpp14) $(cpp15) $(CC)   $(rootlibs) -lTMVA -I../include -I$(DBARPATH)/include -L$(DBARPATH)/lib/ -lntp 
	$(Compiler) $(FLAGS) $(DBARPATH)/obj/RichBDT.o  -o  $(exe1) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6)  $(cpp7)  $(cpp8)  $(cpp9)  $(cpp10) $(cpp11) $(cpp12)  $(cpp13) $(cpp14) $(cpp15) $(CC1)  $(rootlibs) -lTMVA -I../include  -I$(DBARPATH)/include -L$(DBARPATH)/lib/ -lntp 
#	$(Compiler) $(FLAGS) $(DBARPATH)/obj/RichBDT.o  -o  $(exe2) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6)  $(cpp7)  $(cpp8)  $(cpp9)  $(cpp10) $(cpp11) $(cpp12)  $(cpp13) $(cpp14) $(cpp15) $(CC2)  $(rootlibs) -lTMVA -I../include -I$(DBARPATH)/include -L$(DBARPATH)/lib/ -lntp
#	$(Compiler) $(FLAGS) $(DBARPATH)/obj/RichBDT.o  -o  $(exe3) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6)  $(cpp7)  $(cpp8)  $(cpp9)  $(cpp10) $(cpp11) $(cpp12)  $(cpp13) $(cpp14) $(cpp15) $(CC3)  $(rootlibs) -lTMVA -I../include  -L$(DBARPATH)/lib/ -lntp
#	$(Compiler) $(FLAGS) $(DBARPATH)/obj/RichBDT.o  -o  $(exe4) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6)  $(cpp7)  $(cpp8)  $(cpp9)  $(cpp10) $(cpp11) $(cpp12)  $(cpp13) $(cpp14) $(cpp15) $(CC4)  $(rootlibs) -lTMVA -I../include -I$(DBARPATH)/include -L$(DBARPATH)/lib/ -lntp
#	$(Compiler) $(FLAGS) $(DBARPATH)/obj/RichBDT.o  -o  $(exe6) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6)  $(cpp7)  $(cpp8)  $(cpp9)  $(cpp10) $(cpp11) $(cpp12)  $(cpp13) $(cpp14) $(cpp15) $(CC6)  $(rootlibs) -lTMVA -I../include  -L$(DBARPATH)/lib/ -lntp
#	$(Compiler) $(FLAGS) $(DBARPATH)/obj/RichBDT.o  -o  $(exe7) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6)  $(cpp7)  $(cpp8)  $(cpp9)  $(cpp10) $(cpp11) $(cpp12)  $(cpp13) $(cpp14) $(cpp15) $(CC7)  $(rootlibs) -lTMVA -I../include -I$(DBARPATH)/include -L$(DBARPATH)/lib/ -lntp
#	$(Compiler) $(FLAGS) $(DBARPATH)/obj/RichBDT.o  -o  $(exe8) $(cpp) $(cpp2) $(cpp3) $(cpp4) $(cpp5) $(cpp6)  $(cpp7)  $(cpp8)  $(cpp9)  $(cpp10) $(cpp11) $(cpp12)  $(cpp13) $(cpp14) $(cpp15) $(CC8)  $(rootlibs) -lTMVA -I../include  -L$(DBARPATH)/lib/ -lntp
