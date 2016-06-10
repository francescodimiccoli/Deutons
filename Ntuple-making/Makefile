#define compiler
Compiler = g++
#exe
exe = AnalyzeMC.exe
exean = AnalyzeDATA.exe
#file to compile
CC = AnalyzeMC5D.C
CCan = AnalyzeDATA5D.C

#root libs
rootlibs=`root-config --libs --cflags`

all:
	$(Compiler) -o $(exe) $(CC)  $(rootlibs) -lTMVA
	$(Compiler) -o $(exean) $(CCan) $(rootlibs) -lTMVA

#  g++ -g -std=c++11 $(root-config --libs --cflags)  -lTMVA -o AnalyzeMC.exe AnalyzeMC5D.C
# include/GCR_data.cpp  -Wno-pointer-arith  -fmax-errors=5
