#define compiler
Compiler = g++
FLAGS = -std=c++11 -g
#exe
exe = AnalyzeMC.exe
exean = AnalyzeDATA.exe
#file to compile
CC = AnalyzeMC5D.C
CCan = AnalyzeDATA5D.C

GCC = ../include/GCR_data.cpp
#root libs
rootlibs=`root-config --libs --cflags`

all:
	$(Compiler) $(FLAGS) -o $(exe) $(CC) $(GCC)  $(rootlibs) -lTMVA -I../include/ 
	$(Compiler) $(FLAGS) -o $(exean) $(CCan) $(GCC) $(rootlibs) -lTMVA -I../include/

#  g++ -g -std=c++11 $(root-config --libs --cflags)  -lTMVA -o AnalyzeMC.exe AnalyzeMC5D.C
# include/GCR_data.cpp  -Wno-pointer-arith  -fmax-errors=5
#
# g++ -g -std=c++11 $(root-config --libs --cflags)  -lTMVA -o AnalyzeMC.exe AnalyzeMC5D.C ../include/GCR_data.cpp  -I../include/ -Wno-pointer-arith
