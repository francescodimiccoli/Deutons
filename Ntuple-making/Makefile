#define compiler
Compiler = g++
FLAGS = -std=c++11 -g
#exe
exe = Analyzer.exe
#file to compile
CC = Analyzer.cpp

GCC = ../include/GCR_data.cpp
#root libs
rootlibs=`root-config --libs --cflags`

all:
	$(Compiler) $(FLAGS) -o $(exe) $(CC) $(GCC)  $(rootlibs) -lTMVA -I../include/ 

#  g++ -g -std=c++11 $(root-config --libs --cflags)  -lTMVA -o AnalyzeMC.exe AnalyzeMC5D.C
# include/GCR_data.cpp  -Wno-pointer-arith  -fmax-errors=5
#
# g++ -g -std=c++11 $(root-config --libs --cflags)  -lTMVA -o AnalyzeMC.exe AnalyzeMC5D.C ../include/GCR_data.cpp  -I../include/ -Wno-pointer-arith
#
# g++ -g -std=c++11 $(root-config --libs --cflags)  -lTMVA -o AnalyzeDATA.exe AnalyzeDATA5D.C ../include/GCR_data.cpp  -I../include/ -Wno-pointer-arith
# 
