#define compiler
Compiler = g++
FLAGS = -std=c++11 -g

#exe
exe = CountsExtraction
exe1= MCEfficiency
#file to compile
CC = CountsExtraction.cpp
CC1= MCEfficiency.cpp
#root libs
rootlibs=`root-config --libs --cflags` 
all:
	$(Compiler) $(FLAGS)  -o $(exe) $(CC)  $(rootlibs) -lTMVA
	$(Compiler) $(FLAGS)  -o $(exe1) $(CC1)  $(rootlibs) -lTMVA
