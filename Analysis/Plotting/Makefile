#define compiler
Compiler = g++
FLAGS = -std=c++11 -g

#exe
exe = CountsExtraction
#file to compile
CC = CountsExtraction.cpp
#root libs
rootlibs=`root-config --libs --cflags` 
all:
	$(Compiler) $(FLAGS)  -o $(exe) $(CC)  $(rootlibs) -lTMVA
