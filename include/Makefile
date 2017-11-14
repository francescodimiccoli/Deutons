#define compiler
Compiler = g++
FLAGS = -std=c++11 -g

#exe
exe = Cuts
#file to compile
CC = Cuts.cpp
#root libs
rootlibs=`root-config --libs --cflags` 
all:
	$(Compiler) $(FLAGS)  -o $(exe) $(CC)  $(rootlibs) -lTMVA
