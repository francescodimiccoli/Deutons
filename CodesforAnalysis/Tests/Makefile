#define compiler
Compiler = g++
FLAGS = -std=c++11 -g
#exe
exe = System_test.exe
#file to compile
CC = System_test.cpp

#root libs
rootlibs=`root-config --libs --cflags`

all:
	$(Compiler) $(FLAGS) -o $(exe) $(CC) $(rootlibs) -lTMVA -I../include/ 
