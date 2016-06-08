#define compiler
Compiler = g++
#exe
exe = Preliminar.exe
#file to compile
CC = Preliminar.C

#root libs
rootlibs=`root-config --libs --cflags` 
all:
	$(Compiler) -o $(exe) $(CC)  $(rootlibs) -lTMVA
