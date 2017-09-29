#define compiler
Compiler = g++
FLAGS = -std=c++11 -g
#exe
exe = System_test.exe
exe2 = MCFineTuning
exe3 = PrettyPlotting
exe4 = EfficiencyTests
#file to compile
CC = System_test.cpp
CC2 = MCFineTuning.cpp
CC3 = PrettyPlotting.cpp
CC4 = EfficiencyTests.cpp
#root libs
rootlibs=`root-config --libs --cflags`

all:
	#$(Compiler) $(FLAGS) -o $(exe) $(CC) $(rootlibs) -lTMVA -I../include/ 
	$(Compiler) $(FLAGS) -o $(exe2) $(CC2) $(rootlibs) -lTMVA -I../include/
	$(Compiler) $(FLAGS) -o $(exe3) $(CC3) $(rootlibs) -lTMVA -I../include/	
	$(Compiler) $(FLAGS) -o $(exe4) $(CC4) $(rootlibs) -lTMVA -I../include/	

