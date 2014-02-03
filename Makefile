#define compiler
Compiler = g++
#exe
exe = AnalyzeMC.exe
exean = AnalyzeDATA.exe
#file to compile
CC = AnalyzeMC.C
CCan = AnalyzeDATA.C
#root libs
rootlibs=`root-config --libs --cflags`
all:
	$(Compiler) -o $(exe) $(CC) $(rootlibs) 
	#$(Compiler) -o $(exean) $(CCan) $(rootlibs)
clean:
	rm $(exe) $(exean)

