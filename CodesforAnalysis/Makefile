#Compiler options
Compiler = g++  
FLAGS	:=-std=c++11  -Wfatal-errors

# Source files
EXEC	= Analysis QualityLKLHOOD QualityLKLHOOD_NaF QualityLKLHOOD_Agl Preliminar TimeAnalysis
CLASS = TemplateFITClass
BIN = $(EXEC:%=%.exe)

# Libs
rootlibs=`root-config --libs --cflags`


#$?  List of dependencies changed more recently than current target.
#$^  List of dependencies.
#$@  Name of current target.
#$<  Name of current dependency.
#$*  Name of current dependency without extension.


all: $(CLASS) $(BIN) 


$(BIN): $(EXEC)
	$(Compiler) $(FLAGS) -o $@  $(<:%=obj/%.o)  $(CLASS:%=obj/%.o) $(rootlibs) -lTMVA

$(EXEC):
	$(Compiler) $(FLAGS) -o obj/$@.o -c $@.cpp  $(rootlibs) -lTMVA

$(CLASS):
	$(Compiler) $(FLAGS) -o obj/$@.o -c Parte2/$@.cpp  $(rootlibs) -lTMVA

clean:
	rm -rf obj/*.o
