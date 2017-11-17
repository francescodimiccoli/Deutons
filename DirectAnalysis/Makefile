all: CountsExtraction_Parallel EffCorr_Parallel

DBARPATH=../../dbar_ntp

obj/%.o: src/%.cpp
	g++ -std=c++11 -o $@ -c -I$(DBARPATH)/include/ -Iinclude $(shell root-config --libs --cflags) $^                  

obj/%.o: src-main/%.cpp 
	g++ -std=c++11 -o $@ -c -I$(DBARPATH)/include/ -Iinclude $(shell root-config --libs --cflags) $^                  

CountsExtraction_Parallel: obj/GCR_data.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/CountsExtraction_Parallel.o
	g++ $^ -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp

EffCorr_Parallel: obj/GCR_data.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/EffCorr_Parallel.o
	g++ $^ -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp



clean:
	rm -f obj/*

