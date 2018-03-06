all: HeliumContamination_Parallel CountsExtraction_Parallel MCEfficiency_Parallel EffCorr_Parallel Distributions_Plotter Ntuple_Maker Resolution_Calibrations

DBARPATH=../data_reduction

obj/%.o: src/%.cpp
	g++ -std=c++11 -o $@ -c -I$(DBARPATH)/include/ -Iinclude $(shell root-config --libs --cflags) $^                  

obj/%.o: src-main/%.cpp 
	g++ -std=c++11 -o $@ -c -I$(DBARPATH)/include/ -Iinclude $(shell root-config --libs --cflags) $^                  

HeliumContamination_Parallel: obj/GCR_data.o obj/filesaver.o obj/Cuts.o obj/Livetime.o obj/GlobalBinning.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/HeliumContamination_Parallel.o
	g++ $^ -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp

CountsExtraction_Parallel: obj/GCR_data.o obj/filesaver.o obj/Cuts.o obj/Livetime.o obj/GlobalBinning.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/CountsExtraction_Parallel.o
	g++ $^ -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp

MCEfficiency_Parallel: obj/GCR_data.o obj/filesaver.o obj/Cuts.o obj/Livetime.o obj/GlobalBinning.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/MCEfficiency_Parallel.o
	g++ $^ -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp

EffCorr_Parallel: obj/GCR_data.o obj/filesaver.o obj/Cuts.o obj/Livetime.o obj/GlobalBinning.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/EffCorr_Parallel.o
	g++ $^ -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp

Distributions_Plotter: obj/GCR_data.o obj/filesaver.o obj/Cuts.o obj/Livetime.o obj/GlobalBinning.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/HistoBooker.o obj/Distributions_Plotter.o
	g++ $^ -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp

Ntuple_Maker: obj/GCR_data.o obj/filesaver.o obj/Cuts.o obj/Livetime.o obj/GlobalBinning.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/Ntuple_Maker.o
	g++ $^ -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp

Resolution_Calibrations: obj/GCR_data.o obj/filesaver.o obj/Cuts.o obj/Livetime.o obj/GlobalBinning.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/Resolution_Calibrations.o
	g++ $^ -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp

clean:
	rm -f obj/*

