all: Analysis

DBARPATH=/afs/cern.ch/work/f/fdimicco/private/Deutons/data_reduction

obj/%.o: src/%.cpp
	g++ -std=c++11 -o $@ -c -I$(DBARPATH)/include/ -Iinclude $(shell root-config  --cflags) $^                  

obj/%.o: src-main/%.cpp 
	g++ -std=c++11 -o $@ -c -I$(DBARPATH)/include/ -Iinclude $(shell root-config  --cflags) $^                  

Analysis: obj/GCR_data.o obj/filesaver.o obj/Cuts.o obj/particle.o obj/binning.o obj/Livetime.o  obj/GlobalBinning.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/LatReweighter.o obj/InputFileReader.o obj/RangeMerger.o  obj/BetaSmearing.o obj/Efficiency.o obj/TemplateFITbetasmear.o obj/FluxRebin.o obj/effcorr.o  obj/rundb.o obj/Analyzer.o obj/CountExtraction.o  obj/Flux.o obj/AcceptanceAnalysis.o obj/Fluxes.o obj/Acceptance.o obj/TFractionFitter.o obj/Analysis.o 
	g++ $^ $(DBARPATH)/obj/RichBDT.o -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp -lgsl -lgslcblas 

LatitudeReweighter: obj/GCR_data.o obj/filesaver.o obj/Cuts.o obj/particle.o obj/binning.o obj/Livetime.o  obj/GlobalBinning.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/LatReweighter.o obj/InputFileReader.o obj/RangeMerger.o  obj/BetaSmearing.o obj/Efficiency.o obj/TemplateFITbetasmear.o obj/FluxRebin.o obj/effcorr.o  obj/rundb.o obj/Analyzer.o obj/CountExtraction.o  obj/Flux.o obj/AcceptanceAnalysis.o obj/Fluxes.o obj/Acceptance.o obj/TFractionFitter.o  obj/LatitudeReweighter.o
	 g++ $^ -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp

HeliumContamination_Parallel: obj/GCR_data.o obj/filesaver.o obj/Cuts.o obj/particle.o obj/binning.o obj/Livetime.o  obj/GlobalBinning.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/LatReweighter.o obj/InputFileReader.o obj/RangeMerger.o obj/BetaSmearing.o obj/Efficiency.o obj/TemplateFITbetasmear.o obj/FluxRebin.o obj/effcorr.o  obj/rundb.o obj/Analyzer.o obj/CountExtraction.o  obj/AcceptanceAnalysis.o obj/Acceptance.o obj/TFractionFitter.o obj/LatReweighter.o obj/HeliumContamination_Parallel.o
	g++ $^ $(DBARPATH)/obj/RichBDT.o -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp

Distributions_Plotter:obj/GCR_data.o obj/filesaver.o obj/Cuts.o obj/Livetime.o obj/particle.o obj/binning.o obj/GlobalBinning.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/LatReweighter.o obj/InputFileReader.o obj/RangeMerger.o obj/BetaSmearing.o obj/Efficiency.o obj/FluxRebin.o obj/TemplateFITbetasmear.o obj/effcorr.o  obj/rundb.o obj/Analyzer.o obj/CountExtraction.o  obj/AcceptanceAnalysis.o obj/Acceptance.o obj/HistoBooker.o obj/Distributions_Plotter.o
	g++ $^ $(DBARPATH)/obj/RichBDT.o -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp

Ntuple_Maker: obj/GCR_data.o obj/filesaver.o obj/Cuts.o obj/Livetime.o obj/GlobalBinning.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/Ntuple_Maker.o
	g++ $^ -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp

FastCheckSmearing: obj/GCR_data.o obj/filesaver.o obj/Cuts.o obj/particle.o obj/binning.o obj/Livetime.o obj/GlobalBinning.o obj/DBarReader.o obj/Variables.o obj/RangeMerger.o obj/InputFileReader.o obj/reweight.o obj/FastCheckSmearing.o
	g++ $^ -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp

ExtractInfos: obj/GCR_data.o obj/filesaver.o obj/Cuts.o obj/Livetime.o obj/GlobalBinning.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/ExtractInfos.o
	g++ $^ -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp

Resolution_Calibrations: obj/GCR_data.o obj/filesaver.o obj/Cuts.o obj/Livetime.o obj/GlobalBinning.o obj/DBarReader.o obj/Variables.o obj/reweight.o obj/Resolution_Calibrations.o
	g++ $^ -o $@ $(shell root-config --libs --cflags) -lTMVA -L$(DBARPATH)/lib/ -lntp


clean:
	rm -f obj/*

