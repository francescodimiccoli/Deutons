all: dstmake

CASTORLIBS=/opt/exp_software/ams/additional_libs/lib64

%.o: %.cpp
	g++ -g -c $^ -I$(AMSWD)/include $(shell root-config --cflags --libs)

dstmake: InputFileReader.o selezioni.o dstmake.o HitMCTruth.o
	LD_LIBRARY_PATH=$(CASTORLIBS):$(LD_LIBRARY_PATH) g++ -o $@ $^ $(AMSLIBso) $(shell root-config --cflags --libs)

clean:
	rm -f InputFileReader.o selezioni.o mainMCtree_new.o maindatatree_new.o HitMCTruth.o
