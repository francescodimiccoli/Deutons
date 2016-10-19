all: maindatatree

CASTORLIBS=/opt/exp_software/ams/additional_libs/lib64

%.o: %.C
	g++ -c $^ -I$(AMSWD)/include $(shell root-config --cflags --libs)

mainMCtree: selezioni.o mainMCtree_new.o
	LD_LIBRARY_PATH=$(CASTORLIBS):$(LD_LIBRARY_PATH) g++ -o $@ $^ $(AMSLIBso) $(shell root-config --cflags --libs)

maindatatree: selezioni.o maindatatree_new.o
	LD_LIBRARY_PATH=$(CASTORLIBS):$(LD_LIBRARY_PATH) g++ -o $@ $^ $(AMSLIBso) $(shell root-config --cflags --libs)

clean:
	rm -f selezioni.o mainMCtree_new.o maindatatree_new.o
