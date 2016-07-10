#GTEST = /storage/gpfs_ams/ams/users/kostya/googletest-release-1.7.0/
GTEST = /home/kanishev/Programming/googletest-release-1.7.0/
CXX = g++
CXXFLAGS = -std=gnu++0x -g -L$(GTEST) -lgtest -lpthread 
INCS = -I./ -I../ -I$(GTEST)/include
OBJS = testBinning.o testHistogram.o testHistUtils.o
#ROOTLIBS = `root-config --libs --cflags`

runtest: testAll
	LD_LIBRARY_PATH=$(GTEST):$(LD_LIBRARY_PATH) ./testAll

testAll: $(OBJS)
	$(CXX) testAll.cpp ../include/GCR_data.cpp $(CXXFLAGS) $(INCS) $(OBJS) -o testAll

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCS) $(ROOTLIBS)

clean:
	rm testAll *.o testAll.xml
