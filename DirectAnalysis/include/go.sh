g++ -g -std=c++11 -L/afs/cern.ch/exp/ams/Offline/root/Linux/root-v5-34-9-gcc64-slc6/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -pthread -m64 -I/afs/cern.ch/exp/ams/Offline/root/Linux/root-v5-34-9-gcc64-slc6/include  -o GCRtest GCRtest.cpp GCR_data.cpp -fmax-errors=10 -Wno-pointer-arith