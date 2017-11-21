export AMSWD=/cvmfs/ams.cern.ch/Offline/vdev
export AMSDataDir=/cvmfs/ams.cern.ch/Offline/AMSDataDir
source /cvmfs/ams.cern.ch/Offline/root/Linux/root-v5-34-9-gcc64-slc6/bin/thisroot.sh
export AMSLIBso=/cvmfs/ams.cern.ch/Offline/vdev/lib/linuxx8664gcc5.34/ntuple_slc6_PG.so
export LD_LIBRARY_PATH=/opt/exp_software/ams/additional_libs/lib64:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/afs/cern.ch/work/f/fdimicco/private/Deutons/data_reduction/lib:$LD_LIBRARY_PATH
#source /storage/gpfs_ams/ams/users/fdimicco/Deutons/Setup/x86_64-slc6/setup.sh
source /cvmfs/sft.cern.ch/lcg/external/gcc/4.9.1/x86_64-slc6/setup.sh 
