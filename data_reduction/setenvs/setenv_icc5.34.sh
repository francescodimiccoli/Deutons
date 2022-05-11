# compilers
source /afs/cern.ch/sw/lcg/external/gcc/4.9.3/x86_64-slc6/setup.sh /afs/cern.ch/sw/lcg/external
source /cvmfs/projects.cern.ch/intelsw/psxe/linux/all-setup.sh
source /cvmfs/projects.cern.ch/intelsw/psxe/linux/x86_64/2017/compilers_and_libraries_2017.2.174/linux/bin/compilervars.sh intel64
export INTELDIR=/cvmfs/projects.cern.ch/intelsw/psxe/linux/x86_64/2017
export INTELVER=compilers_and_libraries_2017.2.174
export CC=$INTELDIR/$INTELVER/linux/bin/intel64/icc
export CXX=$INTELDIR/$INTELVER/linux/bin/intel64/icpc
export LD=$INTELDIR/$INTELVER/linux/bin/intel64/icpc

# ROOT
export XRDLIB=/afs/cern.ch/ams/local2/opt/xrootd-icc64-17
export ROOTSYS=/afs/cern.ch/exp/ams/Offline/root/Linux/root-v5-34-9-icc64.17-slc6
export PATH=$ROOTSYS/bin:$PATH:/afs/cern.ch/ams/opt/intel/Compiler/11.1/073/bin/intel64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$XRDLIB/lib64:/afs/cern.ch/ams/local2/opt/xrootd-icc64-11/lib64:/afs/cern.ch/ams/opt/intel/Compiler/11.1/073/idb/lib/intel64/:/afs/cern.ch/ams/local2/opt/intel/Compiler/11.1/073/lib/intel64/:$ROOTSYS/lib:/opt/intel/compiler80/lib

# AMS env
export CVSROOT=:ext:`whoami`@lxplus.cern.ch:/afs/cern.ch/ams/Offline/CVS
export CVS_RSH=ssh
export CVSEDITOR=emacs
export Offline=/cvmfs/ams.cern.ch/Offline
export AMSDataDir=$Offline/AMSDataDir
export AMSDataDirRW=$Offline/AMSDataDirRW
export AMSDataBaseEnv=$Offline/AMSDataDir/DataBase/
export CASTORSTATIC=1
export AMSDynAlignment=$AMSDataDirRW/ExtAlig/AlignmentFiles/
export AMSICC=1
export AMSP=1
ulimit -c 0
ulimit -s unlimited
ulimit -d unlimited
export PGTRACK=1
export G4MULTITHREADED=
export G4AMS=
export GENFIT=1
export AMSWD=$MYEOS/AMS
export PATH=/eos/ams/group/dbar/installed/bin:$PATH
