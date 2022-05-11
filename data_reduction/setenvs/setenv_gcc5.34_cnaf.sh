# compilers
. /cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/setup.sh
#. /cvmfs/sft.cern.ch/lcg/releases/gcc/7.3.0/x86_64-slc6/setup.sh
gcc --version
export PATH=/cvmfs/sft.cern.ch/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/bin/:$PATH
python --version

# ROOT
echo "Exporting ROOT vars"
#export ROOTSYS=/cvmfs/ams.cern.ch/Offline/root/Linux/root-v5-34-9-gcc64-slc6
export ROOTSYS=/opt/exp_software/ams/AMSsoft_modified/linux_slc6_gcc64/root_v5.34ams/
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib

# AMS env
echo "Exporting AMS vars"
export CVSROOT=:ext:`whoami`@lxplus.cern.ch:/afs/cern.ch/ams/Offline/CVS
export CVS_RSH=ssh
export CVSEDITOR=emacs
export Offline=/cvmfs/ams.cern.ch/Offline
export AMSDataDir=$Offline/AMSDataDir
export AMSDataDirRW=$Offline/AMSDataDirRW
export AMSDataBaseEnv=$Offline/AMSDataDir/DataBase/
#export CASTORSTATIC=1
export AMSDynAlignment=$AMSDataDirRW/ExtAlig/AlignmentFiles/
export AMSP=1
export PGTRACK=1
export G4MULTITHREADED=
export G4AMS=
export GENFIT=1
export LD_LIBRARY_PATH=$Offline/AMSsoft/linux_slc6_gcc64/xrootd/lib64/:$LD_LIBRARY_PATH
export PATH=$Offline/AMSsoft/linux_slc6_gcc64/xrootd/bin/:$PATH
export PATH=/cvmfs/ams.cern.ch/Offline/dbar/public/installed/bin:$PATH
export EOS_MGM_URL=root://eosams.cern.ch

echo "Tweaking memory management"
ulimit -c 0
ulimit -s unlimited
ulimit -d unlimited
