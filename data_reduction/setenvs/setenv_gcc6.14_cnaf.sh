opt="$1"

# compilers
export EOS_MGM_URL=root://eosams.cern.ch

echo "Selecting gcc"
. /cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/setup.sh
gcc --version
export PATH=/cvmfs/sft.cern.ch/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/bin/:$PATH
python --version
export PATH=/cvmfs/sft.cern.ch/lcg/releases/CMake/3.11.1-daf3a/x86_64-centos7-gcc48-opt/bin:$PATH
cmake --version

# ROOT
echo "Exporting ROOT vars"
export ROOTSYS=/opt/exp_software/ams/AMSsoft_modified/linux_centos7_gcc64/root_v6.14.04ams/
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
if [ "$opt" == "exportRIP" ]; then
  export ROOT_INCLUDE_PATH=/cvmfs/ams.cern.ch/Offline/dbar/public/release_v5/AMS_vdev_190318/include/
fi
echo "ROOT version" `root-config --version`

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
export NORFIOD=1
export AMSDynAlignment=$AMSDataDirRW/ExtAlig/AlignmentFiles/
export AMSP=1
export PGTRACK=1
export G4MULTITHREADED=
export G4AMS=
export GENFIT=1
export PATH=/cvmfs/ams.cern.ch/Offline/dbar/public/installed/bin:$PATH

if [ -z "$AMSWD" ]; then
  export AMSWD=/cvmfs/ams.cern.ch/Offline/dbar/public/release_v5/AMS_vdev_190318
fi

echo "Tweaking memory management"
ulimit -c 0
ulimit -s unlimited
ulimit -d unlimited
