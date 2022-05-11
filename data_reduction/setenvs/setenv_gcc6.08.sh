#!/bin/bash

export EOS_MGM_URL=root://eosams.cern.ch

. /cvmfs/sft.cern.ch/lcg/releases/gcc/4.9.3/x86_64-slc6-gcc49-opt/setup.sh
gcc --version
export PATH=/cvmfs/sft.cern.ch/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt/bin/:$PATH
python --version

#export Offline=/afs/cern.ch/exp/ams/Offline
export Offline=/cvmfs/ams.cern.ch/Offline
export NAGLIBRARY=$Offline/CERN/NagLib
export NAGDIR=$NAGLIBRARY
export AMSVERYBASE=$Offline
export AMSBASE=$AMSVERYBASE/AMSsoft
#export AMSSYS=linux_slc5_gcc64
export AMSSYS=linux_slc6_gcc64
#export ROOTSYS=$AMSBASE/$AMSSYS/root_v6.04.08ams/
export ROOTSYS=/eos/ams/group/dbar/root/root-6.08.06_gcc
export CERN=$AMSBASE/$AMSSYS/
export CERN_LEVEL=2005/
export CERN_ROOT=$CERN/$CERN_LEVEL
export CERNDIR=$CERN/$CERN_LEVEL
export AMSLIB=$CERNDIR
export AMSDataDir=$AMSVERYBASE/AMSDataDir
#export PATH=$ROOTSYS/bin:$CERN_ROOT/bin:/opt/intel/Compiler/11.1/073/bin/intel64:$PATH
#export LD_LIBRARY_PATH=$ROOTSYS/lib:/opt/intel/Compiler/11.1/073/lib/intel64:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$CERN_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export CLHEP_BASE_DIR=$AMSBASE/$AMSSYS/CLHEP
export CLHEP_INCLUDE_DIR=$AMSBASE/$AMSSYS/CLHEP/include
export CLHEP_LIB_DIR=$AMSBASE/$AMSSYS/CLHEP/lib
export G4INSTALL=$AMSBASE/$AMSSYS/geant4_ams

export G4SYSTEM=Linux-g++

# Environment variables needed to find geant4 data files:
#
# Data for neutron scattering processes,
#    distributed in a separate tar file, then placed under data
export NeutronHPCrossSections=$AMSBASE/share/geant4/data/G4NDL
#
#  Nuclear Photon evaporation data,
#    distributed with the source files under data
export G4LEVELGAMMADATA=$AMSBASE/share/geant4/data/PhotonEvaporation
#
# Data for radiative decay hadronic processes under data,
#    distributed in a separate tar file
export G4RADIOACTIVEDATA=$AMSBASE/share/geant4/data/RadiativeDecay
#
# Data for low energy electromagnetic processes,
#    distributed in a separate tar file, then placed under data
export G4LEDATA=$AMSBASE/share/geant4/data/G4EMLOW
#
# Data for elastic scattering processes,
#    distributed in a separate tar file, then placed under data
export G4ELASTICDATA=$AMSBASE/share/geant4/data/G4ELASTIC

export AMSP=1
export PGTRACK=1
export AMSSRC=$MYEOS/AMS
export AMSGeoDir=$AMSSRC/display/ams02/
export CVS_RSH=ssh
export CVSROOT=:ext:vformato@ams.cern.ch:/afs/cern.ch/exp/ams/Offline/CVS
export AMSWD=$AMSSRC
export PATH=$AMSSRC/bin:$PATH
export PATH=$AMSSRC/exe/linuxx8664gcc`root-config --version | cut -b1-4`/:$PATH

export NOCASTOR=1 #otehrwise AMSDisplay can be not compiling
export NOXROOTD=1 #added my Matt
export NORFIOD=1 #added by Matt

#export AMSROOT=/storage/gpfs_ams/ams/users/vformato/myAMSRoot/AMSRoot
#export LD_LIBRARY_PATH=$AMSROOT/lib/:$LD_LIBRARY_PATH

export ECALBDT=1
export G4AMS=
export GENFIT=1

#TrdQT
#export AMS_ACQT_INTERFACE=1
export ACROOTSOFTWARE=$AMSDataDir/v5.00/TRD
export ACROOTLOOKUPS=$ACROOTSOFTWARE/acroot/data

export DATA=/data/

export LD_LIBRARY_PATH=$AMSBASE/$AMSSYS/xrootd/lib64/:$LD_LIBRARY_PATH
export PATH=$AMSBASE/$AMSSYS/xrootd/bin/:$PATH

#Custom variables
export AMS_INC=$AMSSRC/include
export AMS_LIB=$AMSSRC/lib/linuxx8664gcc`root-config --version | cut -b1-4`
export LD_LIBRARY_PATH=$AMS_LIB:$LD_LIBRARY_PATH

echo "AMS env ready"
