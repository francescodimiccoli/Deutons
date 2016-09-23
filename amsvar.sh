#export GCCDIR=/opt/rh/devtoolset-3/root/usr/bin/

export AMSBASE=/storage/gpfs_ams/ams/users/basara/AMSSoft/
export AMSSYS=linux_slc6_gcc64

export CERN=$AMSBASE/$AMSSYS/
export CERN_LEVEL=2005/
export CERN_ROOT=$CERN/$CERN_LEVEL
export CERNDIR=$CERN/$CERN_LEVEL
export AMSLIB=$CERNDIR
export AMSDataDir=/opt/exp_software/ams/AMSDataDir
export PATH=$ROOTSYS/bin:$CERN_ROOT/bin:$GCCDIR:$PATH
export CLHEP_BASE_DIR=$AMSBASE/$AMSSYS/CLHEP
export CLHEP_INCLUDE_DIR=$AMSBASE/$AMSSYS/CLHEP/include
export CLHEP_LIB_DIR=$AMSBASE/$AMSSYS/CLHEP/lib
export G4INSTALL=$AMSBASE/$AMSSYS/geant4_ams

export QTDIR=/opt/exp_software/ams/qt-4.8.3/
export QTLIB=$QTDIR/lib/
export QTINC=$QTDIR/include/
export QTDIR_STATIC=$QTDIR
export LD_LIBRARY_PATH=$QTLIB:$ROOTSYS/lib:$LD_LIBRARY_PATH
export G4SYSTEM=Linux-g++
#export AMS_ACQT_INTERFACE=1

export AMSP=1
export PGTRACK=1
export ECALBDT=1
export AMSSRC=/storage/gpfs_ams/ams/users/basara/AMSSoft/
export AMSWD=$AMSSRC
export AMSGeoDir=$AMSSRC/display/ams02/
export CVSROOT=:ext:lbasara@ams.cern.ch:/afs/cern.ch/ams/Offline/CVS
export CVS_RSH=ssh
export NOCASTOR=1
export VERBOSE=1
