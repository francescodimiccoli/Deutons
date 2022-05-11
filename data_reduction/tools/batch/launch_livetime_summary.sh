#/bin/sh 

# inputs 
SETENV=/afs/cern.ch/work/o/oliva/dbar/data_reduction_production/setenvs/setenv_gcc6.14_cc7.sh
OUTDIR=/eos/ams/user/o/oliva/inbox_check
EXEC=/afs/cern.ch/work/o/oliva/dbar/release_v5/data_reduction/tools/livetime_summary 
LIST=/afs/cern.ch/work/o/oliva/dbar/release_v5/data_reduction/analysis/batch/runlist/release_v5_e1_vdev_181025_neg_ISS.B1130_pass7.txt
NAME=`basename $LIST`
PREFIX=`printf "livetime_summary_%s" $NAME`

source launch_file_list_htc.sh $SETENV $OUTDIR $EXEC $LIST 100 longlunch $PREFIX 

