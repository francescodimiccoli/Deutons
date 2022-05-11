#/bin/sh 

launch_file_check() {
  SETENV=/afs/cern.ch/work/o/oliva/dbar/data_reduction_production/setenvs/setenv_gcc6.14_cc7.sh
  OUTDIR=/eos/ams/user/o/oliva/inbox_check
  EXEC=$PWD/file_check.sh 
  LIST=$1
  NAME=`basename $LIST`
  PREFIX=`printf "file_check_%s" $NAME`
  source launch_file_list_htc.sh $SETENV $OUTDIR $EXEC $LIST 200 longlunch $PREFIX
}

#launch_file_check /afs/cern.ch/work/o/oliva/dbar/release_v5/data_reduction/analysis/batch/runlist/release_v5_e1_vdev_181025_neg_ISS.B1130_pass7.txt
launch_file_check /afs/cern.ch/work/o/oliva/dbar/release_v5/data_reduction/analysis/batch/runlist/release_v5_e1_vdev_190107_full_Pr.B1200_pr.pl1.05100.4_00.txt
launch_file_check /afs/cern.ch/work/o/oliva/dbar/release_v5/data_reduction/analysis/batch/runlist/release_v5_e1_vdev_190107_full_Pr.B1200_pr.pl1.l1.054000.4_00.txt
launch_file_check /afs/cern.ch/work/o/oliva/dbar/release_v5/data_reduction/analysis/batch/runlist/release_v5_e1_vdev_190107_full_el.B1200_el.pl1.0_25200.txt
launch_file_check /afs/cern.ch/work/o/oliva/dbar/release_v5/data_reduction/analysis/batch/runlist/release_v5_e1_vdev_190107_full_el.B1200_el.pl1.2004000.txt
launch_file_check /afs/cern.ch/work/o/oliva/dbar/release_v5/data_reduction/analysis/batch/runlist/release_v5_e1_vdev_190107_full_He.B1200_he4.pl1.21000.4_00.txt
launch_file_check /afs/cern.ch/work/o/oliva/dbar/release_v5/data_reduction/analysis/batch/runlist/release_v5_e1_vdev_190107_full_He.B1200_he4.pl1.l1.150T.4_00.txt
launch_file_check /afs/cern.ch/work/o/oliva/dbar/release_v5/data_reduction/analysis/batch/runlist/release_v5_e1_vdev_190107_full_He.B1200_he4.pl1.l1.24000.4_00.txt
launch_file_check /afs/cern.ch/work/o/oliva/dbar/release_v5/data_reduction/analysis/batch/runlist/release_v5_e1_vdev_190107_full_He.B1200_he4.pl1.l19.216000.4_00.txt
launch_file_check /afs/cern.ch/work/o/oliva/dbar/release_v5/data_reduction/analysis/batch/runlist/release_v5_e1_vdev_190107_full_Pr.B1128_pr.pl1phpsa.5200.txt

