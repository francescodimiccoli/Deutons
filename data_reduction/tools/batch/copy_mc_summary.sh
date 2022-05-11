#!/bin/sh
BASE=/eos/ams/group/dbar
RUNLISTDIR=/afs/cern.ch/work/o/oliva/dbar/release_v5/data_reduction/analysis/batch/runlist
RUNLISTS="
release_v5_e1_vdev_190107_full_Pr.B1200_pr.pl1.05100.4_00.txt
release_v5_e1_vdev_190107_full_Pr.B1200_pr.pl1.l1.054000.4_00.txt
release_v5_e1_vdev_190107_full_el.B1200_el.pl1.0_25200.txt
release_v5_e1_vdev_190107_full_el.B1200_el.pl1.2004000.txt
release_v5_e1_vdev_190107_full_He.B1200_he4.pl1.21000.4_00.txt
release_v5_e1_vdev_190107_full_He.B1200_he4.pl1.l1.150T.4_00.txt
release_v5_e1_vdev_190107_full_He.B1200_he4.pl1.l1.24000.4_00.txt
release_v5_e1_vdev_190107_full_He.B1200_he4.pl1.l19.216000.4_00.txt
release_v5_e1_vdev_190107_full_Pr.B1128_pr.pl1phpsa.5200.txt
"

for RUNLIST in $RUNLISTS
do
  DIR=`echo $RUNLIST | sed -e 's/.txt//g' | sed -e 's:_:/:2' -e 's:_:/:4' -e 's:_:/:4' -e 's:_:/:4' `
  cat /eos/ams/user/o/oliva/inbox_check/mc_summary_$RUNLIST* | grep -v "          0" | sort | uniq > mc_summary.txt
  N1=`ls $BASE/$DIR/*root | wc -l | cut -d" " -f1 `
  N2=`wc -l mc_summary.txt | cut -d" " -f1 `
  N3=`wc -l $RUNLISTDIR/$RUNLIST`
  echo $N1 $N2 $N3 ... 
  echo mc_summary.txt to $BASE/$DIR ... 
  rm $BASE/$DIR/mc_summary.txt 
  cp mc_summary.txt $BASE/$DIR
done
