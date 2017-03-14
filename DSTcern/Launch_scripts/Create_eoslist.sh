#!/bin/bash 
#source /afs/cern.ch/project/eos/installation/ams/etc/setup.sh

rm eos_data.txt
ls /storage/gpfs_ams/ams/Rec/2014/ISS.B950/pass6/ | grep 'root' | awk '{print $1 }' > eos_data.txt
