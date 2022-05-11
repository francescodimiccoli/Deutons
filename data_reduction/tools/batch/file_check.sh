#!/bin/sh

# encapsulation of the file_check executable 
EXEC=/afs/cern.ch/work/o/oliva/dbar/release_v5/data_reduction/tools/file_check

LIST=$1
START=$2
END=$3

LISTNAME=`basename $LIST `
OUTNAME=`printf "file_check_%s_%06d_%06d.txt" $LISTNAME $START $END `

for FILE in `head -$((END+1)) $LIST | tail -$((END-START+1)) `
do
  $EXEC $FILE 
  RETURN=$?
  echo $FILE $RETURN >> $OUTNAME 
done
