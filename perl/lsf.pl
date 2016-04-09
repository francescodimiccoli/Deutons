#/usr/bin/perl
#source \$WORKDIR/amsvar.csh
#
use warnings;
for($j=0;$j<100;$j++)
{
if($ARGV[1]==0){
open(OUT,">","/storage/gpfs_ams/ams/users/fdimicco/Deutons/lsf/lsf$j.tcsh");

print OUT "#!/bin/bash

export WORKDIR=/storage/gpfs_ams/ams/users/fdimicco/Deutons/
source \$WORKDIR/amsvar.sh
sh /storage/gpfs_ams/ams/users/fdimicco/MAIN/Sommaisto$j.sh;
\$WORKDIR/Ntuple-making/AnalyzeDATA.exe $j $ARGV[0]>> \$WORKDIR/logs/log$j.log;\n
\$WORKDIR/Ntuple-making/AnalyzeMC.exe $j $ARGV[0]>> \$WORKDIR/logs/logMC$j.log;\n
\$WORKDIR/CodesforAnalysis/Analysis $ARGV[0] 0 $j >>\$WORKDIR/logs/logAnalysis$j.log;\n";
}
if($ARGV[1]==1){
open(OUT,">","/storage/gpfs_ams/ams/users/fdimicco/Deutons/lsf/lsf$j.tcsh");

print OUT "#!/bin/bash

export WORKDIR=/storage/gpfs_ams/ams/users/fdimicco/Deutons
source \$WORKDIR/amsvar.sh
\$WORKDIR/CodesforAnalysis/Analysis $ARGV[0] 0 $j >>\$WORKDIR/logs/logAnalysis$j.log;\n";
}

close (OUT);

}


