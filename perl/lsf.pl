#/usr/bin/perl
#source \$WORKDIR/amsvar.csh
#
use warnings;
for($j=0;$j<100;$j++)
{
if($ARGV[1]==0){
open(OUT,">","/storage/gpfs_ams/ams/users/fdimicco/lsf/lsf$j.tcsh");

print OUT "#!/bin/bash

export WORKDIR=/storage/gpfs_ams/ams/users/fdimicco/Deutons/
source \$WORKDIR/amsvar.sh
sh \$WORKDIR/MAIN/Sommaisto$j.sh;
\$WORKDIR/AnalyzeDATA.exe $j $ARGV[0]>> \$WORKDIR/logs/log$j.log;\n
\$WORKDIR/AnalyzeMC.exe $j $ARGV[0]>> \$WORKDIR/logs/logMC$j.log;\n
\$WORKDIR/CodesforAnalysis/Analisi $ARGV[0] 0 $j >>\$WORKDIR/logs/logAnalisi$j.log;\n";
}
if($ARGV[1]==1){
open(OUT,">","/storage/gpfs_ams/ams/users/fdimicco/lsf/lsf$j.tcsh");

print OUT "#!/bin/bash

export WORKDIR=/storage/gpfs_ams/ams/users/fdimicco/Deutons
source \$WORKDIR/amsvar.sh
\$WORKDIR/CodesforAnalysis/Analisi $ARGV[0] 0 $j >>\$WORKDIR/logs/logAnalisi$j.log;\n";
}

close (OUT);

}


