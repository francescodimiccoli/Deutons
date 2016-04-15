#/usr/bin/perl

use warnings;
print "Printing LSF scripts..\n\n";
chomp($workdir =`pwd -P |sed 's\\perl\\\\g'|sed 's\\Lancia\\\\g '`);
print "Printed: Work Dir. = ".$workdir."\n\n";
 
for($j=0;$j<100;$j++)
{
if($ARGV[1]==0){
open(OUT,">","../lsf/lsf$j.tcsh");

print OUT "#!/bin/bash

export WORKDIR=$workdir
source \$WORKDIR/amsvar.sh
sh \$WORKDIR/../MAIN/Sommaisto$j.sh;
\$WORKDIR/Ntuple-making/AnalyzeDATA.exe $j $ARGV[0]>> \$WORKDIR/logs/log$j.log;\n
\$WORKDIR/Ntuple-making/AnalyzeMC.exe $j $ARGV[0]>> \$WORKDIR/logs/logMC$j.log;\n
\$WORKDIR/CodesforAnalysis/Analysis $ARGV[0] 0 $j \$WORKDIR >>\$WORKDIR/logs/logAnalysis$j.log;\n";
}
if($ARGV[1]==1){
open(OUT,">","$workdir/lsf/lsf$j.tcsh");

print OUT "#!/bin/bash

export WORKDIR=$workdir
source \$WORKDIR/amsvar.sh
\$WORKDIR/CodesforAnalysis/Analysis $ARGV[0] 0 $j \$WORKDIR >>\$WORKDIR/logs/logAnalysis$j.log;\n";
}

close (OUT);

}


