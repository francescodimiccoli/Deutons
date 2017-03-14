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
			sh \$WORKDIR/../MAIN/SumScripts/Sommaisto$j.sh;
			\$WORKDIR/Ntuple-making/Analyzer.exe $ARGV[0] /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommadati/sommadati$j.root /storage/gpfs_ams/ams/users/fdimicco/Deutons/Risultati/$ARGV[0]/NtupleData$j.root   >> \$WORKDIR/logs/log$j.log;\n
		
			sh \$WORKDIR/../MAIN/SumScripts/SommaistoMC$j.sh;	
			\$WORKDIR/Ntuple-making/Analyzer.exe $ARGV[0] /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/temp/sommaMC$j.root /storage/gpfs_ams/ams/users/fdimicco/Deutons/Risultati/$ARGV[0]/NtupleMC$j.root >> \$WORKDIR/logs/logMC$j.log;\n";

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


