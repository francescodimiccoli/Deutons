#/usr/bin/perl

use warnings;
print "Printing LSF scripts for Analysis...\n\n";
chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
print "Printed: Work Dir. = ".$workdir."\n\n";


#creating sum scripts and output directory
#

system ("rm -r $workdir/AnalysisFiles/$ARGV[0]");
system ("mkdir $workdir/AnalysisFiles/$ARGV[0]");
system ("rm -r $workdir/logs/$ARGV[0]");
system ("mkdir $workdir/logs/$ARGV[0]");



for($j=0;$j<100;$j++)
{
		open(OUT,">","../lsf/lsf$j.tcsh");

		print OUT "#!/bin/bash

			export WORKDIR=$workdir;
			source \$WORKDIR/../amsvar.sh;
			
			\$WORKDIR/CountsExtraction  /storage/gpfs_ams/ams/users/fdimicco/Deutons/Ntuple-making/Ntuples/$ARGV[0]/NtupleData$j.root /storage/gpfs_ams/ams/users/fdimicco/Deutons/Ntuple-making/Ntuples/$ARGV[0]/NtupleMC$j.root  \$WORKDIR/AnalysisFiles/$ARGV[0]/Result$j.root >> \$WORKDIR/logs/$ARGV[0]/log$j.log;\n

			\$WORKDIR/HeliumContamination  /storage/gpfs_ams/ams/users/fdimicco/Deutons/Ntuple-making/Ntuples/$ARGV[0]/NtupleData$j.root /storage/gpfs_ams/ams/users/fdimicco/Deutons/Ntuple-making/Ntuples/$ARGV[0]/NtupleMC$j.root  \$WORKDIR/AnalysisFiles/$ARGV[0]/Result$j.root  >> \$WORKDIR/logs/$ARGV[0]/log$j.log;\n";

				

			

		close (OUT);

}

print "Launching Calibration jobs...";

for($j=0;$j<100;$j++)
{

system("chmod +x $workdir/lsf/lsf$j.tcsh");
system("bsub -q ams -o $workdir/lsf/lsf$j.out -e $workdir/err/lsf$j.err $workdir/lsf/lsf$j.tcsh >>$workdir/lsf/lsf$j.log\n");

}