#/usr/bin/perl

use warnings;
print "Printing LSF scripts..\n\n";
chomp($workdir =`pwd -P |sed 's\\perl\\\\g'`);
print "Printed: Work Dir. = ".$workdir."\n\n";

#creating sum scripts and output directory
#

system ("perl $workdir/perl/Sommaisto2.pl $ARGV[0] $ARGV[1]");
system ("rm -r $workdir/Ntuples/$ARGV[0]");
system ("mkdir $workdir/Ntuples/$ARGV[0]");



for($j=0;$j<100;$j++)
{
		open(OUT,">","../lsf/lsf$j.tcsh");

		print OUT "#!/bin/bash

			export WORKDIR=$workdir
			source \$WORKDIR/../amsvar.sh
			sh \$WORKDIR/../../MAIN/SumScripts/Sommaisto$j.sh;
			\$WORKDIR/Analyzer.exe $ARGV[0] /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommadati/sommadati$j.root $workdir/Ntuples/$ARGV[0]/NtupleData$j.root   >> \$WORKDIR/logs/log$j.log;\n";
}


for($j=0;$j<100;$j++)
{
                open(OUT,">","../lsf/lsfMC$j.tcsh");

                print OUT "#!/bin/bash
			
			export WORKDIR=$workdir
                        source \$WORKDIR/../amsvar.sh
			sh \$WORKDIR/../../MAIN/SumScripts/SommaistoMC$j.sh;	
			\$WORKDIR/Analyzer.exe $ARGV[0] /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/temp/sommaMC$j.root $workdir/Ntuples/$ARGV[0]/NtupleMC$j.root >> \$WORKDIR/logs/logMC$j.log;\n";

}

print "Launching NTuple-making jobs...";

for($j=0;$j<100;$j++)
{

system("chmod +x $workdir/lsf/lsf$j.tcsh");
system("bsub -q ams -o $workdir/lsf/lsf$j.out -e $workdir/err/lsf$j.err $workdir/lsf/lsf$j.tcsh >>$workdir/lsf/lsf$j.log\n");

} 



for($j=0;$j<100;$j++)
{

system("chmod +x $workdir/lsf/lsfMC$j.tcsh");
system("bsub -q ams -o $workdir/lsf/lsfMC$j.out -e $workdir/err/lsfMC$j.err $workdir/lsf/lsfMC$j.tcsh >>$workdir/lsf/lsfMC$j.log\n");

} 
