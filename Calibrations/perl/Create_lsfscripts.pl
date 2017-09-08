#/usr/bin/perl

use warnings;
print "Printing LSF scripts for Calibrations...\n\n";
chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
print "Printed: Work Dir. = ".$workdir."\n\n";

$njobs=300;

#creating sum scripts and output directory
#

system ("perl $workdir/perl/Sommaisto2.pl $ARGV[0] $ARGV[1]");
system ("rm -r $workdir/Calibfiles/$ARGV[0]");
system ("mkdir $workdir/Calibfiles/$ARGV[0]");


for($j=0;$j<$njobs;$j++)
{
		open(OUT,">","../lsf/lsf$j.tcsh");

		print OUT "#!/bin/bash

			export WORKDIR=$workdir;
			source \$WORKDIR/../amsvar.sh;
			sh \$WORKDIR/../../MAIN/SumScripts/SommaistoMC$j.sh;
			sh \$WORKDIR/../../MAIN/SumScripts/Sommaisto$j.sh;
			\$WORKDIR/Resolution_Calibrations  /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommadati/sommadati$j.root  /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/temp/sommaMC$j.root  \$WORKDIR/Calibfiles/$ARGV[0]/Calib_data$j.root  >> \$WORKDIR/logs/log$j.log;\n";
		
       	 	print OUT  "hadd -k \$WORKDIR/Calibfiles/$ARGV[0]/tmp.root \$WORKDIR/Calibfiles/$ARGV[0]/Calib_data$j.root \$WORKDIR/Calibfiles/$ARGV[0]/Calib.root;\n\n";
                print OUT  "mv \$WORKDIR/Calibfiles/$ARGV[0]/tmp.root \$WORKDIR/Calibfiles/$ARGV[0]/Calib.root;\n\n";
			

		close (OUT);

}

print "Launching Calibration jobs...";

for($j=0;$j<$njobs;$j++)
{

system("chmod +x $workdir/lsf/lsf$j.tcsh");
system("bsub -q ams -o $workdir/lsf/lsf$j.out -e $workdir/err/lsf$j.err $workdir/lsf/lsf$j.tcsh >>$workdir/lsf/lsf$j.log\n");

}
