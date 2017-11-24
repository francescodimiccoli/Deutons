#/usr/bin/perl

use warnings;
print "Printing LSF scripts for Analysis...\n\n";
chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
print "Printed: Work Dir. = ".$workdir."\n\n";
$njobs=300;

#creating sum scripts and output directory
#

system ("perl SumInputFilesScripts.pl $ARGV[0] $ARGV[1] $njobs");
system ("rm  $workdir/AnalysisFiles/$ARGV[0]/Result.root");
system ("mkdir $workdir/AnalysisFiles/$ARGV[0]");
system ("rm -r $workdir/logs/$ARGV[0]");
system ("mkdir $workdir/logs/$ARGV[0]");


for($j=0;$j<$njobs;$j++)
{
		open(OUT,">","$workdir/lsf/lsf$j.tcsh");

		print OUT "#!/bin/bash

			export WORKDIR=$workdir;
			source \$WORKDIR/../amsvar_cvmfs.sh;\n";

		print OUT  "\$WORKDIR/CountsExtraction_Parallel \$WORKDIR/InputFileLists/FileListDT$j.root \$WORKDIR/InputFileLists/FileListMC$j.root   \$WORKDIR/AnalysisFiles/$ARGV[0]/Result$j.root 1 >> \$WORKDIR/logs/$ARGV[0]/log$j.log;\n\n";

		close (OUT);

}

print "Launching Calibration jobs...";

for($j=0;$j<$njobs;$j++)
{

system("chmod +x $workdir/lsf/lsf$j.tcsh");
#system("bsub -q ams -o $workdir/lsf/lsf$j.out -e $workdir/err/lsf$j.err $workdir/lsf/lsf$j.tcsh >>$workdir/lsf/lsf$j.log\n");
}


