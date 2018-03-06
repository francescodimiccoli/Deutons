#/usr/bin/perl

use warnings;
print "Printing LSF scripts for Analysis...\n\n";
chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
print "Printed: Work Dir. = ".$workdir."\n\n";
$njobs=$ARGV[2];
#$outdir="/eos/ams/group/dbar/TrentoNTuples/";
$outdir="/eos/ams/user/f/fdimicco/AnalysisFiles";
#creating sum scripts and output directory
#

system ("perl SumInputFilesScripts.pl $ARGV[0] $ARGV[1] $ARGV[2] $ARGV[3]");
system ("rm  $outdir/$ARGV[0]-$ARGV[1]/Result.root");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]");
system ("rm -r $workdir/logs/$ARGV[0]-$ARGV[1]");
system ("mkdir $workdir/logs/$ARGV[0]-$ARGV[1]");


for($j=0;$j<$njobs;$j++)
{
		open(OUT,">","$workdir/lsf/lsf$j.tcsh");

		print OUT "#!/bin/bash

			export WORKDIR=$workdir;
			source \$WORKDIR/../amsvar_cvmfs.sh;\n";

#		print OUT  "\$WORKDIR/Ntuple_Maker \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt   $outdir/$ARGV[0]-$ARGV[1]/Result$j.root 1 >> \$WORKDIR/logs/$ARGV[0]-$ARGV[1]/log$j.log;\n\n";
#		print OUT  "\$WORKDIR/HeliumContamination_Parallel \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt   $outdir/$ARGV[0]-$ARGV[1]/Result$j.root 1 >> \$WORKDIR/logs/$ARGV[0]-$ARGV[1]/log$j.log;\n\n";
		print OUT  "\$WORKDIR/CountsExtraction_Parallel \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt   $outdir/$ARGV[0]-$ARGV[1]/Result$j.root 1 >> \$WORKDIR/logs/$ARGV[0]-$ARGV[1]/log$j.log;\n\n";
#		print OUT  "\$WORKDIR/MCEfficiency_Parallel \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt   $outdir/$ARGV[0]-$ARGV[1]/Result$j.root 1 >> \$WORKDIR/logs/$ARGV[0]-$ARGV[1]/log$j.log;\n\n";
#		print OUT  "\$WORKDIR/Distributions_Plotter \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt   $outdir/$ARGV[0]-$ARGV[1]/Result$j.root 1 >> \$WORKDIR/logs/$ARGV[0]-$ARGV[1]/log$j.log;\n\n";

		close (OUT);


}

print "Launching Calibration jobs...";

for($j=0;$j<$njobs;$j++)
{
system("chmod +x $workdir/lsf/lsf$j.tcsh");
system("bsub -q 1nh -o $workdir/lsf/lsf$j.out -e $workdir/err/lsf$j.err $workdir/lsf/lsf$j.tcsh >>$workdir/lsf/lsf$j.log\n");
}


