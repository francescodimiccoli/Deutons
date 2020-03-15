#!/usr/bin/perl

use warnings;
print "Printing CONDOR scripts for Analysis...\n\n";
#chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
$workdir = "/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/";
print "Printed: Work Dir. = ".$workdir."\n\n";
$njobs=$ARGV[2];
$outdir="/eos/ams/user/f/fdimicco/AnalysisFiles";
#creating sum scripts and output directory

$queue = "ams1nd";
$offset = 0;

print "Resetting work directories...\n";


system ("mkdir $workdir/perl/SumScripts/");

print "Creating job files...";



		open(OUT,">","$workdir/perl/SumScripts/script$ARGV[0]-$ARGV[1].sh");

		print OUT "#!/bin/bash

			export WORKDIR=$workdir;
			source \$WORKDIR/mysetenv.sh;\n";

				print OUT  "$workdir/perl/SumResults.pl $ARGV[0]-$ARGV[1] \$1";
			

		close (OUT);


system("chmod +x $workdir/perl/SumScripts/script$ARGV[0]-$ARGV[1].sh");


open(OUT,">","$workdir/perl/SumScripts/Condor_sum$ARGV[0]-$ARGV[1].sub");
print OUT "executable	= $workdir/perl/SumScripts/script$ARGV[0]-$ARGV[1].sh\narguments	= \$(ProcId)\noutput	= $workdir/perl/SumScripts/script$ARGV[0]-$ARGV[1].\$(ProcId).out\nerror	= $workdir/perl/SumScripts/script$ARGV[0]-$ARGV[1].\$(ProcId).err\nlog	= $workdir/perl/SumScripts/script$ARGV[0]-$ARGV[1].\$(ProcId).log\n+JobFlavour = \"microcentury\"\nqueue $njobs"; 

close (OUT);
chdir "$workdir/perl/SumScripts/";
system("condor_submit Condor_sum$ARGV[0]-$ARGV[1].sub");
chdir "$workdir/perl/";
