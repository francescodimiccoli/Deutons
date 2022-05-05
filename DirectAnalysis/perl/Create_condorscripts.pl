#!/usr/bin/perl

use warnings;
print "Printing CONDOR scripts for Analysis...\n\n";
#chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
$workdir = "/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/";
print "Printed: Work Dir. = ".$workdir."\n\n";
$njobs=$ARGV[2];
$outdir="/eos/ams/user/f/fdimicco/AnalysisFiles";
$outdir_ntuples="/eos/ams/user/f/fdimicco/AnalysisNTuples";
$outdir_plots="/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/AnalysisPlots";
#creating sum scripts and output directory
#

$queue = "ams1nd";
$offset = 0;

print "Resetting work directories...\n";

system ("$workdir/perl/SumInputFilesScripts.pl $ARGV[0] $ARGV[1] $ARGV[2] $offset");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]");
system ("mkdir $outdir_ntuples/$ARGV[0]-$ARGV[1]");
system ("mkdir $outdir_plots/$ARGV[0]-$ARGV[1]");


system ("rm -rf $workdir/perl/AnalysisScripts/$ARGV[0]-$ARGV[1]");
system ("mkdir  $workdir/perl/AnalysisScripts/$ARGV[0]-$ARGV[1]");

system ("mkdir $workdir/perl/AnalysisScripts/$ARGV[0]-$ARGV[1]/logs");
system ("mkdir $workdir/perl/AnalysisScripts/$ARGV[0]-$ARGV[1]/out");
system ("mkdir $workdir/perl/AnalysisScripts/$ARGV[0]-$ARGV[1]/err");


system ("rm -rf $workdir/perl/PlottingScripts/$ARGV[0]-$ARGV[1]");
system ("mkdir  $workdir/perl/PlottingScripts/$ARGV[0]-$ARGV[1]");

system ("mkdir $workdir/perl/PlottingScripts/$ARGV[0]-$ARGV[1]/logs");
system ("mkdir $workdir/perl/PlottingScripts/$ARGV[0]-$ARGV[1]/out");
system ("mkdir $workdir/perl/PlottingScripts/$ARGV[0]-$ARGV[1]/err");


print "Creating job files...";

$bookntuples=0;
$bookplots=1;
$booklat=0;
$bookhecont=0;
$bookflux=1;
$bookinfos=0;

if($bookflux) {
	open(OUT,">","$workdir/perl/AnalysisScripts/script$ARGV[0]-$ARGV[1].sh");

	print OUT "#!/bin/bash

		WDIR=`mktemp -d -p \$PWD`;
	cd \$WDIR;
	export WORKDIR=$workdir;
	source \$WORKDIR/mysetenv.sh;\n";

	print OUT  "\$WORKDIR/Analysis \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListDT\$1.txt \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListMC\$1.txt  \$WDIR/Analysis\$1-$ARGV[0]-$ARGV[1].root 1;  \n\n";

	print OUT "cp *.root* $outdir/$ARGV[0]-$ARGV[1]/;\n
		rm -fr \$WDIR;\n";

	close (OUT);
	system("chmod +x $workdir/perl/AnalysisScripts/script$ARGV[0]-$ARGV[1].sh");


	open(OUT,">","$workdir/perl/AnalysisScripts/Condor_script$ARGV[0]-$ARGV[1].sub");
	print OUT "executable	= $workdir/perl/AnalysisScripts/script$ARGV[0]-$ARGV[1].sh\narguments	= \$(ProcId)\noutput	= $workdir/perl/AnalysisScripts/$ARGV[0]-$ARGV[1]/out/script$ARGV[0]-$ARGV[1].\$(ProcId).out\nerror	= $workdir/perl/AnalysisScripts/$ARGV[0]-$ARGV[1]/err/script$ARGV[0]-$ARGV[1].\$(ProcId).err\nlog	= $workdir/perl/AnalysisScripts/$ARGV[0]-$ARGV[1]/logs/script$ARGV[0]-$ARGV[1].\$(ProcId).log\n+JobFlavour = \"testmatch\"\nqueue $njobs"; 


}


if($bookplots) {
	open(OUT2,">","$workdir/perl/PlottingScripts/script$ARGV[0]-$ARGV[1].sh");

	print OUT2 "#!/bin/bash

		WDIR=`mktemp -d -p \$PWD`;
	cd \$WDIR;
	export WORKDIR=$workdir;
	source \$WORKDIR/mysetenv.sh;\n";

	print OUT2  "\$WORKDIR/Distributions_Plotter \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListDT\$1.txt \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListMC\$1.txt  \$WDIR/Analysis\$1-$ARGV[0]-$ARGV[1].root 1;  \n\n";

	print OUT2 "cp *.root* $outdir_plots/$ARGV[0]-$ARGV[1];\n
		rm -fr \$WDIR;\n";

	close (OUT2);
	system("chmod +x $workdir/perl/PlottingScripts/script$ARGV[0]-$ARGV[1].sh");


	open(OUT2,">","$workdir/perl/PlottingScripts/Condor_script$ARGV[0]-$ARGV[1].sub");
	print OUT2 "executable	= $workdir/perl/PlottingScripts/script$ARGV[0]-$ARGV[1].sh\narguments	= \$(ProcId)\noutput	= $workdir/perl/PlottingScripts/$ARGV[0]-$ARGV[1]/out/script$ARGV[0]-$ARGV[1].\$(ProcId).out\nerror	= $workdir/perl/PlottingScripts/$ARGV[0]-$ARGV[1]/err/script$ARGV[0]-$ARGV[1].\$(ProcId).err\nlog	= $workdir/perl/PlottingScripts/$ARGV[0]-$ARGV[1]/logs/script$ARGV[0]-$ARGV[1].\$(ProcId).log\n+JobFlavour = \"workday\"\nqueue $njobs"; 


}




