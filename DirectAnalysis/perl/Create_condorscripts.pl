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


system ("mkdir $outdir_ntuples/$ARGV[0]-$ARGV[1]/Ntuples");
system ("mkdir $outdir_plots/$ARGV[0]-$ARGV[1]/Plots");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/LatRew");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/HeliumFragm");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/MCEfficiency");		
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/Counts");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/EffCorr");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/Flux");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/Infos");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/Tests");

system ("rm -r $workdir/perl/AnalysisScripts/*.log");
system ("rm -r $workdir/perl/AnalysisScripts/*.out");
system ("rm -r $workdir/err/");
system ("mkdir $workdir/err/");	
system ("rm -r $workdir/logs/");
system ("mkdir $workdir/logs/");	
system ("rm -r $outdir/../logs/$ARGV[0]-$ARGV[1]");
system ("mkdir $outdir/../logs/$ARGV[0]-$ARGV[1]");

print "Creating job files...";

$bookntuples=0;
$bookplots=0;
$booklat=0;
$bookhecont=0;
$bookflux=1;
$bookinfos=0;


		open(OUT,">","$workdir/perl/AnalysisScripts/script$ARGV[0]-$ARGV[1].sh");

		print OUT "#!/bin/bash

			export WORKDIR=$workdir;
			source \$WORKDIR/mysetenv.sh;\n";

			if($bookinfos){
				print OUT  "\$WORKDIR/ExtractInfos \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListDT\$1.txt \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListMC\$1.txt  $outdir/$ARGV[0]-$ARGV[1]/Infos$j-$offset.txt 1 >> /tmp/fdimicco/log\$1.log;\n\n";
				print OUT "cp /tmp/fdimicco/Infos\$1-$offset.txt $outdir/$ARGV[0]-$ARGV[1]/Infos;\n";
			}
		

			if($bookntuples){
				print OUT  "\$WORKDIR/Ntuple_Maker \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListDT\$1.txt \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListMC\$1.txt  $outdir/$ARGV[0]-$ARGV[1]/ 1 >> $workdir/logs/log\$1.log;\n\n";
			}
		
			if($bookplots){

				print OUT  "\$WORKDIR/Distributions_Plotter \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListDT\$1.txt \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListMC\$1.txt  $outdir/$ARGV[0]-$ARGV[1]//Plots\$1-$offset.root 1 >> $workdir/logs/log\$1.log;\n\n";
			}               
		
			if($booklat){

				print OUT  "\$WORKDIR/LatitudeReweighter \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListDT\$1.txt \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListMC\$1.txt  $outdir/$ARGV[0]-$ARGV[1]//LatRew\$1-$offset.root 1 >> $workdir/logs/log\$1.log;\n\n";
			}
		
			if($bookhecont){

				print OUT  "\$WORKDIR/HeliumContamination_Parallel \$WORKDIR/InputNtupleLists/FileListDT\$1.txt \$WORKDIR/InputNtupleLists/FileListMC\$1.txt  $outdir/$ARGV[0]-$ARGV[1]/HeliumFragm\$1-$offset.root 1 >> $workdir/logs/log\$1.log;\n\n";
			}
		

			if($bookflux){

				print OUT  "\$WORKDIR/Analysis \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListDT\$1.txt \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListMC\$1.txt  $outdir/$ARGV[0]-$ARGV[1]//Analysis\$1-$ARGV[0]-$ARGV[1].root 1 > $workdir/logs/log\$1.log 2> $workdir/err/err\$1.log \n\n";
			}
			

		close (OUT);


system("chmod +x $workdir/perl/AnalysisScripts/script$ARGV[0]-$ARGV[1].sh");


#open(OUT,">","$outdir/$ARGV[0]-$ARGV[1]/Counts/Condor_script.sub");
open(OUT,">","$workdir/perl/AnalysisScripts/Condor_script$ARGV[0]-$ARGV[1].sub");
print OUT "executable	= $workdir/perl/AnalysisScripts/script$ARGV[0]-$ARGV[1].sh\narguments	= \$(ProcId)\noutput	= $workdir/perl/AnalysisScripts/script$ARGV[0]-$ARGV[1].\$(ProcId).out\nerror	= $workdir/perl/AnalysisScripts/script$ARGV[0]-$ARGV[1].\$(ProcId).err\nlog	= $workdir/perl/AnalysisScripts/script$ARGV[0]-$ARGV[1].\$(ProcId).log\n+JobFlavour = \"tomorrow\"\nqueue $njobs"; 




