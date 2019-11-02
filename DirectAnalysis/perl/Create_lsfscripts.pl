#!/usr/bin/perl

use warnings;
print "Printing CONDOR scripts for Analysis...\n\n";
#chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
chomp($workdir = "/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/");
print "Printed: Work Dir. = ".$workdir."\n\n";
$njobs=$ARGV[2];
$outdir="/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/AnalysisFiles";
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
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/Fluxes");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/Infos");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/Tests");


system ("rm -r $workdir/perl/AnalysisScripts/");
system ("mkdir $workdir/perl/AnalysisScripts/");
system ("rm -r $workdir/err/");
system ("mkdir $workdir/err/");	
system ("rm -r $outdir/../logs/$ARGV[0]-$ARGV[1]");
system ("mkdir $outdir/../logs/$ARGV[0]-$ARGV[1]");

print "Creating job files...";

$bookntuples=0;
$bookplots=0;
$booklat=0;
$bookhecont=0;
$bookflux=1;
$bookinfos=0;


for($j=0;$j<$njobs;$j++)
{
		open(OUT,">","$workdir/perl/AnalysisScripts/script$j.sh");

		print OUT "#!/bin/bash

			export WORKDIR=$workdir;
			source \$WORKDIR/mysetenv.sh;\n";

			if($bookinfos){
				print OUT "cp $outdir/$ARGV[0]-$ARGV[1]/Infos/Infos$j.txt /tmp/fdimicco/;\n";
				print OUT  "\$WORKDIR/ExtractInfos \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/Infos$j-$offset.txt 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "cp /tmp/fdimicco/Infos$j-$offset.txt $outdir/$ARGV[0]-$ARGV[1]/Infos;\n";
			}
		

			if($bookntuples){
				print OUT "rm /tmp/fdimicco/*.root;\n";
				print OUT  "\$WORKDIR/Ntuple_Maker \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/ 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "(cd /tmp/fdimicco/ && ls *.root*) >> /tmp/fdimicco/Ntuple$j.check;\n";
				print OUT "cp /tmp/fdimicco/*.root* $outdir_ntuples/$ARGV[0]-$ARGV[1]/Ntuples;\n";
				print OUT "cp /tmp/fdimicco/Ntuple$j.check  $outdir_ntuples/$ARGV[0]-$ARGV[1]/Ntuples;\n";
				print OUT "rm -r /tmp/fdimicco/*.root*;\n";	
			}
		
			if($bookplots){

				print OUT  "\$WORKDIR/Distributions_Plotter \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/Plots$j-$offset.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "cp /tmp/fdimicco/Plots$j-$offset.root $outdir_plots/$ARGV[0]-$ARGV[1]/Plots;\n";
			}               
		
			if($booklat){

				print OUT  "\$WORKDIR/LatitudeReweighter \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/LatRew$j-$offset.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "cp /tmp/fdimicco/LatRew$j-$offset.root $outdir/$ARGV[0]-$ARGV[1]/LatRew;\n";
			}
		
			if($bookhecont){

				print OUT  "\$WORKDIR/HeliumContamination_Parallel \$WORKDIR/InputNtupleLists/FileListDT$j.txt \$WORKDIR/InputNtupleLists/FileListMC$j.txt  /tmp/fdimicco/HeliumFragm$j-$offset.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "cp /tmp/fdimicco/HeliumFragm$j-$offset.root $outdir/$ARGV[0]-$ARGV[1]/HeliumFragm;\n";
			}
		

			if($bookflux){

				print OUT  "\$WORKDIR/Analysis \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/Analysis$j-$offset.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "cp /tmp/fdimicco/Analysis$j-$offset.root_Counts $outdir/$ARGV[0]-$ARGV[1]/Counts;\n";
				print OUT "cp /tmp/fdimicco/Analysis$j-$offset.root_Eff    $outdir/$ARGV[0]-$ARGV[1]/MCEfficiency;\n";
				print OUT "cp /tmp/fdimicco/Analysis$j-$offset.root_Corr   $outdir/$ARGV[0]-$ARGV[1]/EffCorr;\n";
				print OUT "cp /tmp/fdimicco/Analysis$j-$offset.root_Flux   $outdir/$ARGV[0]-$ARGV[1]/Fluxes;\n";
				print OUT "cp /tmp/fdimicco/Analysis$j-$offset.root_Test   $outdir/$ARGV[0]-$ARGV[1]/Tests;\n";
			}
			
			print OUT "cp /tmp/fdimicco/log$j.log $outdir/../logs/$ARGV[0]-$ARGV[1]/;\n";


		close (OUT);


}

for($j=0;$j<$njobs;$j++)
{
system("chmod +x $workdir/perl/AnalysisScripts/script$j.sh");
}

open(OUT2,">","$workdir/Jobs.sh");

for ($n=0;$n<$njobs; $n++)
{
	open(OUT,">","$workdir/perl/AnalysisScripts/Condor_script$n.sub");
	print OUT "executable	= $workdir/perl/AnalysisScripts/script$n.sh\narguments	=\noutput	= $workdir/perl/AnalysisScripts/script$n.out\nerror	= $workdir/perl/AnalysisScripts/script$n.err\nlog	= $workdir/perl/AnalysisScripts/script$n.log\n+JobFlavour = \"workday\"\nqueue" 
}

for ($n=0;$n<$njobs; $n++)
{
	print OUT2 "condor_submit $workdir/perl/AnalysisScripts/Condor_script$n.sub\n";
}





