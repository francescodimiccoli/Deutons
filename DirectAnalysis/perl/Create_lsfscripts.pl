#!/usr/bin/perl

use warnings;
print "Printing LSF scripts for Analysis...\n\n";
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


system ("rm -r $workdir/lsf/");
system ("mkdir $workdir/lsf/");
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
		open(OUT,">","$workdir/lsf/lsf$j.tcsh");

		print OUT "#!/bin/bash

			export WORKDIR=$workdir;
			source \$WORKDIR/../amsvar_cvmfs.sh;\n";

			if($bookinfos){
				print OUT "xrdcp -f $outdir/$ARGV[0]-$ARGV[1]/Infos/Infos$j.txt /tmp/fdimicco/;\n";
				print OUT  "\$WORKDIR/ExtractInfos \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/Infos$j-$offset.txt 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/fdimicco/Infos$j-$offset.txt $outdir/$ARGV[0]-$ARGV[1]/Infos;\n";
			}
		

			if($bookntuples){
				print OUT "rm /tmp/fdimicco/*.root;\n";
				print OUT  "\$WORKDIR/Ntuple_Maker \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/ 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "(cd /tmp/fdimicco/ && ls *.root*) >> /tmp/fdimicco/Ntuple$j.check;\n";
				print OUT "xrdcp -f /tmp/fdimicco/*.root* $outdir_ntuples/$ARGV[0]-$ARGV[1]/Ntuples;\n";
				print OUT "xrdcp -f /tmp/fdimicco/Ntuple$j.check  $outdir_ntuples/$ARGV[0]-$ARGV[1]/Ntuples;\n";
				print OUT "rm -r /tmp/fdimicco/*.root*;\n";	
			}
		
			if($bookplots){

				print OUT  "\$WORKDIR/Distributions_Plotter \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/Plots$j-$offset.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/fdimicco/Plots$j-$offset.root $outdir_plots/$ARGV[0]-$ARGV[1]/Plots;\n";
			}               
		
			if($booklat){

				print OUT  "\$WORKDIR/LatitudeReweighter \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/LatRew$j-$offset.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/fdimicco/LatRew$j-$offset.root $outdir/$ARGV[0]-$ARGV[1]/LatRew;\n";
			}
		
			if($bookhecont){

				print OUT  "\$WORKDIR/HeliumContamination_Parallel \$WORKDIR/InputNtupleLists/FileListDT$j.txt \$WORKDIR/InputNtupleLists/FileListMC$j.txt  /tmp/fdimicco/HeliumFragm$j-$offset.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/fdimicco/HeliumFragm$j-$offset.root $outdir/$ARGV[0]-$ARGV[1]/HeliumFragm;\n";
			}
		

			if($bookflux){

				print OUT  "\$WORKDIR/Analysis \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/Analysis$j-$offset.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/fdimicco/Analysis$j-$offset.root_Counts $outdir/$ARGV[0]-$ARGV[1]/Counts;\n";
				print OUT "xrdcp -f /tmp/fdimicco/Analysis$j-$offset.root_Eff    $outdir/$ARGV[0]-$ARGV[1]/MCEfficiency;\n";
				print OUT "xrdcp -f /tmp/fdimicco/Analysis$j-$offset.root_Corr   $outdir/$ARGV[0]-$ARGV[1]/EffCorr;\n";
				print OUT "xrdcp -f /tmp/fdimicco/Analysis$j-$offset.root_Flux   $outdir/$ARGV[0]-$ARGV[1]/Fluxes;\n";
				print OUT "xrdcp -f /tmp/fdimicco/Analysis$j-$offset.root_Test   $outdir/$ARGV[0]-$ARGV[1]/Tests;\n";
			}
			
			print OUT "xrdcp -f /tmp/fdimicco/log$j.log $outdir/../logs/$ARGV[0]-$ARGV[1]/;\n";


		close (OUT);


}

for($j=0;$j<$njobs;$j++)
{
system("chmod +x $workdir/lsf/lsf$j.tcsh");
}

print "Launching job files...";
$joblaunched = 0;
$jobrunning = 0;	
$jobs = 0;
	
for($loop=0; $loop<5; $loop++)	{
	$joblaunched = 0;
	system("kinit -R");
	$jobs = `bjobs -q $queue| wc -l`;
	$control=0;
	
	while($joblaunched<$njobs){
		$jobrunning = `bjobs -q $queue |grep RUN|wc -l`;
		$jobs = `bjobs -q $queue| wc -l`;
		$check;
		$control=0;
		while($jobs<600 and $joblaunched<$njobs) {

			$control=0;
			if($bookinfos) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/Infos/Infos$joblaunched.txt| awk '{print\$5}'`; if($check<1000) {print $check; $control=1}} 
			if($bookntuples) { 
				$checkfile = `ls -la $outdir_ntuples/$ARGV[0]-$ARGV[1]/Ntuples/Ntuple$joblaunched.check`;
				if($checkfile ne ""){
					print "Check file found!!\n";
					$filename = `cat $outdir_ntuples/$ARGV[0]-$ARGV[1]/Ntuples/Ntuple$joblaunched.check`;
					system("root -l -q $outdir_ntuples/$ARGV[0]-$ARGV[1]/Ntuples/$filename 2> $outdir_ntuples/$ARGV[0]-$ARGV[1]/Ntuples/Ntuple$joblaunched.test"); 
					$test = `cat $outdir_ntuples/$ARGV[0]-$ARGV[1]/Ntuples/Ntuple$joblaunched.test`;
					if($test ne "") {print $test; $control=1}
				}
				if($checkfile eq "") {$control=1}	
			} 
			if($bookplots) { $check = `ls -la $outdir_plots/$ARGV[0]-$ARGV[1]/Plots/Plots$joblaunched-$offset.root| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($booklat) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/LatRew/LatRew$joblaunched-$offset.root| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($bookhecont) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/HeliumFragm/HeliumFragm$joblaunched-$offset.root| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($bookflux) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/MCEfficiency/Analysis$joblaunched-$offset.root_Eff| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($bookflux) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/Counts/Analysis$joblaunched-$offset.root_Counts| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($bookflux) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/EffCorr/Analysis$joblaunched-$offset.root_Corr| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($bookflux) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/Fluxes/Analysis$joblaunched-$offset.root_Flux| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($control){
				print "job result nr $joblaunched not found\n"; 
				system("bsub -q $queue -o $workdir/lsf/lsf$joblaunched.out -e $workdir/err/lsf$joblaunched.err $workdir/lsf/lsf$joblaunched.tcsh\n");
			}
			else {
				print "job $joblaunched ok\n";
			}	
			$joblaunched++;
			$jobs = `bjobs -q $queue|wc -l`;
		}
		print "Total Jobs launched: $joblaunched\n";
		print "Jobs currently running: $jobrunning\n";
		sleep(10);
	}
}



