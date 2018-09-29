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

$queue = "8nh";

print "Resetting work directories...\n";

system ("$workdir/perl/SumInputFilesScripts.pl $ARGV[0] $ARGV[1] $ARGV[2] 0");
system ("$workdir/perl/SumInputFilesScripts.pl $ARGV[0] $ARGV[1] $ARGV[2] 1");
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
$bookeff=0;
$bookcounts=0;
$bookeffcorr=1;
$bookflux=0;
$bookinfos=0;


for($j=0;$j<$njobs;$j++)
{
		open(OUT,">","$workdir/lsf/lsf$j.tcsh");

		print OUT "#!/bin/bash

			export WORKDIR=$workdir;
			source \$WORKDIR/../amsvar_cvmfs.sh;\n";

			if($bookinfos){
				print OUT "xrdcp -f $outdir/$ARGV[0]-$ARGV[1]/Infos/Infos$j.txt /tmp/fdimicco/;\n";
				print OUT  "\$WORKDIR/ExtractInfos \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/Infos$j.txt 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/fdimicco/Infos$j.txt $outdir/$ARGV[0]-$ARGV[1]/Infos;\n";
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

				print OUT "xrdcp -f $outdir_plots/$ARGV[0]-$ARGV[1]/Plots/Plots$j.root /tmp/fdimicco/;\n";
				print OUT  "\$WORKDIR/Distributions_Plotter \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/Plots$j.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/fdimicco/Plots$j.root $outdir_plots/$ARGV[0]-$ARGV[1]/Plots;\n";
			}               
		
			if($booklat){

				print OUT "xrdcp -f $outdir/$ARGV[0]-$ARGV[1]/LatRew/LatRew$j.root /tmp/fdimicco/;\n";
				print OUT  "\$WORKDIR/LatitudeReweighter \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/LatRew$j.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/fdimicco/LatRew$j.root $outdir/$ARGV[0]-$ARGV[1]/LatRew;\n";
			}
		
			if($bookhecont){

				print OUT "xrdcp -f $outdir/$ARGV[0]-$ARGV[1]/HeliumFragm/HeliumFragm$j.root /tmp/fdimicco/;\n";
				print OUT  "\$WORKDIR/HeliumContamination_Parallel \$WORKDIR/InputNtupleLists/FileListDT$j.txt \$WORKDIR/InputNtupleLists/FileListMC$j.txt  /tmp/fdimicco/HeliumFragm$j.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/fdimicco/HeliumFragm$j.root $outdir/$ARGV[0]-$ARGV[1]/HeliumFragm;\n";
			}
		
			if($bookeff){

				print OUT "xrdcp -f $outdir/$ARGV[0]-$ARGV[1]/MCEfficiency/MCEfficiency$j.root /tmp/fdimicco/;\n";
				print OUT  "\$WORKDIR/MCEfficiency_Parallel \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/MCEfficiency$j.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/fdimicco/MCEfficiency$j.root $outdir/$ARGV[0]-$ARGV[1]/MCEfficiency;\n";
			}
		
			if($bookcounts){

				print OUT "xrdcp -f $outdir/$ARGV[0]-$ARGV[1]/Counts/Counts$j.root /tmp/fdimicco/;\n";
				print OUT  "\$WORKDIR/CountsExtraction_Parallel \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/Counts$j.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/fdimicco/Counts$j.root $outdir/$ARGV[0]-$ARGV[1]/Counts;\n";
			}	
			
			if($bookeffcorr){

				print OUT "xrdcp -f $outdir/$ARGV[0]-$ARGV[1]/EffCorr/EffCorr$j.root /tmp/fdimicco/;\n";
				print OUT  "\$WORKDIR/EffCorr_Parallel \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/EffCorr$j.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/fdimicco/EffCorr$j.root $outdir/$ARGV[0]-$ARGV[1]/EffCorr;\n";
			}
		

			if($bookflux){

				print OUT "xrdcp -f $outdir/$ARGV[0]-$ARGV[1]/Fluxes/Fluxes$j.root /tmp/fdimicco/;\n";
				print OUT  "\$WORKDIR/Fluxes_Parallel \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/Fluxes$j.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/fdimicco/Fluxes$j.root $outdir/$ARGV[0]-$ARGV[1]/Fluxes;\n";
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

for($k=0;$k<4;$k++){
	$joblaunched = 0;
	system("kinit -R");
	$jobs = `bjobs -q $queue| wc -l`;
	while($jobs>120){
		$jobs = `bjobs -q $queue| wc -l`;
		sleep(10);
		$jobrunning = `bjobs -q $queue |grep RUN|wc -l`;
		print "Total Jobs launched: $joblaunched\n";
		print "Jobs currently running: $jobrunning\n";
	}
	while($joblaunched<$njobs){
		$jobrunning = `bjobs -q $queue |grep RUN|wc -l`;
		$jobs = `bjobs -q $queue| wc -l`;
		$check;
		$control=0;
		while($jobs<350 and $joblaunched<$njobs) {

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
			if($bookplots) { $check = `ls -la $outdir_plots/$ARGV[0]-$ARGV[1]/Plots/Plots$joblaunched.root| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($booklat) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/LatRew/LatRew$joblaunched.root| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($bookhecont) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/HeliumFragm/HeliumFragm$joblaunched.root| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($bookeff) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/MCEfficiency/MCEfficiency$joblaunched.root| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($bookcounts) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/Counts/Counts$joblaunched.root| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($bookeffcorr) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/EffCorr/EffCorr$joblaunched.root| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($bookflux) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/Fluxes/Fluxes$joblaunched.root| awk '{print\$5}'`; if($check<1000) {$control=1}} 

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



