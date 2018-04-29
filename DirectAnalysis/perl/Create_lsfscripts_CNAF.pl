#!/usr/bin/perl

use warnings;
print "Printing LSF scripts for Analysis...\n\n";
#chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
chomp($workdir = "/storage/gpfs_ams/ams/users/fdimicco/Deutons/DirectAnalysis/DirectAnalysis/");
print "Printed: Work Dir. = ".$workdir."\n\n";
$njobs=$ARGV[2];
$outdir="/storage/gpfs_ams/ams/users/fdimicco/Deutons/DirectAnalysis/AnalysisFiles";
#$outdir="/eos/ams/user/f/fdimicco/AnalysisNTuples";
#creating sum scripts and output directory
#

$queue = "ams1nd";

print "Resetting work directories...\n";

system ("$workdir/perl/SumInputFilesScripts.pl $ARGV[0] $ARGV[1] $ARGV[2] 0");
system ("$workdir/perl/SumInputFilesScripts.pl $ARGV[0] $ARGV[1] $ARGV[2] 1");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]");

system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/Ntuples");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/Plots");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/LatRew");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/HeliumFragm");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/MCEfficiency");		
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/Counts");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/Fluxes");

system ("rm -r $workdir/lsf/");
system ("mkdir $workdir/lsf/");
system ("rm -r $workdir/err/");
system ("mkdir $workdir/err/");	
system ("rm -r $outdir/../logs/$ARGV[0]-$ARGV[1]");
system ("mkdir $outdir/../logs/$ARGV[0]-$ARGV[1]");

print "Creating job files...";

$bookntuples=1;
$bookplots=0;
$booklat=0;
$bookhecont=0;
$bookeff=0;
$bookcounts=0;
$bookflux=0;

for($j=0;$j<$njobs;$j++)
{
		open(OUT,">","$workdir/lsf/lsf$j.tcsh");

		print OUT "#!/bin/bash

			export WORKDIR=$workdir;
			source \$WORKDIR/../amsvar_cvmfs.sh;\n";

			if($bookntuples){
				print OUT "xrdcp -f $outdir/$ARGV[0]-$ARGV[1]/Ntuples/Ntuple$j.root /tmp/;\n";
				print OUT  "\$WORKDIR/Ntuple_Maker \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/Ntuple$j.root 1 >> /tmp/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/Ntuple$j.root $outdir/$ARGV[0]-$ARGV[1]/Ntuples;\n";
			}
		
			if($bookplots){

				print OUT "xrdcp -f $outdir/$ARGV[0]-$ARGV[1]/Plots/Plots$j.root /tmp/;\n";
				print OUT  "\$WORKDIR/Distributions_Plotter \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/Plots$j.root 1 >> /tmp/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/Plots$j.root $outdir/$ARGV[0]-$ARGV[1]/Plots;\n";
			}               
		
			if($booklat){

				print OUT "xrdcp -f $outdir/$ARGV[0]-$ARGV[1]/LatRew/LatRew$j.root /tmp/;\n";
				print OUT  "\$WORKDIR/LatitudeReweighter \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/LatRew$j.root 1 >> /tmp/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/LatRew$j.root $outdir/$ARGV[0]-$ARGV[1]/LatRew;\n";
			}
		
			if($bookhecont){

				print OUT "xrdcp -f $outdir/$ARGV[0]-$ARGV[1]/HeliumFragm/HeliumFragm$j.root /tmp/;\n";
				print OUT  "\$WORKDIR/HeliumContamination_Parallel \$WORKDIR/InputNtupleLists/FileListDT$j.txt \$WORKDIR/InputNtupleLists/FileListMC$j.txt  /tmp/HeliumFragm$j.root 1 >> /tmp/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/HeliumFragm$j.root $outdir/$ARGV[0]-$ARGV[1]/HeliumFragm;\n";
			}
		
			if($bookeff){

				print OUT "xrdcp -f $outdir/$ARGV[0]-$ARGV[1]/MCEfficiency/MCEfficiency$j.root /tmp/;\n";
				print OUT  "\$WORKDIR/MCEfficiency_Parallel \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/MCEfficiency$j.root 1 >> /tmp/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/MCEfficiency$j.root $outdir/$ARGV[0]-$ARGV[1]/MCEfficiency;\n";
			}
		
			if($bookcounts){

				print OUT "xrdcp -f $outdir/$ARGV[0]-$ARGV[1]/Counts/Counts$j.root /tmp/;\n";
				print OUT  "\$WORKDIR/CountsExtraction_Parallel \$WORKDIR/InputNtupleLists/FileListDT$j.txt \$WORKDIR/InputNtupleLists/FileListMC$j.txt  /tmp/Counts$j.root 1 >> /tmp/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/Counts$j.root $outdir/$ARGV[0]-$ARGV[1]/Counts;\n";
			}	
			
			if($bookflux){

				print OUT "xrdcp -f $outdir/$ARGV[0]-$ARGV[1]/Fluxes/Fluxes$j.root /tmp/;\n";
				print OUT  "\$WORKDIR/Fluxes_Parallel \$WORKDIR/InputNtupleLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/Fluxes$j.root 1 >> /tmp/log$j.log;\n\n";
				print OUT "xrdcp -f /tmp/Fluxes$j.root $outdir/$ARGV[0]-$ARGV[1]/Fluxes;\n";
			}
			
			print OUT "xrdcp -f /tmp/log$j.log $outdir/../logs/$ARGV[0]-$ARGV[1]/;\n";


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

for($k=0;$k<3;$k++){
	$joblaunched = 0;

	while($joblaunched<$njobs){
		$jobrunning = `bjobs -q $queue |grep RUN|wc -l`;
		$jobs = `bjobs -q $queue| wc -l`;
		$check;
		$control=0;
		while($jobs<500 and $joblaunched<$njobs) {

			$control=0;
			if($bookntuples) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/Ntuples/Ntuple$joblaunched.root| awk '{print\$5}'`; if($check<1000) {print $check; $control=1}} 
			if($bookplots) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/Plots/Plots$joblaunched.root| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($booklat) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/LatRew/LatRew$joblaunched.root| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($bookhecont) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/HeliumFragm/HeliumFragm$joblaunched.root| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($bookeff) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/MCEfficiency/MCEfficiency$joblaunched.root| awk '{print\$5}'`; if($check<1000) {$control=1}} 
			if($bookcounts) { $check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/Counts/Counts$joblaunched.root| awk '{print\$5}'`; if($check<1000) {$control=1}} 
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



