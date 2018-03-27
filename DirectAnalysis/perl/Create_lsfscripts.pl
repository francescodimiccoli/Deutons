#/usr/bin/perl

use warnings;
print "Printing LSF scripts for Analysis...\n\n";
chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
print "Printed: Work Dir. = ".$workdir."\n\n";
$njobs=$ARGV[2];
$outdir="/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/AnalysisFiles";
#$outdir="/eos/ams/user/f/fdimicco/AnalysisNTuples";
#creating sum scripts and output directory
#

$queue = "ams1nd";

print "Resetting work directories...\n";

system ("perl SumInputFilesScripts.pl $ARGV[0] $ARGV[1] $ARGV[2] 0");
system ("perl SumInputFilesScripts.pl $ARGV[0] $ARGV[1] $ARGV[2] 1");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]");
system ("rm -r $workdir/lsf/");
system ("mkdir $workdir/lsf/");
system ("rm -r $workdir/err/");
system ("mkdir $workdir/err/");	
system ("rm -r $outdir/../logs/$ARGV[0]-$ARGV[1]");
system ("mkdir $outdir/../logs/$ARGV[0]-$ARGV[1]");

print "Creating job files...";

for($j=0;$j<$njobs;$j++)
{
		open(OUT,">","$workdir/lsf/lsf$j.tcsh");

		print OUT "#!/bin/bash

			export WORKDIR=$workdir;
			source \$WORKDIR/../amsvar_cvmfs.sh;\n";

		print OUT "xrdcp $outdir/$ARGV[0]-$ARGV[1]/Result$j.root /tmp/fdimicco/;\n";
		print OUT "xrdcp $outdir/../logs/$ARGV[0]-$ARGV[1]/Result$j.root /tmp/fdimicco/;\n";

		print OUT  "\$WORKDIR/Ntuple_Maker \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/Result$j.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
#		print OUT  "\$WORKDIR/LatitudeReweighter \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt   $outdir/$ARGV[0]-$ARGV[1]/Result$j.root 1 >> \$WORKDIR/logs/$ARGV[0]-$ARGV[1]/log$j.log;\n\n";



#		print OUT  "\$WORKDIR/Distributions_Plotter \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  /tmp/fdimicco/Result$j.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
#		print OUT  "\$WORKDIR/HeliumContamination_Parallel \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt   $outdir/$ARGV[0]-$ARGV[1]/Result$j.root 1 >> \$WORKDIR/logs/$ARGV[0]-$ARGV[1]/log$j.log;\n\n";
		#print OUT  "\$WORKDIR/CountsExtraction_Parallel \$WORKDIR/InputNtupleLists/FileListDT$j.txt \$WORKDIR/InputNtupleLists/FileListMC$j.txt  /tmp/fdimicco/Result$j.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
#		print OUT  "\$WORKDIR/MCEfficiency_Parallel \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt   /tmp/fdimicco/Result$j.root 1 >> /tmp/fdimicco/log$j.log;\n\n";
#		print OUT  "\$WORKDIR/Distributions_Plotter \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt   $outdir/$ARGV[0]-$ARGV[1]/Result$j.root 1 >> \$WORKDIR/logs/$ARGV[0]-$ARGV[1]/log$j.log;\n\n";
#		print OUT  "\$WORKDIR/Fluxes_Parallel \$WORKDIR/InputNtupleLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt   $outdir/$ARGV[0]-$ARGV[1]/Result$j.root 1 >> \$WORKDIR/logs/$ARGV[0]-$ARGV[1]/log$j.log;\n\n";


		print OUT "xrdcp -f /tmp/fdimicco/Result$j.root $outdir/$ARGV[0]-$ARGV[1]/;\n";
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
while($joblaunched<$njobs){
	$jobrunning = `bjobs -q $queue |grep RUN|wc -l`;
	$jobs = `bjobs -q $queue| wc -l`;
	while($jobs<500 and $joblaunched<$njobs) {
			system("bsub -q $queue -o $workdir/lsf/lsf$joblaunched.out -e $workdir/err/lsf$joblaunched.err $workdir/lsf/lsf$joblaunched.tcsh\n");
			$joblaunched++;	
			$jobs = `bjobs -q $queue|wc -l`;
			print $jobs."\n";
	}	
	print "Total Jobs launched: $joblaunched\n";
	print "Jobs currently running: $jobrunning\n";
	sleep(10);
}



$joblaunched = 0;

while($joblaunched<$njobs){
        $jobrunning = `bjobs -q $queue |grep RUN|wc -l`;
        $jobs = `bjobs -q $queue| wc -l`;
 	$check;
        while($jobs<500 and $joblaunched<$njobs) {
               		$check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/Result$joblaunched.root| awk '{print\$5}'`;
			if($check<1000){
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

$joblaunched = 0;

while($joblaunched<$njobs){
        $jobrunning = `bjobs -q $queue |grep RUN|wc -l`;
        $jobs = `bjobs -q $queue| wc -l`;
 	$check;
        while($jobs<500 and $joblaunched<$njobs) {
               		$check = `ls -la $outdir/$ARGV[0]-$ARGV[1]/Result$joblaunched.root| awk '{print\$5}'`;
			if($check<1000){
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



