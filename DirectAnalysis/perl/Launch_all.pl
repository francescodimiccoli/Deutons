#!/usr/bin/perl

use warnings;
print "Launching codes for Analysis...\n\n";
#chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
chomp($workdir = "/data1/home/fdimicco/Deutons/DirectAnalysis/");
print "Printed: Work Dir. = ".$workdir."\n\n";
$njobs=$ARGV[2];
$outdir="/data1/home/fdimicco/Deutons/DirectAnalysis/AnalysisFiles";
$outdir_ntuples="/data1/home/fdimicco/Deutons/DirectAnalysis/AnalysisNTuples";
$outdir_plots="/data1/home/fdimicco/Deutons/DirectAnalysis/AnalysisFiles";
#creating sum scripts and output directory

$offset = 0;

print "Resetting work directories...\n";

system ("$workdir/perl/SumInputFilesScripts_TIFPA.pl $ARGV[0] $ARGV[1] $ARGV[2] $offset");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]");
system ("mkdir $outdir_ntuples/$ARGV[0]-$ARGV[1]");
system ("mkdir $outdir_plots/$ARGV[0]-$ARGV[1]");


system ("mkdir $outdir_ntuples/$ARGV[0]-$ARGV[1]/Ntuples");
system ("mkdir $outdir_plots/$ARGV[0]-$ARGV[1]/Plots");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/LatRew");
system ("mkdir $outdir/$ARGV[0]-$ARGV[1]/HeliumFragm");

system ("rm -r $workdir/logs/");
system ("mkdir $workdir/logs/");

system ("rm -r $workdir/errs/");
system ("mkdir $workdir/errs/");



$bookntuples=0;
$bookplots=1;
$booklat=0;
$bookhecont=0;
$bookflux=0;
$bookinfos=0;


open(OUT,">","$workdir/Jobs.tcsh");

		print OUT "#!/bin/bash

		export WORKDIR=$workdir;
		source /data1/home/work/software/root.6.14/bin/thisroot.sh;
		export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/data1/home/work/data_reduction/v5/lib;\n";



for($j=0;$j<$njobs;$j++)
{
			if($bookinfos){
			}
		

			if($bookntuples){
			}
		
			if($bookplots){
				print OUT  "\$WORKDIR/Distributions_Plotter \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListDT$j.txt \$WORKDIR/InputFileLists/$ARGV[0]-$ARGV[1]/FileListMC$j.txt  $outdir/$ARGV[0]-$ARGV[1]/Plots/Plots$j-$offset.root 1  > \$WORKDIR/logs/log$j.log 2>\$WORKDIR/errs/err$j.err &\n\n";
	
			}               
		
			if($booklat){
			}
		
			if($bookhecont){
				print OUT  "\$WORKDIR/HeliumContamination_Parallel \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  $outdir/$ARGV[0]-$ARGV[1]/HeliumFragm/HeFragm$j-$offset.root 1  > \$WORKDIR/logs/log$j.log 2>\$WORKDIR/errs/err$j.err &\n\n";
			}
		

			if($bookflux){

				#print OUT  "\$WORKDIR/Analysis \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  $outdir/$ARGV[0]-$ARGV[1]/Analysis$j-$offset.root 1  > \$WORKDIR/logs/log$j.log 2>\$WORKDIR/errs/err$j.err &\n\n";
				print OUT  "\$WORKDIR/FastCheckSmearing \$WORKDIR/InputFileLists/FileListDT$j.txt \$WORKDIR/InputFileLists/FileListMC$j.txt  $outdir/$ARGV[0]-$ARGV[1]/FastSmearing$j-$offset.root 1  > \$WORKDIR/logs/log$j.log 2>\$WORKDIR/errs/err$j.err &\n\n";
			}
			

}

close (OUT);





