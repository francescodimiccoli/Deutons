#!/usr/bin/perl

use warnings;
chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
#chomp($workdir = "/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/");
print "Printed: Work Dir. = ".$workdir."\n\n";

#$datapath  = "/eos/ams/group/dbar/release_v6/e2_vdev_190525/neg/ISS.B1130/pass7/";
$datapath  = "/eos/ams/group/dbar/release_v6/e2_vdev_190525/full/Pr.B1200/pr.pl1.05100.4_00/";


#$out_path  = "/eos/ams/group/dbar/TrentoNTuples/v6_pass7/data";
$out_path  = "/eos/ams/group/dbar/TrentoNTuples/v6_pass7/MC";





#use warnings;
system("rm $workdir/perl/CuttingScripts/*"); 

print "Listing All Data Files..\n";

chomp (@Rootuple = `ls $datapath | grep root | grep -v "log"`);
$num_Rootuple = scalar(@Rootuple);


for($n=0;$n<$num_Rootuple;$n++){
	print $Rootuple[$n]."\n";
}

print "Total Files: ".$num_Rootuple."\n";
@rootuple;
$i=0;

$njobs = $ARGV[0];


for ($n=0;$n<$njobs; $n++)
{
	open(OUT,">","$workdir/perl/CuttingScripts/script$n.sh");
	print OUT "#!/bin/bash

			export WORKDIR=$workdir;
			source /cvmfs/sft.cern.ch/lcg/views/LCG_88/x86_64-slc6-gcc49-opt/setup.sh;\n";

	for ($j=($num_Rootuple)/$njobs*$n ; $j<($num_Rootuple)/$njobs*($n+1) ; $j++){
	print OUT "rootcp $datapath/$Rootuple[$j]:Compact $datapath/$Rootuple[$j]:RTI $out_path/$Rootuple[$j]\n";
	}

#	print OUT "scp ";
#	for ($j=($num_Rootuple)/$njobs*$n ; $j<($num_Rootuple)/$njobs*($n+1) ; $j++){
#	print OUT "$out_path/$Rootuple[$j]  ";
#	}	
#	print OUT "fdimicco\@ams.tifpa.infn.it://data1/home/data/v6_pass7/Compact_data\n";

}

open(OUT2,">","$workdir/perl/CuttingScripts/allscripts.sh");

for ($n=0;$n<$njobs; $n++)
{
	open(OUT,">","$workdir/perl/CuttingScripts/Condor_script$n.sub");
	print OUT "executable	= script$n.sh\narguments	=\noutput	= script$n.out\nerror	= script$n.err\nlog	= script$n.log\n+JobFlavour = \"workday\"\nqueue" 
}

for ($n=0;$n<$njobs; $n++)
{
	system("chmod +x $workdir/perl/CuttingScripts/script$n.sh");
	print OUT2 "condor_submit Condor_script$n.sub\n";
}



