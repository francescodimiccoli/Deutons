#!/usr/bin/perl

use warnings;
chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
#chomp($workdir = "/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/");
print "Printed: Work Dir. = ".$workdir."\n\n";

#$datapath  = "/eos/ams/group/dbar/release_v6/e2_vdev_190525/full/Be.B1200/be7.pl1.l1.48000.4_00";
$datapath = "/eos/ams/group/dbar/release_v6/e2_vdev_190525/neg/ISS.B1130/pass7/";

$out_path  = "/eos/ams/group/dbar/TrentoNTuples/v6_pass7/data";




#use warnings;
system("rm $workdir/perl/TransferScripts/*"); 

print "Listing All Data Files..\n";

chomp (@Rootuple = `ls $datapath | grep root | grep -v "log" |  sed s/.root//g`);
$num_Rootuple = scalar(@Rootuple);


for($n=0;$n<$num_Rootuple;$n++){
	print $Rootuple[$n]."\n";
}

print "Total Files: ".$num_Rootuple."\n";
@rootuple;
$i=0;

$njobs = $ARGV[2];

for($n=0;$n<$num_Rootuple;$n=$n+2){
	if($Rootuple[$n]>$ARGV[0]&& $Rootuple[$n]<$ARGV[1]){ 
		$rootuple[$i]=$Rootuple[$n];
		$i++;
	}
}

$num_rootuple = scalar(@rootuple);
print "Files in the requested period: ".$num_rootuple."\n";



for ($n=0;$n<$njobs; $n++)
{
	open(OUT,">","$workdir/perl/TransferScripts/script$n.sh");
	print OUT "#!/bin/bash

			export WORKDIR=$workdir;
			source \$WORKDIR/mysetenv.sh;\n";

	for ($j=($num_rootuple)/$njobs*$n ; $j<($num_rootuple)/$njobs*($n+1) ; $j++){
#		print OUT "root -q -b /afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/perl/filtertree.C\\(\\\"$datapath/$rootuple[$j].root\\\",\\\"/eos/ams/group/dbar/TrentoNTuples/FilteredQ2/$rootuple[$j].root\\\",\\\"Event\\\",\\\"q_inn\\\[0\\\]\\\",\\\"2.50\\\"\\);\n";
		print OUT "scp /eos/ams/group/dbar/TrentoNTuples/FilteredQ2/$rootuple[$j].root fdimicco\@ams.tifpa.infn.it:/data1/home/data/prova\n";
	}

}

open(OUT2,">","$workdir/perl/TransferScripts/allscripts.sh");

for ($n=0;$n<$njobs; $n++)
{
	open(OUT,">","$workdir/perl/TransferScripts/Condor_script$n.sub");
	print OUT "executable	= script$n.sh\narguments	=\noutput	= script$n.out\nerror	= script$n.err\nlog	= script$n.log\n+JobFlavour = \"workday\"\nqueue" 
}

for ($n=0;$n<$njobs; $n++)
{
	system("chmod +x $workdir/perl/TransferScripts/script$n.sh");
	print OUT2 "condor_submit Condor_script$n.sub\n";
}



