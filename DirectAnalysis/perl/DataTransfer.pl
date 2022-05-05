#!/usr/bin/perl

use warnings;
chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
#chomp($workdir = "/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/");
print "Printed: Work Dir. = ".$workdir."\n\n";

#$datapath  = "/eos/ams/group/dbar/release_v6/e2_vdev_190525/full/Be.B1200/be7.pl1.l1.48000.4_00";
$datapath = "\\/eos\\/ams\\/group\\/dbar\\/release_v7\\/e1_vdev_200421\\/neg\\/ISS.B1130\\/pass7\\/";

$out_path  = "\\/eos\\/ams\\/group\\/dbar\\/TrentoNTuples\\/v7_pass7_Q2filtered";


#use warnings;
system("rm $workdir/perl/TransferScripts/*"); 

print "Listing All Data Files..\n";

chomp (@Rootuple = `xrdfs root://eosams.cern.ch/ ls $datapath | grep root | grep -v "log" |  sed s/.root//g | sed 's/$datapath//g'`);
$num_Rootuple = scalar(@Rootuple);

print $num_Rootuple."\n";

#for($n=0;$n<$num_Rootuple;$n++){
#	print $Rootuple[$n]."\n";
#}

print "Total Files: ".$num_Rootuple."\n";
@rootuple;
$i=0;

$njobs = $ARGV[2];

for($n=0;$n<$num_Rootuple;$n=$n+1){
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
		print OUT "root -q -b /afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/perl/filtertree.C\\(\\\"$datapath/$rootuple[$j].root\\\",\\\"$out_path/$rootuple[$j].root\\\",\\\"Event\\\",\\\"q_inn\\\[0\\\]\\\",\\\"2.50\\\"\\);\n";
#		print OUT "scp $out_path/$rootuple[$j].root fdimicco\@ams.tifpa.infn.it:/data1/home/data/v7_pass7/Q2filtered\n";
	}

}

open(OUT2,">","$workdir/perl/TransferScripts/allscripts.sh");

for ($n=0;$n<$njobs; $n++)
{
	open(OUT,">","$workdir/perl/TransferScripts/Condor_script$n.sub");
	print OUT "executable	= $workdir/perl/TransferScripts/script$n.sh\narguments	=\noutput	= script$n.out\nerror	= script$n.err\nlog	= script$n.log\n+JobFlavour = \"workday\"\nqueue" 
}

for ($n=0;$n<$njobs; $n++)
{
	system("chmod +x $workdir/perl/TransferScripts/script$n.sh");
	system( "condor_submit $workdir/perl/TransferScripts/Condor_script$n.sub\n");
}



