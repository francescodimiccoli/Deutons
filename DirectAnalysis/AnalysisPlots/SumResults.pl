#!/usr/bin/perl

chomp($workdir =`pwd -P `);
print "Printed: Work Dir. = ".$workdir."\n\n";

#use warnings;
system("rm $workdir/$ARGV[0]/Result_P*");

print "Listing All Data Files..\n";
chomp (@Rootuple = `ls  $workdir/$ARGV[0] | grep -v "_Results"`);
$num_Rootuple = scalar(@Rootuple);

print "Total Files: ".$num_Rootuple."\n";
@rootuple;
$i=0;
$nparts = $num_Rootuple/10;
$command;

system("rm $workdir/$ARGV[0]/Result_P*");
system("rm $workdir/$ARGV[0]/Result.root");
	
for ($n=0;$n<$nparts; $n++)
{
	$command = "hadd -f $workdir/$ARGV[0]/Result_P$n ";
	for($i=($num_Rootuple/$nparts)*$n;$i<(($num_Rootuple/$nparts)*($n+1));$i++){
		$command = $command." ".$workdir."/".$ARGV[0]."/".$Rootuple[$i];
	}
	system($command);	
}

system("hadd -f $workdir/$ARGV[0]/Result.root $workdir/$ARGV[0]/Result_P*");







