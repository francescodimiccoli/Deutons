#!/usr/bin/perl

chomp($workdir =`pwd -P `);
print "Printed: Work Dir. = ".$workdir."\n\n";

#use warnings;
system("rm $workdir/$ARGV[0]/Result_T*");

print "Listing All Data Files..\n";
chomp (@Rootuple = `ls  $workdir/$ARGV[0]/Counts | grep -v "_Results"`);
$num_Rootuple = scalar(@Rootuple);

print "Total Files: ".$num_Rootuple."\n";
@rootuple;
$i=0;
$nparts = 10;
$nsummed = $num_Rootuple/$nparts;
$command;

system("rm $workdir/$ARGV[0]/Result.root");
	
for ($n=0;$n<$nparts; $n++)
{
	$command = "hadd -f -k $workdir/$ARGV[0]/Result_T$n ";
	for($i=($nsummed)*$n;$i<(($nsummed)*($n+1));$i++){
		$command = $command." ".$workdir."/".$ARGV[0]."/Counts/".$Rootuple[$i];
	}
	print $command."\n\n";
	system("$command &");	
}








