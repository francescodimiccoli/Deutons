#!/usr/bin/perl

chomp($workdir =`pwd -P `);
print "Printed: Work Dir. = ".$workdir."\n\n";

#use warnings;
system("mv $workdir/$ARGV[0]/Result.root_Results PreviousResults.root");
system("rm $workdir/$ARGV[0]/Result*");

print "Listing All Data Files..\n";
chomp (@Rootuple = `ls  $workdir/$ARGV[0]/* |grep .root| grep -v "_Results"|grep -v "Check" `);
$num_Rootuple = scalar(@Rootuple);

print "Total Files: ".$num_Rootuple."\n";
@rootuple;
$i=0;
$nparts = 100;
$nsummed = $num_Rootuple/$nparts;
$command;

system("rm $workdir/$ARGV[0]/Result*");
	
for ($n=0;$n<$nparts; $n++)
{
	$command = "hadd -f -k $workdir/$ARGV[0]/Result_P$n ";
	for($i=($nsummed)*$n;$i<(($nsummed)*($n+1));$i++){
		$command = $command." ".$workdir."/".$ARGV[0]."/*/".$Rootuple[$i];
	}
	print $command."\n";
	print "\n";
	system("bsub -q ams1nd $command");	
}

#system("hadd -f $workdir/$ARGV[0]/Result.root $workdir/$ARGV[0]/Result_P*");

$jobs = `bjobs -q ams1nd| wc -l`;
$running = `bjobs -q ams1nd|grep RUN| wc -l`;
$missing = `ls $workdir/$ARGV[0]/Result_P*|wc -l`;
$missing = $nparts - $missing;

while($jobs>0){
	
	$jobs = `bjobs -q ams1nd| wc -l`;
	$running = `bjobs -q ams1nd|grep RUN| wc -l`;
	$missing = `ls $workdir/$ARGV[0]/Result_P*|wc -l`;
	$missing = $nparts - $missing;

	print "results missing: ".$missing."\n";
	print "jobs: ".$jobs."\n";
	print "running: ".$running."\n";

	sleep 10;
}

system ("perl SumPartials.pl $ARGV[0]");

$jobs = `bjobs -q ams1nd| wc -l`;
$running = `bjobs -q ams1nd|grep RUN| wc -l`;
$missing = ` ls $workdir/$ARGV[0]/Result_T*|wc -l`;
$missing = 10 - $missing;


while($jobs>0){
	$jobs = `bjobs -q ams1nd| wc -l`;
	$running = `bjobs -q ams1nd|grep RUN| wc -l`;
	$missing = ` ls $workdir/$ARGV[0]/Result_P*|wc -l`;
	$missing = $nparts - $missing;

	print "results missing: ".$missing."\n";
	print "jobs: ".$jobs."\n";
	print "running: ".$running."\n";

	sleep 10;
}

system("hadd -f -k $ARGV[0]/Result.root $ARGV[0]/Result_T*");


