$workdir=$ARGV[0];
$inizio =$ARGV[1]//0;
$fine   =$ARGV[2]//120000;
$n	=$inizio;

if( !(-d "$workdir/Launch" ) ) { die "The scripts folder $workdir/Launch doesnt exist. Please create." }

chomp (@rootuple = `more eos_data.txt`);

while($n<=$fine){
    $jobs = `bjobs|wc -l`;
    $jobsrun = `bjobs|grep RUN|wc -l`;
	while($jobs<500){
        $jobs = `bjobs|wc -l`;
        $jobsrun = `bjobs|grep RUN|wc -l`;
		system("sh $workdir/Launch/$rootuple[$n].sh");
		$n++;
	}
	print "jobs launched: ".$n."\n";
	print "tot jobs: ".$jobs." running: ".$jobsrun."\n";
    system("sleep 10");

}






