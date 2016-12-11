$inizio =$ARGV[0]//0;
$fine   =$ARGV[1]//120000;
$n	=$inizio;

chomp (@rootuple = `more eos_data.txt`);

while($n<=$fine){
	$jobs = `bjobs|wc -l`;
        $jobsrun = `bjobs|grep RUN|wc -l`;
	while($jobs<500){
        	$jobs = `bjobs|wc -l`;
        	$jobsrun = `bjobs|grep RUN|wc -l`;
		system("sh ../Launch/$rootuple[$n].sh");
		$n++;
	}
	print "jobs launched: ".$n."\n";
	print "tot jobs: ".$jobs." running: ".$jobsrun."\n";
        system("sleep 10");

}






