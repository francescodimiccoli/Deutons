#/usr/bin/perl

use warnings;



@date=(1306886400,1309478400,1312156800,1314835200,1317427200,1320105600,1322697600,1325376000,1328054400,1330560000,1333238400,1335830400,1338508800,1341100800,1343779200,1346457600,1349049600,1351728000,1354320000,1356998400,1359676800,1362096000,1364774400,1367366400,1370044800,1372636800,1375315200,1377993600,1380585600,1383264000,1385856000);

$monthsdone=0;

for($i=$ARGV[0];$i<$ARGV[1];$i++){

	$jobs = `bjobs|wc -l`;
	$jobsrun = `bjobs|grep RUN|wc -l`;
	$seconds=0;
	#if($jobs<1){
	#	system("perl Create_lsfscripts.pl $date[$i] $date[$i+1]");
	#}
	if($monthsdone>0) { system("hadd -f ../Calibfiles/$date[$i-1]/Calib.root ../Calibfiles/$date[$i-1]/Calib_data*"); }
	while($jobs>0){
		$jobs = `bjobs|wc -l`;
		$jobsrun = `bjobs|grep RUN|wc -l`;
		print "tot jobs: ".$jobs." running: ".$jobsrun."\n";
		system("sleep 20");
		if($jobs>0&&$jobs<5){ $seconds++;}
		if($seconds>90) {system("bkill 0");}
	}
	$monthsdone++;
}

