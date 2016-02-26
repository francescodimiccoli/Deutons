#/usr/bin/perl

$inizio=7;
$fine=8;
$n=$inizio;
$secondi=0;
print $jobs;
$workdir="/home/AMS/fdimicco/fdimicco";
@date=(1306886400,1309478400,1312156800,1314835200,1317427200,1320105600,1322697600,1325376000,1328054400,1330560000,1333238400,1335830400,1338508800,1341100800,1343779200,1346457600,1349049600,1351728000,1354320000,1356998400,1359676800,1362096000,1364774400,1367366400,1370044800,1372636800,1375315200,1377993600,1380585600,1383264000,1385856000);
@mesi=("2011_05","2011_06","2011_07","2011_08","2011_09","2011_10","2011_11","2011_12","2012_01","2012_02","2012_03","2012_04","2012_05","2012_06","2012_07","2012_08","2012_09","2012_10","2012_11","2012_12","2013_01","2013_02","2013_03","2013_04","2013_05","2013_06","2013_07","2013_08","2013_09","2013_10");
while($n<=$fine){
	$jobs = `bjobs|wc -l`;
	$jobsrun = `bjobs|grep RUN|wc -l`;
	if($jobs==0){
		if($ARGV[0]==0||$ARGV[0]==1){
			$time = `date`;
			print $jobs." ".$time;
			$secondi=0;
			print "creating $mesi[$n] directory...";
			if($n!=$inizio) {system("rm -r $workdir/Risultati/$mesi[$n]");}
			system("mv $workdir/Risultati/risultati $workdir/Risultati/$mesi[$n]");
			print "created\n";
			print "creating new risultati directory...";
			system("mkdir $workdir/Risultati/risultati");
			print "created\n";
			if($n!=$fine){
				system("perl $workdir/perl/lsf.pl $mesi[$n+1] 0");
				system("perl $workdir/Lancia/Lancia.pl $date[$n] $date[$n+1] 0");
				system("bjobs|wc -l");
				print "Sto Lanciando $mesi[$n+1] : da $date[$n] a $date[$n+1]";
			}
		}
		if($ARGV[0]==3){
			if($n!=$fine){
				system("perl $workdir/perl/lsf.pl $mesi[$n+1] 1");
                                system("perl $workdir/Lancia/Lancia.pl $date[$n] $date[$n+1] 1");
                                system("bjobs|wc -l");
                                print "Sto Lanciando $mesi[$n+1] : da $date[$n] a $date[$n+1]";				
			}
		}
		if($n!=$inizio) {
			if($ARGV[0]==1){
				{system("perl $workdir/Sommarisultati.pl $mesi[$n]");
			    	system("$workdir/CodesforAnalysis/Preliminar.exe $mesi[$n]");
				}
			}
			if($ARGV[0]==2||$ARGV[0]==0||$ARGV[0]==3) {system("perl $workdir/SommaParte1.pl $mesi[$n]");}	
		}	
		$n++;
	}
	system("sleep 20");
	print $jobs." jobs: ".$jobsrun." running...";
	if($jobs>0&&$jobs<5){ $secondi++;}
	if($secondi>90) {system("bkill -u fdimicco 0");}

}

