#!/usr/bin/perl

use warnings;

$njobs=300;
$grouping = 4;

$start=$ARGV[0];
$stop=$ARGV[1];

$launchana=$ARGV[2];
$launchsum=$ARGV[2]-1;
$organizeoutput=$ARGV[2]-2;
$finalanalysis=$ARGV[2]-3;


@bartels = (
1307750400.0,
1310083200.0,
1312416000.0,
1314748800.0,
1317081600.0,
1319414400.0,
1321747200.0,
1324080000.0,
1326412800.0,
1328745600.0,
1331078400.0,
1333411200.0,
1335744000.0,
1338076800.0,
1340409600.0,
1342742400.0,
1345075200.0,
1347408000.0,
1349740800.0,
1352073600.0,
1354406400.0,
1356739200.0,
1359072000.0,
1361404800.0,
1363737600.0,
1366070400.0,
1368403200.0,
1370736000.0,
1373068800.0,
1375401600.0,
1377734400.0,
1380067200.0,
1382400000.0,
1384732800.0,
1387065600.0,
1389398400.0,
1391731200.0,
1394064000.0,
1396396800.0,
1398729600.0,
1401062400.0,
1403395200.0,
1405728000.0,
1408060800.0,
1410393600.0,
1412726400.0,
1415059200.0,
1417392000.0,
1419724800.0,
1422057600.0,
1424390400.0,
1426723200.0,
1429056000.0,
1431388800.0,
1433721600.0,
1436054400.0,
1438387200.0,
1440720000.0,
1443052800.0,
1445385600.0,
1447718400.0,
1450051200.0,
1452384000.0,
1454716800.0,
1457049600.0,
1459382400.0,
1461715200.0,
1464048000.0,
1466380800.0,
1468713600.0,
1471046400.0,
1473379200.0,
1475712000.0,
1478044800.0,
1480377600.0,
1482710400.0,
1485043200.0,
1487376000.0,
1489708800.0,
1492041600.0,
1494417600,
1496750400,
1499083200,
1501416000,
1503748800,
1506081600,
1508414400,
1510747200,
1513080000,
1515412800,
1517745600,
1520078400,
1522411200,
1524744000,
1527076800,
1529409600,
1531915200,
1534075200,
1536408000,
1538740800,
1541073600,
1543406400,
1545739200,
1548072000,
1550358000,

1552690800,
1555020000,
1557352800,
1559685600,
1562018400,
1564351200,
1566684000,
1569016800,
1571349600,
1573686000,
1576018800,
1578351600,
1580684400,
1583017200,
1585350000,
1587679200,
1590012000,
1592344800,
1594677600,
1597010400,
1599343200,
1601676000,
1604012400,
1606345200,
1608678000,
1611010800,
1613343600,
1615676400,
1618005600,
1620338400,
1622671200
);

#$outdir="/eos/ams/user/f/fdimicco/AnalysisFiles";
#$workdir="/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/";
$outdir="/storage/gpfs_ams/ams/users/fdimicco/Deutons/DirectAnalysis/AnalysisFiles";
$workdir="/storage/gpfs_ams/ams/users/fdimicco/Deutons/DirectAnalysis";



$size = scalar(@bartels);
print "Time bins: ".scalar(@bartels)."\n";

$jobsrunning = 0;

print $bartels[$start]." ". $bartels[$stop]."\n"; 

for($i=$start;$i<$stop;$i=$i+$grouping){

	if($launchana==1){
		print $bartels[$i]."\n";
		system("perl $workdir/perl/Create_condorscripts.pl $bartels[$i] $bartels[$i+$grouping] $njobs 0");
	}
	if($launchsum==1) {
		system("perl $workdir/perl/SumPartials.pl $bartels[$i] $bartels[$i+$grouping] 10" );


	}
}

if($launchana==1){
	system("rm $workdir/perl/*.log");
	system("rm $workdir/perl/*.err");
	system("rm $workdir/perl/*.out");
		
	for($i=$start;$i<$stop;$i=$i+$grouping){
		#chdir "$outdir/$bartels[$i]-$bartels[$i+$grouping]/Counts";
		system("mkdir $outdir/$bartels[$i]-$bartels[$i+$grouping]");
		system("condor_submit -name sn-02.cr.cnaf.infn.it $workdir/perl/AnalysisScripts/Condor_script$bartels[$i]-$bartels[$i+$grouping].sub");
	}	
}


if($organizeoutput==1){
	open(OUT,">","$workdir/perl/SumScripts/DoAllPartials.sh");
	for($i=$start;$i<$stop;$i=$i+$grouping){
		print OUT "hadd -f -k $outdir/$bartels[$i]-$bartels[$i+$grouping]/../Grouped/$bartels[$i]-$bartels[$i+$grouping].root  $outdir/$bartels[$i]-$bartels[$i+$grouping]/Partial* &"."\n";
	}
	close (OUT);
}


if($finalanalysis==1){
	for($i=$start;$i<=$stop;$i=$i+$grouping){

		open(OUT,">","$workdir/perl/AnalysisResults/Condor_Ana$bartels[$i]-$bartels[$i+$grouping].sub");
		print OUT "executable	= $workdir/Analysis \narguments	= $workdir/InputFileLists/$bartels[$i]-$bartels[$i+$grouping]/FileListDT12,txt $workdir/InputFileLists/$bartels[$i]-$bartels[$i+$grouping]/FileListMC12,txt $workdir/AnalysisFiles/Grouped/$bartels[$i]-$bartels[$i+$grouping].root \noutput	= $workdir/perl/AnalysisResults/$bartels[$i]-$bartels[$i+$grouping].out\nerror	= $workdir/perl/AnalysisResults/$bartels[$i]-$bartels[$i+$grouping].err\nlog	= $workdir/perl/AnalysisResults/$bartels[$i]-$bartels[$i+$grouping].log\n+JobFlavour = \"microcentury\"\nqueue 1"; 
		close (OUT);
#	chdir "$workdir/perl/SumScripts/";
#	system("condor_submit Condor_sum$ARGV[0]-$ARGV[1].sub");
#	chdir "$workdir/perl/";
	}
}



if($finalanalysis==1){

	open(OUT3,">","$workdir/perl/AnalysisResults/FinalFlux.sh");
	$fixed=36;
	for($i=$start;$i<=$stop;$i=$i+$grouping){

		print OUT3 "$workdir/Analysis $workdir/InputFileLists/$bartels[$i]-$bartels[$i+$grouping]/FileListDT0.txt $workdir/InputFileLists/$bartels[$fixed]-$bartels[$fixed+$grouping]/FileListMC0.txt $workdir/AnalysisFiles/Grouped/$bartels[$i]-$bartels[$i+$grouping].root\n";

		open(OUT,">","$workdir/perl/AnalysisResults/script$bartels[$i]-$bartels[$i+$grouping].sh");

		print OUT "#!/bin/bash

			export WORKDIR=$workdir;
		source \$WORKDIR/mysetenv.sh;\n";

		print OUT  "\$WORKDIR/Analysis \$WORKDIR/InputFileLists/$bartels[$i]-$bartels[$i+$grouping]/FileListDT0.txt \$WORKDIR/InputFileLists/$bartels[$fixed]-$bartels[$fixed+$grouping]/FileListMC0.txt $workdir/AnalysisFiles/Grouped/$bartels[$i]-$bartels[$i+$grouping].root  > $workdir/perl/AnalysisResults/$bartels[$i]-$bartels[$i+$grouping].out 2> $workdir/perl/AnalysisResults/$bartels[$i]-$bartels[$i+$grouping].err \n\n";

		close (OUT);
		
		system("chmod +x $workdir/perl/AnalysisResults/script$bartels[$i]-$bartels[$i+$grouping].sh");

		open(OUT2,">","$workdir/perl/AnalysisResults/Condor_Ana$bartels[$i]-$bartels[$i+$grouping].sub");
		print OUT2 "executable	= $workdir/perl/AnalysisResults/script$bartels[$i]-$bartels[$i+$grouping].sh \narguments	= \noutput	= $workdir/perl/AnalysisResults/$bartels[$i]-$bartels[$i+$grouping].out\nerror	= $workdir/perl/AnalysisResults/$bartels[$i]-$bartels[$i+$grouping].err\nlog	= $workdir/perl/AnalysisResults/$bartels[$i]-$bartels[$i+$grouping].log\n+JobFlavour = \"microcentury\"\nqueue 1"; 
		close (OUT2);
	chdir "$workdir/perl/AnalysisResults/";
	system("condor_submit $workdir/perl/AnalysisResults/Condor_Ana$bartels[$i]-$bartels[$i+$grouping].sub");
	chdir "$workdir/perl/";
	}
	close(OUT3)
}


