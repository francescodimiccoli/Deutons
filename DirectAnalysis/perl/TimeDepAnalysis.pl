#!/usr/bin/perl

use warnings;

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
);

$outdir="/data1/home/fdimicco/Deutons/DirectAnalysis/AnalysisFiles";


$size = scalar(@bartels);
print "Time bins: ".scalar(@bartels)."\n";

$jobsrunning = 0;

for($i=10;$i<25;$i++){

	print $bartels[$i]."\n";
	system("perl Launch_all.pl $bartels[$i] $bartels[$i+1] 40 0");
	system ("cd $outdir/..");
	system("sh $outdir/../Jobs.tcsh");
	system ("cd $outdir/../perl");
	$jobsrunning = `ps ux| grep "Analysis "|grep -v "grep"|wc -l`;
	while($jobsrunning > 1) {
		$jobsrunning = `ps ux| grep Analysis|grep -v "grep"|wc -l`;
		print "bartel nr. $i: Jobs running: $jobsrunning\n";
		sleep(5);
	}	
#	system("hadd -f $outdir/$bartels[$i]-$bartels[$i+1]/Counts/Partial_Counts.root $outdir/$bartels[$i]-$bartels[$i+1]/*/Analysis*Counts*");
#	system("hadd -f $outdir/$bartels[$i]-$bartels[$i+1]/Fluxes/Partial_Flux.root   $outdir/$bartels[$i]-$bartels[$i+1]/*/Analysis*Flux*");
#	system("hadd -f $outdir/$bartels[$i]-$bartels[$i+1]/Result.root   $outdir/$bartels[$i]-$bartels[$i+1]/*/Partial_*");
#	system("$outdir/../Analysis  $outdir/$bartels[$i]-$bartels[$i+1]/Result.root");
}

