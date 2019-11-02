#!/usr/bin/perl

use warnings;

$njobs=40;

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
1700000000
);

$outdir="/eos/ams/user/f/fdimicco/AnalysisFiles";
$workdir="/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/";



$size = scalar(@bartels);
print "Time bins: ".scalar(@bartels)."\n";

$jobsrunning = 0;


$start=$ARGV[0];
$stop=$ARGV[1];

$launchana=1;
$launchsum=0;
$organizeoutput=0;
$finalanalysis=0;

for($i=$start;$i<$stop;$i++){

	if($launchana==1){
		print $bartels[$i]."\n";
		system("perl $workdir/perl/Create_condorscripts.pl $bartels[$i] $bartels[$i+1] $njobs 0");
		$jobsrunning = `condor_q| grep idle|awk '{print\$1}'`;
	}
	if($launchsum==1) {
		system("perl $workdir/perl/SumPartials.pl $bartels[$i] $bartels[$i+1] 10" );


	}
}

if($launchana==1){
	system("rm $workdir/perl/*.log");
	system("rm $workdir/perl/*.err");
	system("rm $workdir/perl/*.out");
		
	for($i=$start;$i<$stop;$i++){
		chdir "$outdir/$bartels[$i]-$bartels[$i+1]/Counts";
		system("condor_submit Condor_script.sub");
	}	
}

open(OUT,">","$workdir/perl/SumScripts/DoAllPartials.sh");

if($launchsum==1){
	for($i=$start;$i<$stop;$i++){
		chdir "$outdir/$bartels[$i]-$bartels[$i+1]";
			print $i."\n";
			print OUT "$workdir/perl/SumScripts/script$bartels[$i]-$bartels[$i+1].sh \$1\n";
		}	
}
close(OUT);

if($organizeoutput==1){
#	for($i=$start;$i<$stop;$i++){
#		chdir "$outdir/$bartels[$i]-$bartels[$i+1]";
#		system("mv $outdir/$bartels[$i]-$bartels[$i+1]/Counts/*_Flux $outdir/$bartels[$i]-$bartels[$i+1]/Flux");
#		system("mv $outdir/$bartels[$i]-$bartels[$i+1]/Counts/*_Corr $outdir/$bartels[$i]-$bartels[$i+1]/EffCorr");
#		system("hadd -f -k Result.root Partial*");
#
#	}
	for($i=$start;$i<=$stop;$i=$i+4){
#		system("hadd -f $outdir/Grouped/$bartels[$i]-$bartels[$i+4].root $outdir/$bartels[$i]-$bartels[$i+1]/Result.root  $outdir/$bartels[$i+1]-$bartels[$i+2]/Result.root $outdir/$bartels[$i+2]-$bartels[$i+3]/Result.root $outdir/$bartels[$i+3]-$bartels[$i+4]/Result.root");
	}
}

if($finalanalysis==1){
	for($i=$start;$i<=$stop;$i=$i+4){
	system("$workdir/Analysis  $workdir/AnalysisFiles/Grouped/$bartels[$i]-$bartels[$i+4].root >> $workdir/perl/AnalysisResults/$bartels[$i]-$bartels[$i+4].out 2> $workdir/perl/AnalysisResults/$bartels[$i]-$bartels[$i+4].err &");
	}
}

