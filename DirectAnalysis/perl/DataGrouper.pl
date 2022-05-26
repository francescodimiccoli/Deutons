#!/usr/bin/perl

use warnings;
chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
#chomp($workdir = "/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/");
print "Printed: Work Dir. = ".$workdir."\n\n";

#$datapath  = "/eos/ams/group/dbar/release_v6/e2_vdev_190525/full/Be.B1200/be7.pl1.l1.48000.4_00";
$datapath =  "\\/eos\\/ams\\/group\\/dbar\\/TrentoNTuples\\/v7_pass7_Q2filtered\\/";
$out_path  = "\\/eos\\/ams\\/group\\/dbar\\/TrentoNTuples\\/v7_pass7_Q2filtered\\/bartels\\/";

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


#use warnings;
system("rm $workdir/perl/GroupingScripts/*"); 

print "Listing All Data Files..\n";

chomp (@Rootuple = `xrdfs root://eosams.cern.ch/ ls $datapath | grep root | grep -v "log" |  sed s/.root//g | sed 's/$datapath//g'`);
$num_Rootuple = scalar(@Rootuple);

print $num_Rootuple."\n";

#for($n=0;$n<$num_Rootuple;$n++){
#	print $Rootuple[$n]."\n";
#}

print "Total Files: ".$num_Rootuple."\n";

$njobs = scalar(@bartels);



for ($n=0;$n<$njobs; $n++)
{
	@rootuple;
	$i=0;

	for($m=0;$m<$num_Rootuple;$m=$m+1){
		if($Rootuple[$m]>$bartels[$n]&& $Rootuple[$m]<$bartels[$n+1]){ 
			$rootuple[$i]=$Rootuple[$m];
			$i++;
		}
	}

	$num_rootuple = scalar(@rootuple);
	print "Files in the requested period: ".$num_rootuple."\n";


	open(OUT,">","$workdir/perl/GroupingScripts/script$n.sh");
	print OUT "#!/bin/bash

		export WORKDIR=$workdir;
	source \$WORKDIR/mysetenv.sh;\n
	hadd -f -k $out_path$rootuple[0]_1.root";

	for ($j=1; $j<($num_rootuple) ; $j=$j+2){
		print OUT " $datapath$rootuple[$j].root";
#		print OUT "scp $out_path/$rootuple[$j].root fdimicco\@ams.tifpa.infn.it:/data1/home/data/v7_pass7/Q2filtered\n";
	}

	print OUT "\n

		export WORKDIR=$workdir;
	source \$WORKDIR/mysetenv.sh;\n
	hadd -f -k $out_path$rootuple[0]_2.root";

	for ($j=0; $j<($num_rootuple) ; $j=$j+2){
		print OUT " $datapath$rootuple[$j].root";
#		print OUT "scp $out_path/$rootuple[$j].root fdimicco\@ams.tifpa.infn.it:/data1/home/data/v7_pass7/Q2filtered\n";
	}



}

open(OUT2,">","$workdir/perl/GroupingScripts/allscripts.sh");

for ($n=0;$n<$njobs; $n++)
{
	open(OUT,">","$workdir/perl/GroupingScripts/Condor_script$n.sub");
	print OUT "executable	= $workdir/perl/GroupingScripts/script$n.sh\narguments	=\noutput	= script$n.out\nerror	= script$n.err\nlog	= script$n.log\n+JobFlavour = \"workday\"\nqueue" 
}

for ($n=0;$n<$njobs; $n++)
{
	system("chmod +x $workdir/perl/GroupingScripts/script$n.sh");
	system( "condor_submit $workdir/perl/GroupingScripts/Condor_script$n.sub\n");
}



