#!/usr/bin/perl

use warnings;
chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
#chomp($workdir = "/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/");
print "Printed: Work Dir. = ".$workdir."\n\n";

$datapath  = "\\/eos\\/ams\\/group\\/dbar\\/release_v4\\/e1_vdev_180213\\/neg\\/ISS.B950\\/pass6\\/";
$mcP_path  = "\\/eos\\/ams\\/group\\/dbar\\/release_v4\\/e1_vdev_180213\\/full\\/Pr.B1082\\/pr.pl1.l1.05100.2_01\\/";
$mcD_path  = "\\/eos\\/ams\\/group\\/dbar\\/release_v4\\/e1_vdev_180213\\/full\\/D.B1081\\/d.pl1.l1.0_5200.2_01\\/";
$mcHe_path = "\\/eos\\/ams\\/group\\/dbar\\/release_v4\\/e1_vdev_180213\\/full\\/He.B1081\\/he4.pl1.l1.2200.2_01\\/";
$mcT_path  = "\\/eos\\/ams\\/group\\/dbar\\/release_v4\\/e1_vdev_180213\\/full\\/T.B1059\\/t.pl1.0_520_GG_BlicDPMJet\\/";

$out_path  = "\/eos\/ams\/user\/f\/fdimicco\/";

#$ntuplepath  = "/eos/ams/group/dbar/TrentoNTuples/$ARGV[0]-$ARGV[1]";
$ntuplepath  = "\\/eos\\/ams\\/user\\/f\\/fdimicco\\/AnalysisNTuples\\/$ARGV[0]-$ARGV[1]\\/Ntuples";



#use warnings;
if($ARGV[3]==0) { system("rm $workdir/InputFileLists/*"); }
else { system("rm $workdir/InputNtupleLists/*");}


print "Listing All Data Files..\n";
chomp (@Rootuple = `xrdfs root://eosams.cern.ch/ ls  $datapath | grep -v "log" |  sed s/.root//g  | sed 's/$datapath//g'`);
$num_Rootuple = scalar(@Rootuple);


print "xrdfs root://eosams.cern.ch/ ls  $ntuplepath | grep -v \"log\" |  sed s/.root//g| sed 's/$ntuplepath//g'";
chomp (@NTuple = `xrdfs root://eosams.cern.ch/ ls  $ntuplepath | grep -v "log" |  sed s/.root//g| sed 's/$ntuplepath//g'`);
$num_NTuple = scalar(@NTuple);

print "Total Files: ".$num_Rootuple."\n";
print "Total Ntuples: ".$num_NTuple."\n";


$datapath  =~ s/\\//g;
$ntuplepath  =~ s/\\//g;


@rootuple;
$i=0;

$njobs = $ARGV[2];

for($n=0;$n<$num_Rootuple;$n++){
	if($Rootuple[$n]>$ARGV[0]&& $Rootuple[$n]<$ARGV[1]){ 
		$rootuple[$i]=$Rootuple[$n];
		$i++;
	}
}

$num_rootuple = scalar(@rootuple);
print "Files in the requested period: ".$num_rootuple."\n";

for ($n=0;$n<$njobs; $n++)
{
	if($ARGV[3]==1){
		open(OUT,">","$workdir/InputNtupleLists/FileListDT$n.txt");
		for ($j=($num_NTuple)/$njobs*$n ; $j<($num_NTuple)/$njobs*($n+1) ; $j++){
			print OUT  "$ntuplepath/$NTuple[$j].root\n";
		}
	}
	else{
		open(OUT,">","$workdir/InputFileLists/FileListDT$n.txt");
		for ($j=($num_rootuple)/$njobs*$n ; $j<($num_rootuple)/$njobs*($n+1) ; $j++)
		{
			print OUT  "$datapath/$rootuple[$j].root\n";
		}
	}
}



print "Listing All MC Files..\n";
chomp (@MC_P = `xrdfs root://eosams.cern.ch/ ls  $mcP_path | grep -v "log" |  sed s/.root//g | sed 's/$mcP_path//g'`);
$num_MC_P = scalar(@MC_P);

print "Total Files MC P: ".$num_MC_P."\n";

chomp (@MC_D = `xrdfs root://eosams.cern.ch/ ls  $mcD_path | grep -v "log" |  sed s/.root//g | sed 's/$mcD_path//g'`);
$num_MC_D = scalar(@MC_D);

print "Total Files MC D: ".$num_MC_D."\n";

chomp (@MC_He = `xrdfs root://eosams.cern.ch/ ls $mcHe_path  | grep -v "log" |  sed s/.root//g | sed 's/$mcHe_path//g'`);
$num_MC_He = scalar(@MC_He);

print "Total Files MC He: ".$num_MC_He."\n";


chomp (@MC_T = `xrdfs root://eosams.cern.ch/ ls $mcT_path  | grep -v "log" |  sed s/.root//g  | sed 's/$mcT_path//g'`);
$num_MC_T = scalar(@MC_T);


print "Total Files MC T: ".$num_MC_T."\n";

$mcP_path  =~ s/\\//g;
$mcD_path  =~ s/\\//g;
$mcHe_path =~ s/\\//g;
$mcT_path  =~ s/\\//g;

for ($n=0;$n<$njobs; $n++)
{

	if($ARGV[3]==1){
		open(OUT,">","$workdir/InputNtupleLists/FileListMC$n.txt");
		for ($j=($num_NTuple)/$njobs*$n ; $j<($num_NTuple)/$njobs*($n+1) ; $j++){
			print OUT  "$ntuplepath/$NTuple[$j].root\n";
		}
	}
	else{
		open(OUT,">","$workdir/InputFileLists/FileListMC$n.txt");

		for ($j=($num_MC_P)/$njobs*$n ; $j<($num_MC_P)/$njobs*($n+1) ; $j++)
		{
			print OUT  "$mcP_path/$MC_P[$j].root\n"
		}
		for ($j=($num_MC_D)/$njobs*$n ; $j<($num_MC_D)/$njobs*($n+1) ; $j++)
		{
			print OUT  "$mcD_path/$MC_D[$j].root\n"
		}
		for ($j=($num_MC_He)/$njobs*$n ; $j<($num_MC_He)/$njobs*($n+1) ; $j++)
		{
			print OUT  "$mcHe_path/$MC_He[$j].root\n"
		}
		for ($j=($num_MC_T)/$njobs*$n ; $j<($num_MC_T)/$njobs*($n+1) ; $j++)
		{
			print OUT  "$mcT_path/$MC_T[$j].root\n"
		}


	}
}








