#usr/bin/perl

chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
print "Printed: Work Dir. = ".$workdir."\n\n";

$datapath  = "/eos/ams/group/dbar/release_1507300526/neg/ISS.B950/pass6/";
$mcP_path  = "/eos/ams/group/dbar/release_1507300526/full/Pr.B1116/pr.pl1.l1.054000.3_00/";
$mcD_path  = "/eos/ams/group/dbar/release_1507300526/full/D.B1081/d.pl1.l1.0_5200.2_01/";
$mcHe_path = "/eos/ams/group/dbar/release_1507300526/full/He.B1116/he4.pl1.l1.24000.3_02/";
$out_path  = "/eos/ams/user/f/fdimicco/";


#use warnings;
system("rm $workdir/InputFileLists/*");

print "Listing All Data Files..\n";
chomp (@Rootuple = `ls  $datapath | grep -v "log" |  sed s/.root//g`);
$num_Rootuple = scalar(@Rootuple);

print "Total Files: ".$num_Rootuple."\n";
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
	open(OUT,">","$workdir/InputFileLists/FileListDT$n.txt");

	for ($j=($num_rootuple)/$njobs*$n ; $j<($num_rootuple)/$njobs*($n+1) ; $j++)
	{
		print OUT  "$datapath/$rootuple[$j].root\n";
	}
}



print "Listing All MC Files..\n";
chomp (@MC_P = `ls  $mcP_path | grep -v "log" |  sed s/.root//g`);
$num_MC_P = scalar(@MC_P);

print "Total Files MC P: ".$num_MC_P."\n";

chomp (@MC_D = `ls  $mcD_path | grep -v "log" |  sed s/.root//g`);
$num_MC_D = scalar(@MC_D);

print "Total Files MC D: ".$num_MC_D."\n";

chomp (@MC_He = `ls $mcHe_path  | grep -v "log" |  sed s/.root//g`);
$num_MC_He = scalar(@MC_He);

print "Total Files MC He: ".$num_MC_He."\n";



for ($n=0;$n<$njobs; $n++)
{

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
	
}








