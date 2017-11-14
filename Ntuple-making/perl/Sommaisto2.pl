#usr/bin/perl


#use warnings;
system("rm /storage/gpfs_ams/ams/users/fdimicco/MAIN/SumScripts/Sommaisto*");

print "Listing All Data Files..\n";
chomp (@Rootuple = `ls  /storage/gpfs_ams/ams/users/fdimicco/MAIN/istogrammidatiCNAF | grep -v "log" |  sed s/.root//g`);
$num_Rootuple = scalar(@Rootuple);

print "Total Files: ".$num_Rootuple."\n";
@rootuple;
$i=0;
$njobs = 400;

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

	open(OUT,">","/storage/gpfs_ams/ams/users/fdimicco/MAIN/SumScripts/Sommaisto$n.sh");
	print OUT  "#!/bin/bash

		hadd -f /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommadati/sommadati$n.root ";

	for ($j=($num_rootuple)/$njobs*$n ; $j<($num_rootuple)/$njobs*($n+1) ; $j++)
	{
				print OUT  " /storage/gpfs_ams/ams/users/fdimicco/MAIN/istogrammidatiCNAF/$rootuple[$j].root "
	}
	print OUT "\n";
}



print "Listing All MC Files..\n";
chomp (@MC_P = `ls  /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/L1MC/protons | grep -v "log" |  sed s/.root//g`);
$num_MC_P = scalar(@MC_P);

print "Total Files MC P: ".$num_MC_P."\n";

chomp (@MC_D = `ls  /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/L1MC/deuterons | grep -v "log" |  sed s/.root//g`);
$num_MC_D = scalar(@MC_D);

print "Total Files MC D: ".$num_MC_D."\n";

chomp (@MC_He = `ls  /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/L1MC/He | grep -v "log" |  sed s/.root//g`);
$num_MC_He = scalar(@MC_He);

print "Total Files MC He: ".$num_MC_He."\n";



for ($n=0;$n<$njobs; $n++)
{

        open(OUT,">","/storage/gpfs_ams/ams/users/fdimicco/MAIN/SumScripts/SommaistoMC$n.sh");
        print OUT  "#!/bin/bash

        hadd -f /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/temp/sommaMC_P$n.root ";

	for ($j=($num_MC_P)/$njobs*$n ; $j<($num_MC_P)/$njobs*($n+1) ; $j++)
        {
                                print OUT  " /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/L1MC/protons/$MC_P[$j].root "
        }
        print OUT "\n

        hadd -f /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/temp/sommaMC_D$n.root ";

	for ($j=($num_MC_D)/$njobs*$n ; $j<($num_MC_D)/$njobs*($n+1) ; $j++)
        {
                                print OUT  " /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/L1MC/deuterons/$MC_D[$j].root "
        }
        print OUT "\n

	 hadd -f /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/temp/sommaMC_He$n.root ";

	for ($j=($num_MC_He)/$njobs*$n ; $j<($num_MC_He)/$njobs*($n+1) ; $j++)
        {
                                print OUT  " /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/L1MC/He/$MC_He[$j].root "
        }
        print OUT "\n

	hadd -f /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/temp/sommaMC$n.root ";
	print OUT "/storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/temp/sommaMC_P$n.root ";
	print OUT "/storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/temp/sommaMC_D$n.root ";
	print OUT "/storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/temp/sommaMC_He$n.root ";
	print OUT "\n";
	
}




open(OUT,">","/storage/gpfs_ams/ams/users/fdimicco/MAIN/SumScripts/Sommaisto.sh");

print OUT  "#!/bin/bash";
print OUT "\n";
for ($n=0;$n<$njobs; $n++)
{
	print OUT  "sh /storage/gpfs_ams/ams/users/fdimicco/MAIN/SumScripts/Sommaisto$n.sh\n";
}










