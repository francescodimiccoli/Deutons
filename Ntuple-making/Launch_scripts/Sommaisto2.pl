#usr/bin/perl


#use warnings;
system("rm /storage/gpfs_ams/ams/users/fdimicco/MAIN/SumScripts/Sommaisto*");

print "Listing All Files..\n";
chomp (@Rootuple = `ls  /storage/gpfs_ams/ams/users/fdimicco/MAIN/istogrammidatiCNAF | grep -v "log" |  sed s/.root//g`);
$num_Rootuple = scalar(@Rootuple);

print "Total Files: ".$num_Rootuple."\n";
@rootuple;
$i=0;

for($n=0;$n<$num_Rootuple;$n++){
	if($Rootuple[$n]>$ARGV[0]&& $Rootuple[$n]<$ARGV[1]){ 
		$rootuple[$i]=$Rootuple[$n];
		$i++;
	}
}

$num_rootuple = scalar(@rootuple);
print "Files in the requested period: ".$num_rootuple."\n";

for ($n=0;$n<100; $n++)
{

	open(OUT,">","/storage/gpfs_ams/ams/users/fdimicco/MAIN/SumScripts/Sommaisto$n.sh");
	print OUT  "#!/bin/bash

		hadd -f /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommadati/sommadati$n.root ";

	for ($j=($num_rootuple)/100*$n ; $j<($num_rootuple)/100*($n+1) ; $j++)
	{
				print OUT  " /storage/gpfs_ams/ams/users/fdimicco/MAIN/istogrammidatiCNAF/$rootuple[$j].root "
	}
	print OUT "\n";
}




for ($n=0;$n<100; $n++)
{

        open(OUT,">","/storage/gpfs_ams/ams/users/fdimicco/MAIN/SumScripts/SommaistoMC$n.sh");
        print OUT  "#!/bin/bash

                hadd -f /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/temp/sommaMC$n.root ";

        print OUT  " /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/L1MC/protons/sommaMC$n.root";
	print OUT  " /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/L1MC/deuterons/sommaMC$n.root";
	print OUT  " /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/L1MC/He/sommaMC$n.root";
	print OUT "\n";
}




open(OUT,">","/storage/gpfs_ams/ams/users/fdimicco/MAIN/SumScripts/Sommaisto.sh");

print OUT  "#!/bin/bash";
print OUT "\n";
for ($n=0;$n<100; $n++)
{
	print OUT  "sh /storage/gpfs_ams/ams/users/fdimicco/MAIN/SumScripts/Sommaisto$n.sh\n";
}










