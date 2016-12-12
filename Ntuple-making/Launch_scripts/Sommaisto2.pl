#usr/bin/perl


#use warnings;
system("rm /storage/gpfs_ams/ams/users/fdimicco/MAIN/Sommaisto*");

print "Listing All Files..\n";
chomp (@Rootuple = `ls  /storage/gpfs_ams/ams/users/fdimicco/MAIN/istogrammidati | grep -v "log" |  sed s/.root//g`);
chomp (@Control =  `ls -la /storage/gpfs_ams/ams/users/fdimicco/MAIN/istogrammidati |grep "root" |grep -v "log" | awk '{print \$5}'`);
$num_Rootuple = scalar(@Rootuple);
$num_Contr = scalar(@Control);

print "File Totali: ".$num_Rootuple." ".$num_Contr."\n";
@rootuple;
$i=0;

for($n=0;$n<$num_Rootuple;$n++){
	if($Rootuple[$n]>$ARGV[0]&& $Rootuple[$n]<$ARGV[1]){ 
		$rootuple[$i]=$Rootuple[$n];
		$i++;
	}
}

$num_rootuple = scalar(@rootuple);
$controllo;
print "File nel periodo richiesto: ".$num_rootuple."\n";

for ($n=0;$n<100; $n++)
{

	open(OUT,">","/storage/gpfs_ams/ams/users/fdimicco/MAIN/Sommaisto$n.sh");
	print OUT  "#!/bin/bash

		hadd -f /storage/gpfs_ams/ams/users/fdimicco/MAIN/sommadati/sommadati$n.root ";

	for ($j=($num_rootuple)/100*$n ; $j<($num_rootuple)/100*($n+1) ; $j++)
	{
		if(($rootuple[$j]!=$controllo)){
				print OUT  " /storage/gpfs_ams/ams/users/fdimicco/MAIN/istogrammidati/$rootuple[$j].root "
		}
		$controllo=$rootuple[$j];
	}
	print OUT "\n";
}

open(OUT,">","/storage/gpfs_ams/ams/users/fdimicco/MAIN/Sommaisto.sh");

print OUT  "#!/bin/bash";
print OUT "\n";
for ($n=0;$n<100; $n++)
{
	print OUT  "sh /storage/gpfs_ams/ams/users/fdimicco/MAIN//Sommaisto$n.sh\n";
}










