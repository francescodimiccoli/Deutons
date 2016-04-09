#usr/bin/perl


use warnings;

#$percorso="/home/francesco/PhD/LocalCNAF/";
#$percorso="/storage/gpfs_ams/ams/users/fdimicco/Deutons";
system("cp  ../Histos/*/*tot*  ../Histos/TOT");
chomp (@Rootuple = `ls  ../Histos/TOT/|grep tot  |  sed s/_tot_P1.root//g`);
$num_Rootuple = scalar(@Rootuple);
open(OUT,">","./Mesi.h");
print $num_Rootuple."\n";

print OUT "int num_mesi = $num_Rootuple;\n";
print OUT "std::string mesi[]{\n";
for($i=0;$i<$num_Rootuple;$i++) {
	print OUT "\"$Rootuple[$i]\",\n";
}
print OUT "};\n";

for($i=0;$i<$num_Rootuple;$i++) {
	system("rm ./Final_plots/$Rootuple[$i].root");
	print "Analyzing ".$Rootuple[$i]."...\n";
	system("./Analysis $Rootuple[$i] 2 tot");
	system("cp  ../Histos/$Rootuple[$i]/*tot*  ../Histos/TOT");
}











           
