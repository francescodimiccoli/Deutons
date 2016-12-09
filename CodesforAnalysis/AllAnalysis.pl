#usr/bin/perl


use warnings;
$mode =  $ARGV[0]//"0";
$month = $ARGV[1]//"ALL";

system("rm ../Histos/TOT/*");
system("cp  ../Histos/*/*tot*  ../Histos/TOT");

if($ARGV[1]=="ALL"){
	chomp (@Rootuple = `ls  ../Histos/TOT/|grep tot  |  sed s/_tot_P1.root//g`);

}
if($ARGV[1]!="ALL"){
	@Rootuple = $ARGV[1];
}

$num_Rootuple = scalar(@Rootuple);

open(OUT,">","./Mesi.h");
print $num_Rootuple."\n";

print OUT "int num_mesi = $num_Rootuple;\n";
print OUT "std::string mesi[]{\n";
print OUT "\"$Rootuple[$num_Rootuple-1]\",\n";
for($i=0;$i<$num_Rootuple-1;$i++) {
	print OUT "\"$Rootuple[$i]\",\n";
}
print OUT "};\n";
if($ARGV[0]=="0"){
	for($i=0;$i<$num_Rootuple;$i++) {
		system("rm ./Final_plots/$Rootuple[$i].root");
		print "Analyzing ".$Rootuple[$i]."...\n";
		system("./Analysis $Rootuple[$i] 2 tot >out");
	}
}










           
