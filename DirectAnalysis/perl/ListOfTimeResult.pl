#!/usr/bin/perl

use warnings;
print "Listing available result files...\n\n";
#chomp($workdir =`pwd -P |sed 's\\perl\\\\g '`);
chomp($workdir = "/data1/home/fdimicco/Deutons/DirectAnalysis");
print "Printed: Work Dir. = ".$workdir."\n\n";

chomp (@files = `ls $workdir/AnalysisFiles/1*/Result.root_Results`);

for($i=0;$i<scalar(@files);$i++) {
	print "$files[$i]\n";
}

open(OUT,">","$workdir/perl/List.h");

print OUT "#include <string.h>\n";
print OUT "#include <TVector.h>\n";

print OUT "std::vector<std::string> TimeFiles { \n";

for($i=0;$i<scalar(@files);$i++) {
	if($i==0) {
		print OUT "\" $files[$i]\"";
	}
	print OUT ",\n\" $files[$i]\"";
}

print OUT "};\n\n";


chomp (@files_group = `ls $workdir/AnalysisFiles/Grouped/1*_Results`);

for($i=0;$i<scalar(@files_group);$i++) {
	print "$files_group[$i]\n";
}


print OUT "std::vector<std::string> GroupedFiles { \n";

for($i=0;$i<scalar(@files_group);$i++) {
	if($i==0) {
		print OUT "\" $files_group[$i]\"";
	}
	print OUT ",\n\" $files_group[$i]\"";
}

print OUT "};";


