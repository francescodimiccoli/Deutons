#!/usr/bin/perl

use warnings;

$out_path  = "/eos/ams/user/f/fdimicco/";
$njobs = $ARGV[2];

for($n=0;$n<$njobs;$n++){
	$num =`ls  $out_path/AnalysisFiles/$ARGV[0]-$ARGV[1]/Analysis$n-*|wc -l`;
	if($num!=4 && $num!=0){
		print $n." job inconsistent: Only found ".$num."output files\n";
		system("rm $out_path/AnalysisFiles/$ARGV[0]-$ARGV[1]/Analysis$n-*|wc -l");
	}
}
