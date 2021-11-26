#!/usr/bin/perl
#use warnings;


$workdir="/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/perl";
print "Printed: Work Dir. = ".$workdir."\n\n";


$outdir="/eos/ams/user/f/fdimicco/AnalysisFiles/";
print "Printed: Out Dir. = ".$outdir."\n\n";

$n = $ARGV[2];
print "Listing All Data Files..\n";
chomp (@RootupleCounts = `ls  $outdir/$ARGV[0]/ | grep -v "Result"| grep -v "Partial"|grep _Counts`);
chomp (@RootupleFlux = `ls  $outdir/$ARGV[0]/ 	 | grep -v "Result"| grep -v "Partial"|grep _Flux`);
chomp (@RootupleCorr = `ls  $outdir/$ARGV[0]/   | grep -v "Result"| grep -v "Partial"|grep _Corr`);
chomp (@RootupleLat  = `ls  $outdir/$ARGV[0]/   | grep -v "Result"| grep -v "Partial"|grep Lati`);


$num_Rootuplecounts = scalar(@RootupleCounts);
$num_Rootupleflux = scalar(@RootupleFlux);
$num_Rootuplecorr = scalar(@RootupleCorr);
$num_Rootuplelat = scalar(@RootupleLat);


print "Total Files Counts: ".$num_Rootuplecounts."\n";
print "Total Files Flux: ".$num_Rootupleflux."\n";
print "Total Files Eff: ".$num_Rootuplecorr."\n";
print "Total Files Lat: ".$num_Rootuplelat."\n";


@rootuple;$i=0;
$nparts = $ARGV[1];
$command;
$s="_";

system("source /cvmfs/sft.cern.ch/lcg/views/LCG_88/x86_64-slc6-gcc49-opt/setup.sh");


if(scalar(@RootupleCounts)>0){
	if($command eq"") { $command = "hadd -f -k -j $outdir/$ARGV[0]/Partial$n$s ";}
	$nsummed = $num_Rootuplecounts/$nparts;
	for($i=($nsummed)*$n;$i<(($nsummed)*($n+1));$i++){
		$command = $command." ".$outdir."/".$ARGV[0]."/".$RootupleCounts[$i];
	}
}

if(scalar(@RootupleFlux)>0){
	if($command eq "") { $command = "hadd -f -k -j $outdir/$ARGV[0]/Partial$n$s ";}
	$nsummed = $num_Rootupleflux/$nparts;
	for($i=($nsummed)*$n;$i<(($nsummed)*($n+1));$i++){
		$command = $command." ".$outdir."/".$ARGV[0]."/".$RootupleFlux[$i];
	}
}

if(scalar(@RootupleCorr)>0){
	if($command eq "") { $command = "hadd -f -k -j $outdir/$ARGV[0]/Partial$n$s ";}
	$nsummed = $num_Rootuplecorr/$nparts;
	for($i=($nsummed)*$n;$i<(($nsummed)*($n+1));$i++){
		$command = $command." ".$outdir."/".$ARGV[0]."/".$RootupleCorr[$i];
	}
}

if(scalar(@RootupleLat)>0){
	$command = "";
	if($command eq "") { $command = "hadd -f -k -j  $outdir/$ARGV[0]/Partial$n$s ";}
	$nsummed = $num_Rootuplelat/$nparts;
	for($i=($nsummed)*$n;$i<(($nsummed)*($n+1));$i++){
		$command = $command." ".$RootupleLat[$i];
	}
}

print $command;
system("$command");
