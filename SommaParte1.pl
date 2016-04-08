#!/bin/perl

system("source ../amsvar.sh");
$workdir="/home/AMS/fdimicco/fdimicco/Deutons";
$mese = $ARGV[0];
system("rm $workdir/Histos/$mese/*tot*");
for($n=0;$n<10;$n++) {
system("mv $workdir/Histos/$mese/$mese\_$n\_P1.root $workdir/Histos/$mese/$mese\_0$n_\P1.root");
}

system("hadd -f $workdir/Histos/$mese/$mese\_tot_P1.root $workdir/Histos/$mese/$mese\_*;");
