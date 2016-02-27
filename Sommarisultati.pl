#!/bin/perl

system("source ../amsvar.sh");
$workdir="/home/AMS/fdimicco/fdimicco";
$mese = $ARGV[0];
for($n=0;$n<10;$n++) {
system("mv $workdir/Risultati/$mese/RisultatiDATI_$n.root $workdir/Risultati/$mese/RisultatiDATI_0$n.root");
system("mv ./$mese/RisultatiMC_$n.root ./$mese/RisultatiMC_0$n.root");
}

system("hadd -f $workdir/Risultati/$mese/RisultatiDATI.root $workdir/Risultati/$mese/RisultatiDATI_*;");
system("hadd -f $workdir/Risultati/$mese/RisultatiMC.root $workdir/Risultati/$mese/RisultatiMC_*;");

#system("mv $workdir/Risultati/$mese/RisultatiDATI.root $workdir/Risultati/DATI-rootfile/$mese.root");
#system("mv $workdir/Risultati/$mese/ $workdir/Risultati/DATI-datafile");
