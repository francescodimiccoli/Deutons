#/usr/bin/perl

use warnings;
chomp (@LanciaDati = `ls | grep Dati`);
$num_Dati = scalar(@LanciaDati);
chomp($workdir =`pwd -P |sed 's\\Lancia\\\\g'`);

system("rm $workdir/logs/*.log");
system("rm $workdir/err/*.err");

if($ARGV[0]!=0&&$ARGV[2]!=1) {
	system("rm $workdir..//MAIN/sommadati/sommadati*");
	system("perl $workdir/../MAIN/perl/Sommaisto2.pl $ARGV[0] $ARGV[1]");
}


for ($n=0;$n<$num_Dati; $n++){
system("sh $LanciaDati[$n]");
}


