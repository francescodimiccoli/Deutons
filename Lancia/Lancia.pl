#/usr/bin/perl

use warnings;
chomp (@LanciaDati = `ls | grep Dati`);
$num_Dati = scalar(@LanciaDati);
system("rm /home/AMS/fdimicco/fdimicco/fdimicco/Deutons/logs/*.log");
system("rm /home/AMS/fdimicco/fdimicco/fdimicco/Deutons/err/*.err");
if($ARGV[0]!=0&&$ARGV[2]!=1) {
system("rm /home/AMS/fdimicco/fdimicco/fdimicco/MAIN/sommadati/sommadati*");
system("perl /home/AMS/fdimicco/fdimicco/fdimicco/MAIN/perl/Sommaisto2.pl $ARGV[0] $ARGV[1]");
}


for ($n=0;$n<$num_Dati; $n++){
system("sh $LanciaDati[$n]");
}


