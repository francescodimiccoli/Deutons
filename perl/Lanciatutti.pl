#/usr/bin/perl
#
use warnings;
open(OUT,">","/storage/gpfs_ams/ams/users/fdimicco//Lancia/Lancia.sh");

#dati dal 4 maggio (06:00) al 5 maggio (10:00) 2013
#for($j=9553;$j<9753;$j++){
#print OUT "sh Lancia_Code$j.sh\n"

#dati dal 15 giugno (06:00) al 17 luglio (10:00) 2013
#for($j=20000;$j<20700;$j++){
#print OUT "sh Lancia_Code$j.sh\n"

for($j=0;$j<100;$j++){
print OUT "sh Lancia$j.sh\n"

}
