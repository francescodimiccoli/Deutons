#!/usr/bin/perl

use warnings;


@Times = (1305944557,1337566957,1369102957,1400638957,1432174957,1463797357,1495333357);

$size = scalar(@Times);
print "Time bins: ".scalar(@Times)."\n";


for($i=3;$i<$size;$i++){

	print $Times[$i]."\n";
	system ("perl Create_lsfscripts.pl $Times[$i]   $Times[$i+1] 1000");	

}

