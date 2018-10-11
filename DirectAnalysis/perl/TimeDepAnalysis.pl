#!/usr/bin/perl

use warnings;


@Times = (1306886400,1383264000,1446336000,1509494400);
#@Times = {1306886400,1383264000,1446336000,1509494400};

$size = scalar(@Times);
print "Time bins: ".scalar(@Times)."\n";


for($i=2;$i<$size;$i++){

	print $Times[$i]."\n";
	system ("perl Create_lsfscripts.pl $Times[$i-1]   $Times[$i] 1000");	

}

