#/usr/bin/perl

use warnings;
print "Printing Launch scripts..\n\n";

print "Printing Launch scripts..\n\n";
chomp($workdir =`pwd -P |sed 's\\perl\\\\g'|sed 's\\Lancia\\\\g'`);
print "Printed: Work dir. = ".$workdir."\n\n";

for($j=0;$j<100;$j++){

open(OUT,">","../Lancia/LanciaDati$j.sh");

print OUT "#!/bin/bash

export WORKDIR=$workdir

chmod +x \$WORKDIR/lsf/lsf$j.tcsh
bsub -q ams -o \$WORKDIR/lsf/lsf$j.out -e \$WORKDIR/err/lsf$j.err \$WORKDIR/lsf/lsf$j.tcsh >>\$WORKDIR/lsf/lsf$j.log\n";

close (OUT);
}

