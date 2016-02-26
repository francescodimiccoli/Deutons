#/usr/bin/perl

use warnings;
chomp (@rootuple = `ls ../lsf `);
$numero=scalar(@rootuple);


for($j=0;$j<100;$j++){
print $j."\n" ;

open(OUT,">","/storage/gpfs_ams/ams/users/fdimicco/Lancia/LanciaDati$j.sh");

print OUT "#!/bin/bash

export WORKDIR=/storage/gpfs_ams/ams/users/fdimicco

chmod +x \$WORKDIR/lsf/lsf$j.tcsh
bsub -q ams -o \$WORKDIR/lsf/lsf$j.out -e \$WORKDIR/err/lsf$j.err \$WORKDIR/lsf/lsf$j.tcsh >>\$WORKDIR/lsf/lsf$j.log\n";

close (OUT);
}

