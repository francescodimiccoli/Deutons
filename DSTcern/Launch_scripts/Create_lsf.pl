#!/usr/bin/perl


use warnings;

system ("sh Create_eoslist.sh");
chomp (@rootuple = `more eos_data.txt`);


$num_rootuple = scalar(@rootuple);

print "Tot. number of files: ".$num_rootuple."\n";

for($j=0;$j<$num_rootuple;$j++)
{
	open(OUT,">","/storage/gpfs_ams/ams/users/fdimicco/Deutons/DSTcern/lsf/$rootuple[$j].sh");
	print OUT "#!/bin/bash
export WORKDIR=/storage/gpfs_ams/ams/users/fdimicco/Deutons/DSTcern 
export INPUTPATH=/storage/gpfs_ams/ams/Rec/2014/ISS.B950/pass6/
export WORKSPACE=/storage/gpfs_ams/ams/users/fdimicco/MAIN/istogrammidatiCNAF
source \$WORKDIR/amsvar_cvmfs.sh
source \$WORKDIR/amsvar.lxplus.sh	

LD_LIBRARY_PATH=\$WORKDIR/:\$LD_LIBRARY_PATH

\$WORKDIR/dstmake \$INPUTPATH/$rootuple[$j] \$WORKSPACE/$rootuple[$j] >> \$WORKDIR/logs/$rootuple[$j].log\n
		  ";
	close (OUT);
}
for($j=0;$j<$num_rootuple;$j++)
{
	open(OUT2,">","/storage/gpfs_ams/ams/users/fdimicco/Deutons/DSTcern/Launch/$rootuple[$j].sh");	
	print OUT2 "#!/bin/bash
export WORKDIR=/storage/gpfs_ams/ams/users/fdimicco/Deutons/DSTcern 
chmod +x \$WORKDIR/lsf/$rootuple[$j].sh
	   
bsub -q ams -e \$WORKDIR/err/$rootuple[$j].err \$WORKDIR/lsf/$rootuple[$j].sh ";

	close (OUT2);

}


