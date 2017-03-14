#!/usr/bin/perl


use warnings;
use File::Basename;
use Cwd 'abs_path';

my $dirname = dirname(abs_path(__FILE__));


my $outputdir = $ARGV[0];


if( !(-d $outputdir) ){ die "The folder $outputdir does not exist. Please create." }
if( !(-w $outputdir) ){ die "Cannot write to $outputdir, please check your permissions."}

system("mkdir -p $outputdir/lsf");
system("mkdir -p $outputdir/logs");
system("mkdir -p $outputdir/err");
system("mkdir -p $outputdir/Launch");

system ("sh Create_eoslist.sh ");
chomp (@rootuple = `more eos_data.txt`);


$num_rootuple = scalar(@rootuple);

print "Tot. number of files: ".$num_rootuple."\n";

for($j=0;$j<$num_rootuple;$j++)
{
	open(OUT,">","$outputdir/lsf/$rootuple[$j].sh");
	print OUT "#!/bin/bash
export WORKDIR=$dirname/../ 
export INPUTPATH=/storage/gpfs_ams/ams/Rec/2014/ISS.B950/pass6/
export WORKSPACE=$outputdir
source \$WORKDIR/amsvar_cvmfs.sh
source \$WORKDIR/amsvar.lxplus.sh	

LD_LIBRARY_PATH=\$WORKDIR/:\$LD_LIBRARY_PATH

\$WORKDIR/dstmake \$INPUTPATH/$rootuple[$j] \$WORKSPACE/$rootuple[$j] >> \$WORKSPACE/logs/$rootuple[$j].log\n
		  ";
	close (OUT);
}
for($j=0;$j<$num_rootuple;$j++)
{
	open(OUT2,">","$outputdir/Launch/$rootuple[$j].sh");	
	print OUT2 "#!/bin/bash
export WORKDIR=$dirname/../
export WORKSPACE=$outputdir
chmod +x \$WORKSPACE/lsf/$rootuple[$j].sh
	   
bsub -q ams -e \$WORKSPACE/err/$rootuple[$j].err \$WORKSPACE/lsf/$rootuple[$j].sh ";

	close (OUT2);

}


