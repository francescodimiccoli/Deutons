#/usr/bin/perl

use warnings;

system("./HeliumContamination_Parallel dummy.root dummy.root $ARGV[0]");
system("./CountsExtraction_Parallel dummy.root dummy.root $ARGV[0]");
system("./MCEfficiency_Parallel dummy.root dummy.root $ARGV[0]");
system("./EffCorr_Parallel dummy.root dummy.root $ARGV[0]");
system("./Fluxes_Parallel dummy.root dummy.root $ARGV[0]");

