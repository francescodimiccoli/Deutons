#/usr/bin/perl

use warnings;

system("./CountsExtraction $ARGV[0] $ARGV[1]");
system("./MCEfficiency $ARGV[0] $ARGV[1]");
system("./P_Fluxes $ARGV[0] $ARGV[1]");
system("./D_Fluxes  $ARGV[0] $ARGV[1]");
system("./RICHEffCorr $ARGV[0] $ARGV[1]");
system("./QualEffCorr $ARGV[0] $ARGV[1]");
