#/usr/bin/perl

use warnings;

system("./HeliumContamination dummy.root dummy.root $ARGV[0]");
system("./CountsExtraction dummy.root dummy.root $ARGV[0]");
system("./MCEfficiency dummy.root dummy.root $ARGV[0]");
