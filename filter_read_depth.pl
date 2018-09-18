#!usr/bin/perl -w 
use strict;

open(IN, "<$ARGV[0]") || die "ERROR: cannot open $ARGV[0]: $1";
open(OUT, ">$ARGV[1]") || die "ERROR: cannot open $ARGV[1]: $!";

my @row;
my $read_depth = $ARGV[2];
my $quality = $ARGV[3];

while(<IN>){
    chomp;
    if ($_ =~ /^#/) { 
	print OUT $_."\n";
    }
    else {
	@row = split("\t", $_);
	unless ($row[5] < $read_depth){ 
	    foreach my $value (split(";", $row[7])){ 
		my @temp = split("=", $value );
		if ($temp[0] eq 'DP'){
		    if ($temp[1] >= $quality  ) {
			print OUT  $_."\n";
		    }
		}
	    }

	}
    }
}
