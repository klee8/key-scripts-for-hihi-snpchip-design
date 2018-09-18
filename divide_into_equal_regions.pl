# Kate Lee, May 2016
# divide_into_equal_regions.pl
#
# script to take in a fasta file of draft genome from soap-de-novo
# outputs a list of genome regions of approximately equal length for use in GATK
#
# Calculates total length of genome
# Divide for number of jobs [default:50]
# Output stats - total length, number of regions, approximate time for jobs
#
# use: perl get_chr_regions.pl <fastafile.fa> <optional: #output files>
#
#####################################################################################

#!usr/bin/perl -w
use strict;


my $fastafile = $ARGV[0];
my $num_out = $ARGV[1] || 50;


#open fasta file
open (FASTA, "<$fastafile") || die "ERROR, couldn't open fasta file: $!";


# Loop through fasta to get headers and info
my %headers;
my $length = 0;
my $count = 0;

while(<FASTA>){
    if ($_ =~ /^>(\d+)\slength\s(\d+)\scvg/) {
	if ($2 >= 200) {
	    $headers{$1}=$2;
	    $length = $length + $2;
	    $count++;
	}
    }
}

# decide optimal interval size
my $interval = $length/$num_out;


# make output dir
`mkdir interval_files`;
`mkdir interval_bed`;

# create output files
my $counter = 1;
my $length = 0;
my $list = "";
my $bedlist = "";
my $chr;

for $chr (sort keys %headers) {
    if  ($length + $headers{$chr} > $interval) {
	open (OUT, ">interval_files/regions_$counter.intervals") || die "ERROR, couldn't open interval file $counter: $!";
	open (BED, ">interval_bed/regions_$counter.bed") || die "ERROR, couldn't open interval files $counter: $!";
	print OUT "$list$chr\n";
	print BED "$bedlist$chr\t0\t$headers{$chr}\n";
	$counter++;
	$length = 0;
	$list = "";
	$bedlist = "";
	close OUT;
	close BED;
    }
    else { 
	$length = $length + $headers{$chr};
	$bedlist = "$bedlist$chr\t0\t$headers{$chr}\n";
	$list = "$list$chr\n";
    }
}

# get last set of regions
open (OUT, ">interval_files/regions_$counter.intervals") || die "ERROR, couldn't open interval file $counter:$!";
print OUT "$list";
close OUT;
open (BED, ">interval_bed/regions_$counter.bed") || die "ERROR, couldn't open interval files $counter: $!";
print BED "$bedlist";
close BED;
exit;









