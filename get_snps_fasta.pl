# Kate Lee Oct 2016
# get_snps_fasta.pl
# read in annotated and filtered vcf snp info
# takes in filelist
# outputs data for affy
# usage: perl get_snps_fasta.pl <filtered_vcf> <filelist> <outfile>
#
# filelist should be a list of hihi contig files, with a second column denoting the assembly
############################################################################################
#!usr/bin/perl -w
use strict; 


#unless (length(@ARGV) == 3) { die "Needs three inputs <filtered_vcf> <filelist> <outfile>\n filelist must contain filenames and assembly numbers on each line sperated by a tab \n";}

# open filtered vcf file
open(VCF, "<$ARGV[0]") || die "ERROR: couldn't open filtered vcf file: $!";

# open contig filelist
open (FILELIST, "<$ARGV[1]")  || die "ERROR: couldn't open filelist: $!";

# open outfile
open(OUT, ">$ARGV[2]") || die "ERROR: couldn't open outfile: $!";


# put contig filelist into a hash
my %filehash;
while(<FILELIST>) {
	chomp;
	my @row = split("\t", $_);
	$filehash{$row[1]} = $row[0];
}

# loop through filtered vcf file to get list of sequences
my %seqlist;
my @values = qw|chrom pos ref alt|;
my @row;
print "looping through vcf to get positional and assembly info...\n";
while(<VCF>) {
	chomp;
	@row = split("\t", $_);
	if ($row[0] eq '#CHROM') { 
	    $row[0] = 'CHROM';
	    unshift (@row, '#snp_seq');
	    my $writeline = join("\t", @row);
	    print OUT "$writeline\n";
	}
	# seqlist{assembly}{chromosome}{position}{ref/alt}
#	if ($_ =~ /RAD/){
#	    print $_."\n";
#	    $seqlist{$row[82]}{$row[0]}{$row[1]}{'ref'} = $row[3];
#	    $seqlist{$row[82]}{$row[0]}{$row[1]}{'alt'} = $row[4];
#	}
	else {
	$seqlist{$row[82]}{$row[0]}{$row[1]}{'ref'} = $row[3];
	$seqlist{$row[82]}{$row[0]}{$row[1]}{'alt'} = $row[4];
	my $keepline = join("\t", @row);
	$seqlist{$row[82]}{$row[0]}{$row[1]}{'line'} = $keepline;
	#print "row going into hash: $keepline\n\n";
	}
#	print "$row[61], $row[82]\n"
}

for my $assembly (sort keys %seqlist){
    print "processing sequences from assembly $assembly...\n";
    
    # open assembly fasta file and read in sequence info
    open(DB, "$filehash{$assembly}") || die "ERROR: couldn't open assembly fasta for $filehash{$assembly}: $!\n";
    my $seqname;
    my %dbhash = ();
    while (<DB>)
    {
	my $line = $_;
	chomp $line;
	if ($_= /^>([^\s]+)/){
	    #$dbhash{$1} = $line;
	    $seqname = $1;
	    next;
	}
	else {$dbhash{$seqname} .= $line;}
    }
    
    # loop through vcf info and grab fasta sequence (snp plus 36 bp either side)
    for my $chromosome (sort keys %{$seqlist{$assembly}}) {
	print "printing for chr: $chromosome\n";
	for my $position (sort keys %{$seqlist{$assembly}{$chromosome}}){
	   # print "pos: $position\n";
	    my $ref = $seqlist{$assembly}{$chromosome}{$position}{'ref'};
	    my $alt = $seqlist{$assembly}{$chromosome}{$position}{'alt'};
	   # my $rowref = $seqlist{$assembly}{$chromosome}{$position}{'line'};
	   # my $row = @{$rowref};
	    my @line = split("\t", $seqlist{$assembly}{$chromosome}{$position}{'line'});
	    
	   # print "line: @line\n\n";
	   # print "row: @row\n\n";
	    my $seqstart = $position + 35 ;
	    my $seqend = $position - 36;
	    
	   # if ( ( substr($dbhash{$chromosome},($position - 1) ,1) ) eq $ref) { print "ref ok!! $chromosome $position $ref $alt \n"; }
	    my $fastaseq =   substr($dbhash{$chromosome},($position -36),35)."[$ref/$alt]".substr($dbhash{$chromosome},($position),35);
	    unshift(@line, $fastaseq);
	    my $printrow = join("\t", @line);
	    print OUT $printrow."\n";

	}
	
    }
}

