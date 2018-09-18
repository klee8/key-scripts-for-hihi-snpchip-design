# Kate Lee Sept 2016
# Blast parser.pl
# parse a blast file and get positions of contigs for vcf annotation
####################################################################

#!usr/bin/perl -w
use strict; 
use Storable;

my $blast_file = $ARGV[0];
my $blastfile;
# open the blast file
open($blastfile, "<$blast_file") || die "ERROR, couldn't open $blast_file:$!";

open(OUT, ">hashfile.txt")  || die "ERROR, couldn't open hashfile.txt:$!";  

my $linecount = `wc -l $blast_file`;

# read in whole file and chunk by contig hits
my $contig_name="lalala";
my $contig_match_length;
my $old_contig_name = "start";
my $db_chr = "";  # database chromosome
my $db_chr_counter = 0; # counts number of chromosomes a query hits
my $counter = 0; # number of hits
my @blastresults;
my $length_disparity_counter; 
my %blasthash; # hash for output info (and easy searching)
my $linecounter = 0;
my $lineflag = 0;

# ADD header line to hashfile for ref
@{$blasthash{"headers"}} = qw|qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore num_query_hits num_chrs_hit|; 

# read in whole file and chunk by contig hits
# sort and pass best hit to hash for each query
while(<$blastfile>){
    $linecounter++;
    if ( ( ($linecounter >= $linecount*0.2) && ($lineflag == 0) ) ||  ( ($linecounter >= $linecount*0.4) && ($lineflag == 20) ) ||  ( ($linecounter > $linecount*0.4) && ($lineflag == 40) ) ||  ( ($linecounter >= $linecount*0.6) && ($lineflag == 60) ) ||  ( ($linecounter >= $linecount*0.8) && ($lineflag == 80) ) ) { 
	$lineflag = $lineflag + 20;
	print "parsed $lineflag\% of blastfile...\n"
    }
    # capture cutsomised blast output (fmt 6 plust query length [qlen])
    # qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore
    chomp;
    my @blastline = split("\t", $_);
    $contig_name = $blastline[0];

    ### first line of file, or new hits of the same query
    if ( ($contig_name == $old_contig_name) || ($old_contig_name == "start") ) {
		# count query hits
		$counter++;
        $old_contig_name = $contig_name;
		push(@blastresults, \@blastline);
		# count new chromosome hits 
		unless ($blastline[1] == $db_chr){
	    $db_chr_counter++;
	    $db_chr = $blastline[1];
		}
    }


    ### new query
    elsif ($contig_name ne $old_contig_name)  {

	#  GET INFO FOR LAST QUERY ###
        (my $line_array_ref, $length_disparity_counter)  = &hash_data($old_contig_name, \@blastresults, $counter, $db_chr_counter, $length_disparity_counter);
        $blasthash{$old_contig_name} = $line_array_ref;

	print OUT "@{$blasthash{$old_contig_name}}\n"; 

	## RESET VARIABLES  ##
	$contig_name = $blastline[0];
	$contig_match_length;
	$counter = 1; # number of hits
	$db_chr_counter = 1;
	
	## empty blastresults
	@blastresults = ();
	
	## START NEW QUERY CHUNK ##
	$old_contig_name = $contig_name;
	push(@blastresults, \@blastline);
	$db_chr = $blastline[1];
	
    }
    if ( eof($blastfile) ){
        ##  GET INFO FOR LAST QUERY ###
        (my $line_array_ref, $length_disparity_counter)  = &hash_data($contig_name, \@blastresults, $counter, $db_chr_counter, $length_disparity_counter);
		$blasthash{$contig_name} = $line_array_ref;
		print OUT "@{$blasthash{$old_contig_name}}\n";
	next;
    }
	
}

# store blasthash to file for later use by vcf annotation scripts
print "Storing blast information hash to file for use by vcf annotation scripts...\n";
store \%blasthash, 'BlastHash.dat';




#########################################################################################################################
#####         SUBROUTINES
#########################################################################################################################


# sub transpose from http://www.perlmonks.org/?node_id=15209
sub transpose {
    map {
        my $j = $_;
    [ map $_[$_][$j], 0..$#_ ]
    } 0..$#{$_[0]};
}


sub hash_data {
    # takes in blasthash to append to, current block of blast results to work on
    # the counter of how many hits there are, the counter of how many different 
    # db sequences were hit
    # sorts the blast ouput block by evalue, and takes the lowest evalue result,
    # but skipping matches of less than %80 length of the query
    # outputs updated hash, with seleted blast hit for the query and counter info

    (my $contig_name,  my $blastresults, my $counter, my $db_chr_counter, my $length_disparity_counter) = @_; # my @temp, my @sorted) = @_;
    my $line_array_ref;
    # sort by evalue
    my @temp =
        transpose
	sort { $a->[11] <=> $b->[11] } @blastresults;
    my @sorted =
	transpose @temp;
    # loop through sorted query sequence blast results
    foreach my $element (@sorted){
        # skip if the length of the match is less than 80% total query length
	if ($element->[3] < $element->[4]*0.8){
	    $length_disparity_counter++;
	    next;
	}
       # get the next best hit
	else {
	    # add on number of hits for this contig (>10 indicates it may be a repeated region)
	    push (@$element, $counter);
	    # record the number of chromosomes hit by this query
	    push (@$element, $db_chr_counter); 
	    $line_array_ref = $element;
	    last;
	}
    }
    return ($line_array_ref, $length_disparity_counter);
}







