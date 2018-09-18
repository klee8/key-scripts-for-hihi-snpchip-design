# Kate Lee Oct 2016
# overlaps_filtering.pl
# script to look for overlaps in the hihi snps that mapped to ZF and filter them
#
# usage: perl overlaps_filtering.pl <filelist>
###############################################################################3

#!usr/bin/perl -w
use strict;
use Storable;
use List::Util qw( reduce );


my $filelist = $ARGV[0];
my $merged = $ARGV[1] || 'logfile';
my %list;
my $vcf_file;
my %pos_info;


# inherited headers
my @wgsheaders = qw|#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8 sample9 sample10 none none none none none none none none none none none none none none none none none none none none none Tiri.00 Tiri.01 Tiri.11 LBI.00 LBI.01 LBI.11 REF_Tiri REF_LBI ALT_Tiri ALT_LBI REF ALT SNP_type inZF ZFchrom ZFpos lenLHS lenRHS qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore num_query_hits num_chrs_hit ENSgene ENStran chrname gene_start gene_end strand ass_gene_name description ass_with_snp current_ass ass_list gapsize min_flanking_contig max_flanking_contig mapqual RPB MQ0F strand_bias read_depth distL distR |;

my @radheaders = qw |#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 24792_1 59629_1 72002_1 72003_1 72005_1 72006_1 72007_1 72008_1 72009_1 72010_1 72011_1 72012_1 72015_1 72019_1 72020_1 72021_1 75649_1 75667_1 76055_1 79372_1 83382_1 C72016_1 C72017_1 C72018_1 C75621_1 C83383_1 P2753_1 P2802_1 P2825_1 P2883_1 P5007_1 Tiri.00 Tiri.01 Tiri.11 LBI.00 LBI.01 LBI.11 REF_Tiri REF_LBI ALT_Tiri ALT_LBI REF ALT SNP_type inZF ZFchrom ZFpos lenLHS lenRHS qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore num_query_hits num_chrs_hit ENSgene ENStran chrname gene_start gene_end strand ass_gene_name description ass_with_snp current_ass ass_list gapsize min_flanking_contig max_flanking_contig mapqual RPB MQ0F strand_bias read_depth distL distR |;


my $lenwgs = @wgsheaders;
my $lenrad = @radheaders;
print "checking length wgs: $lenwgs   and length rad: $lenrad\n";


# get annotated vcf assembly file location info
open (LIST, "<$filelist ") || die "ERROR: cannot open $filelist: $!";
while(<LIST>){
    chomp;
    my @values = split("\t", $_); 
    $list{$values[0]} = $values[1];
}
close LIST;


# open file for rejected snps
open (REJECT, ">rejected_mapped_snps.txt") || die "ERROR: couldn't open unmapped_snps.txt: $!";

# open outfile for merged data
open (my $merged, ">$merged") || die "ERROR: couldn't open merged outfile: $!";
print "finding overlaps in filetest list...\n";


my $wgs = join("\t", @wgsheaders);
my $rad = join("\t", @radheaders);
print $merged "$wgs\n$rad\n";




# loop through filelist, make lines the same length and parse info to hash  {ZFchrom}{position}{assembly} = (lineinfo)
for my $key (sort keys %list){
    print  "parsing file $key...\n";
    open (IN, "<$key") || die "ERROR: cannot open $key: $!";
    my $linecount = 0;
    # loop through each file
    while(<IN>){
	chomp;
	if ($_=~ /^#/) {next;}
	else {
	    $linecount++;
	    my @row_values = split("\t", $_);
	    my $current_ass = $list{$key};

            #### if not a RAD line, add extra columns to make them the same length
	    my @emptycolumns = qw |none none none none none none none none none none none none none none none none none none none none none|;
	    unless ($current_ass eq 'RAD') {
		my @line1 = @row_values[0..18];
		my @line2 = @row_values[19..59];
		@row_values = ();
		push(@row_values,  (@line1, @emptycolumns, @line2));
	    }


	        # skip lines that don't have blast values ( shouldn't have any unblasted things in here )
	    if  ( $row_values[58] eq "NA" )  { 
		#print "capturing unblasted data....\n";
		print REJECT "noblastvalue:\t@row_values\n"; 
		next; 
	    }
	    else {
		# positional_hash{ZF_chr}{ZF_position}{Assembly} = ref_pos.ann.vcf_line array 
		$pos_info{$row_values[54]}{$row_values[55]}{$list{$key}} = \@row_values;
		
	    }
	    #print "$row_values[54]   $row_values[55]   $list{$key}\n";                # TEST collecting ZF chr and position ;)
	    
	}
    }
    close IN;
}



# loop through postional hash to apply positional info
my $two_pos_back = 'starting';
my $one_pos_back = 'starting';
my $current_pos = 'starting';
my $two_ass_back = 'starting';
my $one_ass_back = 'starting';
my $current_ass = 'starting';
my $gapsize_two_back = 0;
my $gapsize_one_back = 0;
my $gapsize = 0;
my @current_line;
my %info;

# Loop through chromosomes
for my $chr (sort keys %pos_info){
    print "sorting through ZF chromosome $chr...\n";
    my $max_val_key = reduce { ${$pos_info{$chr}}{$a} > ${$pos_info{$chr}}{$b} ? $a : $b } keys %{$pos_info{$chr}};
    print "max snp position in chromosome $chr is $max_val_key\n";
    # loop through ZF positions, calculate and store info for one_pos_back
    for my $pos (sort {$a<=>$b} keys %{$pos_info{$chr}}){                       ## Note: must sort positional information numerically to get distance to closest snp later 
	$two_pos_back = $one_pos_back;
	$one_pos_back = $current_pos;
	$current_pos = $pos;
	$two_ass_back = $one_ass_back;
	$one_ass_back = $current_ass;
	$gapsize_two_back = $gapsize_one_back;
	$gapsize_one_back = $gapsize;
	@current_line = ();


	################      Find 'best version' of each snp across assemblies      #################

	my @assemblies=();
	my $max_qual = 7;
	my $min_flanking_contig = 0;
	my $max_flanking_contig = 0;
	# get best assembly based on qual and contig length
	# collect #assemblies and their names
	for my $ass ( sort keys %{$pos_info{$chr}{$pos}} ){
	    push(@assemblies, $ass);

	    my $lenLHS = ${$pos_info{$chr}{$pos}{$ass}}[56];
	    my $lenRHS = ${$pos_info{$chr}{$pos}{$ass}}[57];
	    $min_flanking_contig = ($lenLHS, $lenRHS)[$lenLHS > $lenRHS];  # returns the smaller number
	    $max_flanking_contig = ($lenLHS, $lenRHS)[$lenLHS < $lenRHS];  # returns larger number
            # chuck anything that doesn't have a minimum flanking contig of 36                     
            ##### if the numbers are too small can use contig length = 90 instead!!!                                          # NOTE    THIS CAN BE CHANGED 
	    # get assembly with required flanking seq length and highest quality score 
	    if ( ( ( $min_flanking_contig > 9) && ($max_flanking_contig > 34)  && (${$pos_info{$chr}{$pos}{$ass}}[5] > $max_qual ) ) || ( $max_qual == 7 ) ) {
		$max_qual = ${$pos_info{$chr}{$pos}{$ass}}[5]; 
		$current_ass = $ass;
	    }
	}




        ###############       Set current line and add info for current line                    ##################

        @current_line = @{$pos_info{$chr}{$pos}{$current_ass}};


        # assembly info
	my $num_ass_with_snp = @assemblies;
	my $ass_list = join(":", @assemblies);
	

        # check gapsize for current pos sqrt( ((qend - qstart) - (send - sstart))**2 )
        $gapsize = sqrt( ( sqrt((${$pos_info{$chr}{$pos}{$current_ass}}[66] - ${$pos_info{$chr}{$pos}{$current_ass}}[65])**2) - sqrt((${$pos_info{$chr}{$pos}{$current_ass}}[67] - ${$pos_info{$chr}{$pos}{$current_ass}}[65])**2) )**2);


	# set snp key value
	my $mergedSNP = "$current_line[3]$current_line[4]";
	if ( ($mergedSNP =~ /A/) &&  ($mergedSNP =~ /T/) ) { $current_line[52] = 'W'; }
        if ( ($mergedSNP =~ /C/) &&  ($mergedSNP =~ /G/) ) { $current_line[52] = 'S'; }
        if ( ($mergedSNP =~ /A/) &&  ($mergedSNP =~ /G/) ) { $current_line[52] = 'R'; }
        if ( ($mergedSNP =~ /C/) &&  ($mergedSNP =~ /T/) ) { $current_line[52] = 'Y'; }
        if ( ($mergedSNP =~ /G/) &&  ($mergedSNP =~ /T/) ) { $current_line[52] = 'K'; }
        if ( ($mergedSNP =~ /A/) &&  ($mergedSNP =~ /C/) ) { $current_line[52] = 'M'; }
     

	# Grab info column
	%info;
	foreach my $value (split(";", @{$pos_info{$chr}{$pos}{$current_ass}}[7])){   
	    my @temp = split("=", $value);
	    my $key = $temp[0];
	    $info{$key} = $temp[1];

	}

        my $mapqual = $info{'MQ'} || 'NA';
        my $RPB = $info{'RPB'} || 'NA';
        my $MQ0F = $info{'MQ0F'} || 'NA';
        my @DP4 = split(",", $info{'DP4'});
        my $strand_bias = ( ($DP4[1] + $DP4[3]) > 0) && ( ($DP4[0] + $DP4[2]) / ($DP4[1] + $DP4[3]) ) || 'NA';
        my $read_depth = $info{'DP'} || 'NA';



        push( @current_line, ($num_ass_with_snp, $current_ass, $ass_list, $gapsize, $min_flanking_contig, $max_flanking_contig, $mapqual, $RPB, $MQ0F, $strand_bias, $read_depth ) );

	@{$pos_info{$chr}{$pos}{$current_ass}} = @current_line;



	
	################################################################################################
	################     Calculate min Distances to next snp for previous line     #################
	################################################################################################

	# Note distances can be negative if the gaps in one or both blast results exceed the distance between assumed snp positions 


        # max gapsize for one pos back to prev
        my $distL = ($one_pos_back - $two_pos_back);  # - ($gapsize_two_back + $gapsize_one_back);

        # max gapsize for one pos back to next
        my $distR = ($pos - $one_pos_back);   # - ($gapsize_one_back + $gapsize);
	

        ###############       Add distance to next SNP (L and R) to the previous info line and filtering                    ##################
	my @prev_line;
	my $printline;
	unless ($one_ass_back eq 'starting') {
	    @prev_line = @{$pos_info{$chr}{$one_pos_back}{$one_ass_back}};
	    #print "current:$pos:$current_ass prev:$one_pos_back:$one_ass_back @prev_line\n";
	    push (@prev_line, ($distL, $distR));
	    #print "put distL:@prev_line[92] distR:@prev_line[93]\n";
	    my $subline = join("\t", @prev_line);
	    my $line = &filter($subline, \*REJECT);
	    print $merged $line;
	}

        if ($pos == $max_val_key) {
            # if last position in chromosome calculate min Distances to next snp for CURRENT line  
            # max distance to previous snp
            my $distL = ($pos - $one_pos_back);   # ($gapsize_one_back + $gapsize for total potential gap length);
            my $distR = 100;
	    push (@current_line, ($distL, $distR));
	    my $subline = join("\t", @current_line);
            my $line = &filter($subline, \*REJECT);
            print $merged $line;

	}


    }
}




#############################################################################################
###############        FILTERING Previous line (with dist measures)      ####################
#############################################################################################




sub filter{
    my $subline = @_[0];
    my $REJECT = $_[1];
    my $reject = 0;
    my $printline = "";
    my @prev_line = split("\t", $subline);
    
    # blast number of chromosomes hit = 1   (snp region maps to one chromosome)
    if ($prev_line[72] > 1) { $reject++; if ($reject == 1) {print $REJECT "rejecting_chromhits:\t$subline\n";}}
    
    # blast number of query hits < 10  (snp not in repetitive region)
    if ($prev_line[71] > 10) { $reject++; if ($reject == 1) {print $REJECT "rejecting_queryhits:\t$subline\n";}}

    # get rid of indels
    if ($prev_line[7] =~ /^INDEL/) { $reject++; if ($reject == 1) {print $REJECT "rejecting_INDELhits:\t$subline\n";}}
    
    # biallelic only
    if (length($prev_line[4]) > 1) { $reject++; if ($reject == 1) {print $REJECT "rejecting_nonbialellic:\t$subline\n";}}

    # distance to nearest snp - minimum 40bp (do this after you have the fasta sequences - check nearby snps aren't the same snp)
    #if ( ($prev_line[92] < 40) || ($prev_line[93] < 40) ) { $reject++; if ($reject == 1) {print $REJECT "rejecting_disttonextsnp:\t$subline\n"; }}
    
    # poly in Tiri (Tiri_REF > 1 or Tiri_ALT > 1)
    if ( ($prev_line[46] < 1) || ($prev_line[48] < 1) ) { $reject++; if ($reject == 1) {print $REJECT "rejecting_notPolyInTiri:\t$subline\n"; }}

    # Read Depth > 19
    #if (length($prev_line[91]) < 19) { $reject++; if ($reject == 1) {print $REJECT "rejecting_readdepth:\t$subline\n";}}
						      
    # min_flanking_contig < 6, max_flanking_contig < 35
    if ( ($prev_line[85] < 6) || ($prev_line[86] < 35) ) { $reject++; if ($reject == 1) {print $REJECT "rejecting_flankingContigLength:\t$subline\n"; }}

    # snp type
    unless ($prev_line[52] =~ /[RYKM]/) { $reject++; if ($reject == 1) {print $REJECT "rejecting_snptype:\t$subline\n"; }}

#    if ($reject > 0) { $printline = join("\t", @prev_line); print "REJECTED:\t$printline";}

    if ($reject == 0) {
	$printline = join("\t", @prev_line);
	$printline = $printline."\n";
    }


    return $printline;
}


