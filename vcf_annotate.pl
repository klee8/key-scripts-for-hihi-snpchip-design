# Kate Lee Feb 2016
# vcf_annotate.pl
###################################################################################################

#!usr/bin/perl -w
use strict;

=head1 NAME

vcf_annotate.pl

=head1 USAGE

perl vcf_annotate.pl input_vcf_file output_vcf_annotated config_file

=head1 DESCRIPTION

reads in vcf file
outputs annotated vcf to help rate snps for snp-chip


=head1 ARGUMENTS

=over 4

=item • input_vcf

An input VCF file containing snp calls.

=item • output_vcf

Name of an output vcf file which will contain additional columns with new annotation for snpchip design.

=item • config_file

A tab file listing sample names in the first colum and population in the second.

=back

=head1 AUTHOR

Kate Lee

=cut

#Verify that there are three arguments.
if(@ARGV!=3)
{die("Perl script failed: must have 3 arguments, exiting;");}

#Set variables from arguments.
my ($infile_vcf, $outfile_vcf, $config_file)  = @ARGV;


#####    Open files    #####
open (my $vcf_fh, "<$infile_vcf") || die "couldn't open vcf file: $!";
open (my $out_vcf, ">$outfile_vcf") || die "couldn't open outfile: $!";
open (my $config, "<$config_file") || die "couldn't open config file: $!";


#####     Read in Config info to hash     #####
my %sample_pop;
my $pop_ref = \%sample_pop;
print "Reading in config file...\n";
while(<$config>){
    chomp;
    (my $sample, my $pop) =  split("\t", $_, 2);
    $sample =~ s/\s//g;
    $pop =~ s/\s//g;
    $sample_pop{$sample}=$pop;
}


#####     Read in VCF file     #####
my @col_headers;
my $col_ref = \@col_headers;
my %col_header_index;
my $index_ref = \%col_header_index;
my @row_values;
my $row_ref = \@row_values;
my @annotation_headers = ('Tiri.00', 'Tiri.01', 'Tiri.11', 'LBI.00', 'LBI.01', 'LBI.11', 'REF_Tiri', 'REF_LBI', 'ALT_Tiri', 'ALT_LBI', 'REF', 'ALT', 'SNP_type');

my $read_depth = 20;
my $genotype_quality = 20;


while(<$vcf_fh>){
    chomp;

# if header line, print straight to annotation output
    if ($_=~ m/^##/){
	print $out_vcf "$_\n";
	next;
    }
# check vcf column headers and create header index
    elsif ($_=~ m/^#CHROM/){
	unless ($_=~ m/#CHROM\tPOS\tID\tREF\tALT\t/){
	    die "expected input vcf to have column headers #CHROM\tPOS\tID\tREF\tALT\tetc...."; 
	}
	# add in header to show this script has run
	print $out_vcf "## vcf_annotate.pl @ARGV\n";
	# add new header line with extra column annotations
	my $extra_headers = join("\t", @annotation_headers);
	print $out_vcf "$_\t$extra_headers\n";
	@col_headers = split("\t|\n", $_);
	my $num_of_cols = @col_headers;
	my $index = 0;
	foreach my $header (@col_headers){
	    $header=~ s/\s//g;
	    $col_header_index{$header} = $index;
	    $index++;
	}
	print "index number for FORMAT header is $col_header_index{'FORMAT'}\n";
	my $num_of_samples = $num_of_cols - $col_header_index{'FORMAT'}-1;
	print "Finished loading headers...\n";
	print "Number of samples calculated = $num_of_samples\n";
	next;
    }
# read in tab delim into array
    else {
	@row_values = split("\t|\n", $_);
	if (&id_monomorphic($row_ref, $index_ref) eq "TRUE") {
	    next;
	}
        my $count_ref = &count_genotypes($row_ref, $index_ref, $pop_ref, $read_depth, $genotype_quality);
	my %counts = %$count_ref;
	my $new_values;
        for my $value (@annotation_headers){
	    $new_values = "$new_values\t$counts{$value}";
	}
	print $out_vcf $_.$new_values."\n";
    }
}




######################################################################################
#                       SUBROUTINES
######################################################################################



# identify monomorphic lines
sub id_monomorphic($row_ref, $index_ref){
    # return true if snp is monomorphic
    my @row_values = @$row_ref;
    my %column_index = %$index_ref;
    my @info_values = split(":", $row_values[$column_index{'INFO'}]);
    unless ($info_values[1] eq "AN=1"){
	return "FALSE";
    }
    else {
	return "TRUE";
    }
}


# counting genotypes
sub count_genotypes($row_ref, $index_ref, $pop_ref, $read_depth, $genotype_quality){
    my @row_values = @$row_ref;
    my %column_index = %$index_ref;
    my %sample_pop = %$pop_ref;
    my %counts;
    my @annotation_headers = qw|Tiri.00 Tiri.01 Tiri.11 LBI.00 LBI.01 LBI.11 REF_Tiri REF_LBI ALT_Tiri ALT_LBI REF ALT SNP_type|;
    foreach my $category (@annotation_headers){
	$counts{$category} = 0;
    }
    print "########### counting genotypes......\n";
    my $i = $column_index{"FORMAT"} + 1;
    my %rhash = reverse %column_index;
   

    ### LOOP THROUGH GENOTYPE CALLS FOR EACH SAMPLE 
    for ($i > $column_index{"FORMAT"} ; $i < scalar @row_values; $i++){
	my $col_name = $rhash{$i};
	my @genotype_info = split(":", $row_values[$i]);
	my $pop = $sample_pop{$col_name};


	### FILTER OUT POOR GENOTYPE CALLS
	if ($genotype_info[0] =~ m/\./){
	    next;
	}

        ### COUNTING STARTS HERE
        else {   
            # GENOTYPE INFO
	    $genotype_info[0] =~ s/\///g;
	    if ($genotype_info[0] == '10') {$genotype_info[0] = '01';}
	    print $col_name.":".$pop.":".$genotype_info[0]."\n";

	    ### COUNT GENOTYPES PER POPULATION
	    $counts{"$pop.$genotype_info[0]"}++;
	    
	    ### COUNT REF and ALT genes
	    if (  ($genotype_info[0] eq '00') && ($pop eq 'Tiri')  ) {$counts{'REF_Tiri'} = $counts{'REF_Tiri'} + 2;}
	    if (  ($genotype_info[0] eq '00') && ($pop eq 'LBI')   ) {$counts{'REF_LBI'} = $counts{'REF_LBI'} + 2;}
	    if (  ($genotype_info[0] eq '01') && ($pop eq 'Tiri')  ) {$counts{'REF_Tiri'}++; $counts{'ALT_Tiri'}++;}
	    if (  ($genotype_info[0] eq '01') && ($pop eq 'LBI')   ) {$counts{'REF_LBI'}++; $counts{'ALT_LBI'}++;}
	    if (  ($genotype_info[0] eq '11') && ($pop eq 'Tiri')  ) {$counts{'ALT_Tiri'} = $counts{'ALT_Tiri'} + 2;}
	    if (  ($genotype_info[0] eq '11') && ($pop eq 'LBI')   ) {$counts{'ALT_LBI'} = $counts{'ALT_LBI'} + 2;}

	}


    }

    ### note for REF and ALT whether it is present in LBI/Tiri/both ... $REF, $ALT
    if ( ($counts{'REF_Tiri'} > 0) && ($counts{'REF_LBI'} == 0) ) {$counts{'REF'} = "TIRI";}
    if ( ($counts{'REF_Tiri'} == 0) && ($counts{'REF_LBI'} > 0) ) {$counts{'REF'} = "LBI";}
    if ( ($counts{'REF_Tiri'} > 0) && ($counts{'REF_LBI'} > 0) ) {$counts{'REF'} = "BOTH";}
    if ( ($counts{'ALT_Tiri'} > 0) && ($counts{'ALT_LBI'} == 0) ) {$counts{'ALT'} = "TIRI";}
    if ( ($counts{'ALT_Tiri'} == 0) && ($counts{'ALT_LBI'} > 0) ) {$counts{'ALT'} = "LBI";}
    if ( ($counts{'ALT_Tiri'} > 0) && ($counts{'ALT_LBI'} > 0) ) {$counts{'ALT'} = "BOTH";}

    return \%counts;
}















