# Kate Lee September 2016
# vcf.pos.ann.pl
#
# takes in blasthash.dat file (from blastparser.pl output)
# takes in vcf file output frm vcf_annotate.pl (with pop and allele summaries)
# iterates through vcf file and calculates
#      ZF position of snp (approximate if there are gaps)
#      lhs flanking length
#      rhs flanking length
# outputs vcf with added annotation
#      ZF postion of snp
#      lhs flanking length
#      rhs flanking length
#      blast line info (fmt 6 plus qlen, num_hits for that query, num_chrs hit in ZF)
#      
#      vcf positional information hash (zf position, no contigs hit by hihi, no chrs hit)
# usage: vcf.pos.ann.pl <ann_vcf_file> <blast.dat> <vcf_pos_output_file>
##########################################################################

#!usr/bin/perl -w
use strict;
use Storable;


# check input is correct

if(@ARGV!=3){
    print " usage: vcf.pos.ann.pl <ann_vcf_file> <blast.dat> <vcf_pos_output_file> ";
    exit;
}

(my $annVcfFile, my $blastDatFile, my $vcfOUT) = @ARGV;


# get blasthash data from .dat file
my $i, my $href = retrieve("$blastDatFile", {binmode=>':raw'});
my %blasthash = %$href; 

# get ZFgenehash data from ZFgeneHash.dat file
my $i, my $href = retrieve("../ZF_genes/ZFgeneHash.dat", {binmode=>':raw'});
my %ZFgeneHash = %$href; 

# create hash for snp ZF positions
my %poshash;


#open input and output vcf files
print "opening input and output files...\n";
open (IN, "<$annVcfFile") || die "ERROR: couldn't open vcf file $annVcfFile: $!";
open (my $out_vcf, ">$vcfOUT") || die "ERROR: couldn't open outfile $vcfOUT: $!";

print "adding header annotation...\n";
my @annotation_headers = qw|inZF ZFchrom ZFpos lenLHS lenRHS qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore num_query_hits num_chrs_hit|;
my @zfheaders = qw|Ensembl_Gene_ID Ensembl_Transcript_ID Chromosome_Name Gene_Start Gene_End Strand Associated_Gene_Name Description|;
my @col_headers;
my $col_ref = \@col_headers;
my %col_header_index;
my $index_ref = \%col_header_index;
my @row_values;
my $row_ref = \@row_values;

my $counter = 0; 

# loop through vcf
print "looping through vcf file...\n";
while(<IN>){
    $counter++;
    chomp;
    # if header line, print straight to annotation output
    if ($_=~ m/^##/){
        print $out_vcf "$_\n";
        next;
    }

    # if vcf column headers, create new headers index to print out
    elsif ($_=~ m/^#CHROM/){
        unless ($_=~ m/#CHROM\tPOS\tID\tREF\tALT\t/){
            die "expected input vcf to have column headers #CHROM\tPOS\tID\tREF\tALT\tetc....";
	}
	# print out that this script has been used
	print $out_vcf "##vcf.positional.annotations.pl @ARGV\n";
	
	# print out new header line with addtional annotations
	my $extra_headers = join("\t", @annotation_headers);
	my $zfann_headers = join("\t", @zfheaders);
        print $out_vcf "$_\t$extra_headers\t$zfann_headers\n";
        @col_headers = split("\t|\n", $_);
        my $num_of_cols = @col_headers;
        my $index = 0;
        foreach my $header (@col_headers){
            $header=~ s/\s//g;
            $col_header_index{$header} = $index;
            $index++;
        }
        print "Finished loading headers...\n";
        next;
    }

    #else get positional info
    else {
	@row_values = split("\t", $_);
	#print "$row_values[0]\n";
	my @new_values = qw|inZF ZFchrom ZFpos lenLHS lenRHS|;
		
	# if there is no blast info for this line, fill extra columns with NA
	unless (exists($blasthash{$row_values[0]})){
	    print "$row_values[0] contig has no blast hit....\n";
	    foreach my $header (@new_values){
		$poshash{$header} = 'NA';
	    }
	    my $emptyblast = "NA\t"x23;
	    my $vcf_values = join("\t", @row_values);
	    print $out_vcf "$vcf_values\t$poshash{'inZF'}\t$poshash{'ZFchrom'}\t$poshash{'ZFpos'}\t$poshash{'lenLHS'}\t$poshash{'lenRHS'}\t$emptyblast\n";
	}

	# if there is a blast hit get pos info
	else {
	    
	    print "$counter\t$row_values[0]\tblasthash_array:$blasthash{$row_values[0]}\tfirst_element:${$blasthash{$row_values[0]}}[0]\n";
	    my $ref = $blasthash{$row_values[0]};
	    my @blastinfo = @$ref;
	    $poshash{'inZF'} = 1;
	    $poshash{'ZFchrom'} = $blastinfo[1];
	    
	    # if the match has both fw sequences
	    if ( ( ($blastinfo[10] - $blastinfo[9]) > 0) && ( ($blastinfo[8] - $blastinfo[7]) > 0) ) {
		$poshash{'ZFpos'} = $blastinfo[9] + $row_values[1] - $blastinfo[7];
		print "fw sequences\n";
	    }
	    # if the blast match has ZF sequence fw and hihi rev
	    elsif ( ( ($blastinfo[10] - $blastinfo[9]) > 0) && ( ($blastinfo[8] - $blastinfo[7]) < 0) ) {
		$poshash{'ZFpos'} = $blastinfo[10] + $row_values[1] - $blastinfo[7];
		print "hihi reverse sequence\n";
	    }
	    # if the blast match has ZF sequence rev and hihi fw
	    elsif ( ( ($blastinfo[10] - $blastinfo[9]) < 0) && ( ($blastinfo[8] - $blastinfo[7]) > 0) ) {
		$poshash{'ZFpos'} = $blastinfo[10] + $row_values[1] - $blastinfo[7];
		print "ZF reverse sequence\n";
	    }
	    # flag if something else is happening
	    else {
		print "someting else happening here? hihi contig: $ref  query start:$blastinfo[7] end:$blastinfo[8]  ref start:$blastinfo[10]  end:$blastinfo[9]\n";
	    }
	    
	    # GENE ANNOTATION snps???
	    my $nearest_gene = 0;
	    my $gene_start;
	    my $gene_ann = "NA\t"x8;
	    my @ann_array; 
	    unless ($poshash{'inZF'} == 0){
		my $position = $poshash{'ZFpos'};
		for $gene_start (sort keys %{$ZFgeneHash{"$poshash{'ZFchrom'}"}}) {
		    if ($gene_start < $poshash{'ZFpos'} ) {
			$nearest_gene = $gene_start; 	
		    }
		    else {last;}
		}
		print "nearest gene = $nearest_gene\n";
		if ( ( $nearest_gene < $poshash{'ZFpos'} ) && ( $poshash{'ZFpos'} < $ZFgeneHash{"$poshash{'ZFchrom'}"}{$nearest_gene}{'Gene_End'} ) ) {
		    foreach my $header (@zfheaders){
			push (@ann_array, $ZFgeneHash{"$poshash{'ZFchrom'}"}{$nearest_gene}{$header});
		    }
		    $gene_ann = join("\t", @ann_array);
		    print "GENE ANNOTATIONS!!!!!!!!!!! $gene_ann....\n";
		}
	    }
	    #length of LHS of contig = position of snp - 1
	    $poshash{'lenLHS'} = $row_values[1] - 1;
	    #length of RHS of contig = length of query - position of snp
	    $poshash{'lenRHS'} = $blastinfo[3] - $row_values[1];
	    my $vcf_values = join("\t", @row_values);
	    my $blast_values = join("\t", @blastinfo);
	    print $out_vcf "$vcf_values\t$poshash{'inZF'}\t$poshash{'ZFchrom'}\t$poshash{'ZFpos'}\t$poshash{'lenLHS'}\t$poshash{'lenRHS'}\t$blast_values\t$gene_ann\n";
	    
	}
    }
    
}    
    
    
