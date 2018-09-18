# Key set of scripts used in the design of the hihi snpchip. 

### Kate D Lee Sept 2018

MIT licenced


### divide_into_equal_regions.pl
Used to divide bam files into regions of equal length to parallelise GATK re-alignment

### vcf_annotate.pl
Used to add information to the end of the vcf files

### blastparse.pl
Used to interrogate BLAST information for SNP positions and gene information

### vcf.pos.ann.pl
Used to integrate BLAST information into the vcf files

### filter_read_depth.pl
Used to do an initial pass of read depth filtering on vcfs

### setup_split_vcfs.sh
Used to set up an awk script to split vcf files by zebra finch chromosome (added during positional annotation)

### merge_filter.pl
Used to merge the SNP information from different assemblies using the zebra finch position information

### get_snps_fasta.pl
Used to get the fasta sequences flanking filtered SNPs


