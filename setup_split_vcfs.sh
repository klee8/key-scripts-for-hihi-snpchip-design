# create awk script to divide each file by chromosome for WGS and RAD data

#!/bin/bash


# split wgs data (chr in column 34 in vcf file)

for i in {1..10} 3in1 10in1; do mkdir assembly_$i; 
echo '''grep -v "#" ../../8_add_ZF_locations/POS_TEMP/assembly.TEMP.ann.pos.vcf | awk REPLACE{print >> $34."_subset.vcf" ; close($34."_subset.vcf")}REPLACE ''' > assembly_$i/split_vcf_$i.sbatch; sed -i "s/TEMP/$i/g" assembly_$i/split_vcf_$i.sbatch; sed -i "s/REPLACE/'/g" assembly_$i/split_vcf_$i.sbatch; done



# split RAD data (chr in column 55 in vcf file)

mkdir assembly_RAD
echo '''grep -v "#" ../../8_add_ZF_locations/POS_RAD/assembly.RAD.ann.pos.vcf | awk REPLACE{print >> $55."_subset.vcf" ; close($55."_subset.vcf")}REPLACE ''' > assembly_RAD/split_vcf_RAD.sbatch; sed -i "s/REPLACE/'/g" assembly_RAD/split_vcf_RAD.sbatch




