#!/bin/bash

#SBATCH -o ~/slurm_jobs/slurm_%j.log
#SBATCH -e ~/slurm_jobs/slurm_%j.err
#SBATCH -J rRNA_remove
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH -c 4
#SBATCH -t 48:00:00
#SBATCH --mem=100G

echo "SLURM Job ID        : ${SLURM_JOB_ID}";
echo "SLURM Job name      : ${SLURM_JOB_NAME}";
echo "SLURM Node list     : ${SLURM_JOB_NODELIST}";
echo "SLURM partition     : ${SLURM_JOB_PARTITION}";

source $HOME/.bashrc
echo "activate conda";
conda activate RNA-seq;
echo "activated conda";

id=$1
cwd_path=$2
chrom_sizes=$3
cores=$4
gtf_file=$5
rRNA=$6

echo "working on $id"
mkdir -p $cwd_path
cd $cwd_path
temp_d=$(mktemp -d --tmpdir=$HOME/tmp)
echo "temporary files go to $temp_d"

echo "=========================================================="
echo "Sorting rRNA bed file... [ $(date) ]"
###Sort rRNA bed file to match order
bedtools sort -i $rRNA -faidx $chrom_sizes > $temp_d/rRNA_coords.bed

echo "=========================================================="
echo "Copying bam files... [ $(date) ]"

#5 Prime
cp $cwd_path/$id.FivePrime_BothStrands.Primary.Multimapping.Read2.sorted.bam $temp_d/Aligned.Primary.Multimapping.Read2.sorted.bam
cp $cwd_path/$id.FivePrime_BothStrands.Primary.Multimapping.Read2.sorted.bam.bai $temp_d/Aligned.Primary.Multimapping.Read2.sorted.bam.bai

cp $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.Read2.sorted.bam $temp_d/Aligned.Primary.Unique.Read2.sorted.bam
cp $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.Read2.sorted.bam.bai $temp_d/Aligned.Primary.Unique.Read2.sorted.bam.bai

#else
cp $cwd_path/$id.Else.unique.sorted.bam $temp_d/else.primary.unique.sorted.bam
cp $cwd_path/$id.Else.unique.sorted.bam.bai $temp_d/else.primary.unique.sorted.bam.bai

cp $cwd_path/$id.Else.all.multimapping.sorted.bam $temp_d/else.all.multimapping.sorted.bam
cp $cwd_path/$id.Else.all.multimapping.sorted.bam.bai $temp_d/else.all.multimapping.sorted.bam.bai

#all
cp $cwd_path/$id.All.unique.sorted.bam $temp_d/all.primary.unique.sorted.bam
cp $cwd_path/$id.All.unique.sorted.bam.bai $temp_d/all.primary.unique.sorted.bam.bai

cp $cwd_path/$id.All.all.multimapping.sorted.bam $temp_d/all.all.multimapping.sorted.bam
cp $cwd_path/$id.All.all.multimapping.sorted.bam.bai $temp_d/all.all.multimapping.sorted.bam.bai

#5 Prime - both strands
cp $cwd_path/$id.FivePrime_BothStrands.Primary.Multimapping.BothStrands.sorted.bam $temp_d/Aligned.Primary.Multimapping.BothReads.sorted.bam
cp $cwd_path/$id.FivePrime_BothStrands.Primary.Multimapping.BothStrands.sorted.bam.bai $temp_d/Aligned.Primary.Multimapping.BothReads.sorted.bam.bai

cp $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.BothStrands.sorted.bam $temp_d/Aligned.Primary.Unique.BothReads.sorted.bam
cp $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.BothStrands.sorted.bam.bai $temp_d/Aligned.Primary.Unique.BothReads.sorted.bam.bai


echo "=========================================================="
echo "Removing rRNA mapping reads from 5' data - both strands... [ $(date) ]"
#5' multimapping
bedtools intersect -a $temp_d/Aligned.Primary.Multimapping.BothReads.sorted.bam -b $temp_d/rRNA_coords.bed -sorted -g $chrom_sizes > $temp_d/rRNAreads.bam
samtools view $temp_d/rRNAreads.bam | awk '{print $1}' > $temp_d/filter.txt
picard FilterSamReads --QUIET true -I $temp_d/Aligned.Primary.Multimapping.BothReads.sorted.bam -O $temp_d/Aligned.Primary.Multimapping.BothReads.sorted.norRNA.bam --READ_LIST_FILE $temp_d/filter.txt --FILTER excludeReadList
samtools index $temp_d/Aligned.Primary.Multimapping.BothReads.sorted.norRNA.bam


#5' unique
bedtools intersect -a $temp_d/Aligned.Primary.Unique.BothReads.sorted.bam -b $temp_d/rRNA_coords.bed -sorted -g $chrom_sizes > $temp_d/rRNAreads.bam
samtools view $temp_d/rRNAreads.bam | awk '{print $1}' > $temp_d/filter.txt
picard FilterSamReads --QUIET -I $temp_d/Aligned.Primary.Unique.BothReads.sorted.bam -O $temp_d/Aligned.Primary.Unique.BothReads.sorted.norRNA.bam --READ_LIST_FILE $temp_d/filter.txt --FILTER excludeReadList
samtools index $temp_d/Aligned.Primary.Unique.BothReads.sorted.norRNA.bam

cp $temp_d/Aligned.Primary.Multimapping.BothReads.sorted.norRNA.bam $cwd_path/$id.FivePrime_BothStrands.Primary.Multimapping.BothReads.sorted.norRNA.bam
cp $temp_d/Aligned.Primary.Multimapping.BothReads.sorted.norRNA.bam.bai $cwd_path/$id.FivePrime_BothStrands.Primary.Multimapping.BothReads.sorted.norRNA.bam.bai

cp $temp_d/Aligned.Primary.Unique.BothReads.sorted.norRNA.bam $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.BothReads.sorted.norRNA.bam
cp $temp_d/Aligned.Primary.Unique.BothReads.sorted.norRNA.bam.bai $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.BothReads.sorted.norRNA.bam.bai



echo "=========================================================="
echo "Removing rRNA mapping reads from 5' data... [ $(date) ]"
#5' multimapping
bedtools intersect -a $temp_d/Aligned.Primary.Multimapping.Read2.sorted.bam -b $temp_d/rRNA_coords.bed -sorted -g $chrom_sizes > $temp_d/rRNAreads.bam
samtools view $temp_d/rRNAreads.bam | awk '{print $1}' > $temp_d/filter.txt
picard FilterSamReads --QUIET -I $temp_d/Aligned.Primary.Multimapping.Read2.sorted.bam --O $temp_d/Aligned.Primary.Multimapping.Read2.sorted.norRNA.bam --READ_LIST_FILE $temp_d/filter.txt --FILTER excludeReadList
samtools index $temp_d/Aligned.Primary.Multimapping.Read2.sorted.norRNA.bam


#5' unique
bedtools intersect -a $temp_d/Aligned.Primary.Unique.Read2.sorted.bam -b $temp_d/rRNA_coords.bed -sorted -g $chrom_sizes > $temp_d/rRNAreads.bam
samtools view $temp_d/rRNAreads.bam | awk '{print $1}' > $temp_d/filter.txt
picard FilterSamReads --QUIET -I $temp_d/Aligned.Primary.Unique.Read2.sorted.bam -O $temp_d/Aligned.Primary.Unique.Read2.sorted.norRNA.bam --READ_LIST_FILE $temp_d/filter.txt --FILTER excludeReadList
samtools index $temp_d/Aligned.Primary.Unique.Read2.sorted.norRNA.bam

cp $temp_d/Aligned.Primary.Multimapping.Read2.sorted.norRNA.bam $cwd_path/$id.FivePrime_BothStrands.Primary.Multimapping.Read2.sorted.norRNA.bam
cp $temp_d/Aligned.Primary.Multimapping.Read2.sorted.norRNA.bam.bai $cwd_path/$id.FivePrime_BothStrands.Primary.Multimapping.Read2.sorted.norRNA.bam.bai

cp $temp_d/Aligned.Primary.Unique.Read2.sorted.norRNA.bam $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.Read2.sorted.norRNA.bam
cp $temp_d/Aligned.Primary.Unique.Read2.sorted.norRNA.bam.bai $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.Read2.sorted.norRNA.bam.bai


#Pileup
echo "=========================================================="
echo "Pileup... [ $(date) ]"
n=$(samtools view -c $temp_d/Aligned.Primary.Unique.Read2.sorted.norRNA.bam)
echo "Found $n Five Prime Reads in BOTH strands"

scale=$(bc <<< "scale=6;1000000/$n")
echo "Five_prime_reads_BOTH_strands $n\n" > $cwd_path/Read_numbers_norRNA.txt

bedtools genomecov -bga -5 -g $chrom_sizes -scale $scale -ibam $temp_d/Aligned.Primary.Unique.Read2.sorted.norRNA.bam > $temp_d/primary.unique.bdg
bedSort $temp_d/primary.unique.bdg $temp_d/primary.unique.sort.bdg
bedGraphToBigWig $temp_d/primary.unique.sort.bdg $chrom_sizes $temp_d/primary.unique.bw
mv $temp_d/primary.unique.bdg $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.Read2.sorted.norRNA.rpm.bdg
mv $temp_d/primary.unique.bw $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.Read2.sorted.norRNA.rpm.bw


echo "=========================================================="
echo "Splitting into PLUS strand... [ $(date) ]"
samtools view -F 16 -@ $(($cores - 1)) -b -o $temp_d/plus.unique.bam $temp_d/Aligned.Primary.Unique.Read2.sorted.norRNA.bam
samtools index $temp_d/plus.unique.bam

n=$(samtools view -c $temp_d/plus.unique.bam)
echo "Found $n Five Prime Reads in PLUS strand"

scale=$(bc <<< "scale=6;1000000/$n")
echo "Five_prime_reads_PLUS_strand $n\n" >> $cwd_path/Read_numbers_norRNA.txt

bedtools genomecov -bga -5 -g $chrom_sizes -scale $scale -ibam $temp_d/plus.unique.bam > $temp_d/plus.unique.bdg
bedSort $temp_d/plus.unique.bdg $temp_d/plus.unique.sort.bdg
bedGraphToBigWig $temp_d/plus.unique.sort.bdg $chrom_sizes $temp_d/plus.unique.bw

mv $temp_d/plus.unique.bam $cwd_path/$id.FivePrime_PlusStrand.Primary.Unique.Read2.sorted.norRNA.bam
mv $temp_d/plus.unique.bam.bai $cwd_path/$id.FivePrime_PlusStrand.Primary.Unique.Read2.sorted.norRNA.bam.bai
mv $temp_d/plus.unique.bw $cwd_path/$id.FivePrime_PlusStrand.Primary.Unique.Read2.sorted.norRNA.rpm.bw

echo "=========================================================="
echo "Splitting into MINUS strand... [ $(date) ]"
samtools view -f 16 -@ $(($cores - 1)) -b -o $temp_d/minus.unique.bam $temp_d/Aligned.Primary.Unique.Read2.sorted.norRNA.bam
samtools index $temp_d/minus.unique.bam

n=$(samtools view -c $temp_d/minus.unique.bam)
echo "Found $n Five Prime Reads in MINUS strand"

scale=$(bc <<< "scale=6;1000000/$n")
echo "Five_prime_reads_MINUS_strand $n\n" >> $cwd_path/Read_numbers_norRNA.txt

bedtools genomecov -bga -5 -g $chrom_sizes -scale $scale -ibam $temp_d/minus.unique.bam > $temp_d/minus.unique.bdg
bedSort $temp_d/minus.unique.bdg $temp_d/minus.unique.sort.bdg
bedGraphToBigWig $temp_d/minus.unique.sort.bdg $chrom_sizes $temp_d/minus.unique.bw

mv $temp_d/minus.unique.bam $cwd_path/$id.FivePrime_MinusStrand.Primary.Unique.Read2.sorted.norRNA.bam
mv $temp_d/minus.unique.bam.bai $cwd_path/$id.FivePrime_MinusStrand.Primary.Unique.Read2.sorted.norRNA.bam.bai
mv $temp_d/minus.unique.bw $cwd_path/$id.FivePrime_MinusStrand.Primary.Unique.Read2.sorted.norRNA.rpm.bw



echo "=========================================================="
echo "=========================================================="
echo "Removing rRNA mapping reads from ALL data... [ $(date) ]"
#all multimapping
bedtools intersect -a $temp_d/all.all.multimapping.sorted.bam -b $temp_d/rRNA_coords.bed -sorted -g $chrom_sizes > $temp_d/rRNAreads.bam
samtools view $temp_d/rRNAreads.bam | awk '{print $1}' > $temp_d/filter.txt
picard FilterSamReads --QUIET -I $temp_d/all.all.multimapping.sorted.bam -O $temp_d/all.all.multimapping.sorted.norRNA.bam --READ_LIST_FILE $temp_d/filter.txt --FILTER excludeReadList
samtools index $temp_d/all.all.multimapping.sorted.norRNA.bam

#all unique
bedtools intersect -a $temp_d/all.primary.unique.sorted.bam -b $temp_d/rRNA_coords.bed -sorted -g $chrom_sizes > $temp_d/rRNAreads.bam
samtools view $temp_d/rRNAreads.bam | awk '{print $1}' > $temp_d/filter.txt
picard FilterSamReads --QUIET -I $temp_d/all.primary.unique.sorted.bam -O $temp_d/all.primary.unique.sorted.norRNA.bam --READ_LIST_FILE $temp_d/filter.txt --FILTER excludeReadList
samtools index $temp_d/all.primary.unique.sorted.norRNA.bam

cp $temp_d/all.all.multimapping.sorted.norRNA.bam $cwd_path/$id.All.all.multimapping.sorted.norRNA.bam
cp $temp_d/all.all.multimapping.sorted.norRNA.bam.bai $cwd_path/$id.All.all.multimapping.sorted.norRNA.bam.bai

cp $temp_d/all.primary.unique.sorted.norRNA.bam $cwd_path/$id.All.unique.sorted.norRNA.bam
cp $temp_d/all.primary.unique.sorted.norRNA.bam.bai $cwd_path/$id.All.unique.sorted.norRNA.bam.bai


#Pileup
echo "=========================================================="
echo "Pileup... [ $(date) ]"
n=$(samtools view -c $temp_d/all.primary.unique.sorted.norRNA.bam)
echo "Found $n Total Reads in BOTH strands"

scale=$(bc <<< "scale=6;1000000/$n")
echo "All_Reads_BOTH_strands $n\n" >> $cwd_path/Read_numbers_norRNA.txt

bedtools genomecov -bga -split -g $chrom_sizes -scale $scale -ibam $temp_d/all.primary.unique.sorted.norRNA.bam > $temp_d/primary.unique.bdg
bedSort $temp_d/primary.unique.bdg $temp_d/primary.unique.sort.bdg
bedGraphToBigWig $temp_d/primary.unique.sort.bdg $chrom_sizes $temp_d/primary.unique.bw


echo "=========================================================="
echo "=========================================================="
echo "Removing rRNA mapping reads from ELSE data... [ $(date) ]"
#else multimapping
bedtools intersect -a $temp_d/else.all.multimapping.sorted.bam -b $temp_d/rRNA_coords.bed -sorted -g $chrom_sizes > $temp_d/rRNAreads.bam
samtools view $temp_d/rRNAreads.bam | awk '{print $1}' > $temp_d/filter.txt
picard FilterSamReads --QUIET -I $temp_d/else.all.multimapping.sorted.bam -O $temp_d/else.all.multimapping.sorted.norRNA.bam --READ_LIST_FILE $temp_d/filter.txt --FILTER excludeReadList
samtools index $temp_d/else.all.multimapping.sorted.norRNA.bam

#else unique
bedtools intersect -a $temp_d/else.primary.unique.sorted.bam -b $temp_d/rRNA_coords.bed -sorted -g $chrom_sizes > $temp_d/rRNAreads.bam
samtools view $temp_d/rRNAreads.bam | awk '{print $1}' > $temp_d/filter.txt
picard FilterSamReads --QUIET -I $temp_d/else.primary.unique.sorted.bam -O $temp_d/else.primary.unique.sorted.norRNA.bam --READ_LIST_FILE $temp_d/filter.txt --FILTER excludeReadList
samtools index $temp_d/else.primary.unique.sorted.norRNA.bam


cp $temp_d/else.all.multimapping.sorted.norRNA.bam $cwd_path/$id.Else.all.multimapping.sorted.norRNA.bam
cp $temp_d/else.all.multimapping.sorted.norRNA.bam.bai $cwd_path/$id.Else.all.multimapping.sorted.norRNA.bam.bai

cp $temp_d/else.primary.unique.sorted.norRNA.bam $cwd_path/$id.Else.unique.sorted.norRNA.bam
cp $temp_d/else.primary.unique.sorted.norRNA.bam.bai $cwd_path/$id.Else.unique.sorted.norRNA.bam.bai


#Pileup
echo "=========================================================="
echo "Pileup... [ $(date) ]"
n=$(samtools view -c $temp_d/else.primary.unique.sorted.norRNA.bam)
echo "Found $n Total Reads in BOTH strands"

scale=$(bc <<< "scale=6;1000000/$n")

bedtools genomecov -bga -split -g $chrom_sizes -scale $scale -ibam $temp_d/else.primary.unique.sorted.norRNA.bam > $temp_d/primary.unique.bdg
bedSort $temp_d/primary.unique.bdg $temp_d/primary.unique.sort.bdg
bedGraphToBigWig $temp_d/primary.unique.sort.bdg $chrom_sizes $temp_d/primary.unique.bw
mv $temp_d/primary.unique.bdg $cwd_path/$id.Else.unique.norRNA.rpm.bdg
mv $temp_d/primary.unique.bw $cwd_path/$id.Else.unique.norRNA.rpm.bw

rm -r $temp_d

echo "done";
