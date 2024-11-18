#!/bin/bash

#SBATCH -o ~/slurm_jobs/slurm_%j.log
#SBATCH -e ~/slurm_jobs/slurm_%j.err
#SBATCH -J RNAprocessing-dump
#SBATCH -p normal_q
#SBATCH -c 8
#SBATCH -t 48:00:00
#SBATCH --mem=120G

echo "SLURM Job ID        : ${SLURM_JOB_ID}"; 
echo "SLURM Job name      : ${SLURM_JOB_NAME}";
echo "SLURM Node list     : ${SLURM_JOB_NODELIST}";
echo "SLURM partition     : ${SLURM_JOB_PARTITION}";

echo "activate conda";
source activate RNA-seq;
echo "activated conda";

id=$1
cwd_path=$2
read1=$3
read2=$4
illumina_adaptors=$5
other_adaptors=$6
genome_dir=$7
chrom_sizes=$8
cores=$9
duplicate_pixel_dist=${10}
ref_genome=${11}
gtf_file=${12}
annotation_bed=${13}


echo "working on $id in $cwd_path"

mkdir -p $cwd_path
cd $cwd_path
temp_d=$(mktemp -d --tmpdir=$HOME/tmp)
echo "temporary files go in $temp_d";

echo "=========================================================="
echo "Copying fastq files... [ $(date) ]"
cp $read1 $temp_d/$id\_1.fastq.gz
cp $read2 $temp_d/$id\_2.fastq.gz

#Fastqc
echo "=========================================================="
echo "Starting FASTQC... [ $(date) ]"
fastqc $temp_d/$id\_1.fastq.gz $temp_d/$id\_2.fastq.gz -t $cores -o $cwd_path


#Clip adaptors and low quality ends
echo "=========================================================="
echo "Trimming adaptors with trimmomatic... [ $(date) ]"
cp $illumina_adaptors $temp_d/adaptors.fa
cp $other_adaptors $temp_d/other_adaptors.fa
trimmomatic PE $temp_d/$id\_1.fastq.gz $temp_d/$id\_2.fastq.gz -baseout $temp_d/$id.clean.fastq.gz -threads $cores ILLUMINACLIP:$temp_d/adaptors.fa:2:30:10:1:true ILLUMINACLIP:$temp_d/other_adaptors.fa:2:30:10:1:true MINLEN:30 &>> $cwd_path/QC_trimmomatic.txt
rm $temp_d/$id\_1.fastq.gz $temp_d/$id\_2.fastq.gz $temp_d/$id.clean_1U.fastq.gz $temp_d/$id.clean_2U.fastq.gz

#Clip adaptors and low quality ends
echo "=========================================================="
echo "Trimming adaptors and splitting based on tags... [ $(date) ]"
adaptor_root="AAGCAGTGGTATCAACGCAGAGTAC"
adaptor_five="AAGCAGTGGTATCAACGCAGAGTACATGGG"
adaptor_three="AAGCAGTGGTATCAACGCAGAGTACTTTTT"
launchTailedRNAseqProcessing.sh $temp_d/$id.clean_1P.fastq.gz $temp_d/$id.clean_2P.fastq.gz $adaptor_root $adaptor_five $adaptor_three $temp_d/$id\_five $temp_d/$id\_three $temp_d/$id\_else $cores

echo "=========================================================="
echo "Copying split files... [ $(date) ]"

cat $temp_d/$id\_five_1.fastq.gz > $temp_d/$id\_all_1.fastq.gz
cat $temp_d/$id\_five_2.fastq.gz > $temp_d/$id\_all_2.fastq.gz
cat $temp_d/$id\_three_1.fastq.gz >> $temp_d/$id\_all_1.fastq.gz
cat $temp_d/$id\_three_2.fastq.gz >> $temp_d/$id\_all_2.fastq.gz
cat $temp_d/$id\_else_1.fastq.gz >> $temp_d/$id\_all_1.fastq.gz
cat $temp_d/$id\_else_2.fastq.gz >> $temp_d/$id\_all_2.fastq.gz

cp $temp_d/$id\_five_1.fastq.gz $cwd_path
cp $temp_d/$id\_five_2.fastq.gz $cwd_path
cp $temp_d/$id\_three_1.fastq.gz $cwd_path
cp $temp_d/$id\_three_2.fastq.gz $cwd_path
cp $temp_d/$id\_else_1.fastq.gz $cwd_path
cp $temp_d/$id\_else_2.fastq.gz $cwd_path
cp $temp_d/$id\_all_1.fastq.gz $cwd_path
cp $temp_d/$id\_all_2.fastq.gz $cwd_path

#Clip adaptors and low quality ends
echo "=========================================================="
echo "=========================================================="
echo "Starting 5' processing... [ $(date) ]"
echo "=========================================================="
echo "Trimming adaptors..."
trimmomatic PE $temp_d/$id\_five_1.fastq.gz $temp_d/$id\_five_2.fastq.gz -quiet -baseout $temp_d/$id.five_clean.fastq.gz -threads $cores TRAILING:10 MINLEN:30 &>> $cwd_path/QC_trimmomatic_5prime.txt
rm $temp_d/$id\_five_1.fastq.gz $temp_d/$id\_five_2.fastq.gz

#align using STAR
echo "=========================================================="
echo "Aligning... [ $(date) ]"
#Alignment: no mismatch filtering. will keep only best scoring in subsequent step.
#Multimappers: all multimapping reads will be output, and only those with the top score will be considered primary alignments. Reads with more than 20 alignments will be discarded.
#Output: not sorted, read pairs are adjacent
STAR --runMode alignReads --runThreadN $cores --genomeDir $genome_dir --readFilesIn $temp_d/$id.five_clean_1P.fastq.gz $temp_d/$id.five_clean_2P.fastq.gz --readFilesCommand zcat --outFilterType BySJout --outFilterMultimapNmax 100 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix $temp_d/ --outMultimapperOrder Random --outSAMprimaryFlag AllBestScore --outSAMstrandField intronMotif --runRNGseed 23 --outSAMtype BAM Unsorted --quantMode GeneCounts --twopassMode Basic
mv $temp_d/Log.final.out $cwd_path/$id.Five_prime_Log.final.out
mv $temp_d/ReadsPerGene.out.tab $cwd_path/$id.Five_prime_ReadsPerGene.out.tab
mv $temp_d/SJ.out.tab $cwd_path/$id.Five_prime.SJ.out.tab
mv $temp_d/Log.out $cwd_path/$id.Five_prime.Log.out
rm -r $temp_d/_STARgenome $temp_d/_STARpass1 $temp_d/$id.five_clean_1P.fastq.gz $temp_d/$id.five_clean_2P.fastq.gz

echo "=========================================================="
echo "Removing secondary alignments, singletons and unmapped reads... [ $(date) ]"
samtools view -F 268 -@ $(($cores - 1)) $temp_d/Aligned.out.bam -o $temp_d/Aligned.Primary.Multimapping.bam

echo "=========================================================="
echo "Sorting, indexing and processing both stranded files... [ $(date) ]"
samtools sort -T $temp_d -@ $(($cores - 1)) -o $temp_d/Aligned.Primary.Multimapping.BothReads.sorted.bam $temp_d/Aligned.Primary.Multimapping.bam
samtools index $temp_d/Aligned.Primary.Multimapping.BothReads.sorted.bam
samtools view -q 255 -@ $(($cores - 1)) -b -o $temp_d/Aligned.Primary.Unique.BothReads.sorted.bam $temp_d/Aligned.Primary.Multimapping.BothReads.sorted.bam
samtools index $temp_d/Aligned.Primary.Unique.BothReads.sorted.bam
cp $temp_d/Aligned.Primary.Multimapping.BothReads.sorted.bam $cwd_path/$id.FivePrime_BothStrands.Primary.Multimapping.BothStrands.sorted.bam
cp $temp_d/Aligned.Primary.Multimapping.BothReads.sorted.bam.bai $cwd_path/$id.FivePrime_BothStrands.Primary.Multimapping.BothStrands.sorted.bam.bai
cp $temp_d/Aligned.Primary.Unique.BothReads.sorted.bam $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.BothStrands.sorted.bam
cp $temp_d/Aligned.Primary.Unique.BothReads.sorted.bam.bai $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.BothStrands.sorted.bam.bai

echo "=========================================================="
echo "Keeping only Read #2... [ $(date) ]"
#Multimapping bam includes all alignments, not only primary ones
samtools view -f 128 -@ $(($cores - 1)) -b -o $temp_d/Aligned.Primary.Multimapping.Read2.bam $temp_d/Aligned.Primary.Multimapping.bam

echo "=========================================================="
echo "Sorting and indexing... [ $(date) ]"
samtools sort -T $temp_d -@ $(($cores - 1)) -o $temp_d/Aligned.Primary.Multimapping.Read2.sorted.bam $temp_d/Aligned.Primary.Multimapping.Read2.bam
samtools index $temp_d/Aligned.Primary.Multimapping.Read2.sorted.bam
cp $temp_d/Aligned.Primary.Multimapping.Read2.sorted.bam $cwd_path/$id.FivePrime_BothStrands.Primary.Multimapping.Read2.sorted.bam
cp $temp_d/Aligned.Primary.Multimapping.Read2.sorted.bam.bai $cwd_path/$id.FivePrime_BothStrands.Primary.Multimapping.Read2.sorted.bam.bai
rm $temp_d/Aligned.out.bam $temp_d/Aligned.Primary.Multimapping.Read2.bam

######starting unique analysis#####
echo "=========================================================="
echo "Filtering out multimapping... [ $(date) ]"
samtools view -q 255 -@ $(($cores - 1)) -b -o $temp_d/Aligned.Primary.Unique.Read2.sorted.bam $temp_d/Aligned.Primary.Multimapping.Read2.sorted.bam
samtools index $temp_d/Aligned.Primary.Unique.Read2.sorted.bam
cp $temp_d/Aligned.Primary.Unique.Read2.sorted.bam $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.Read2.sorted.bam
cp $temp_d/Aligned.Primary.Unique.Read2.sorted.bam.bai $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.Read2.sorted.bam.bai


#Pileup
echo "=========================================================="
echo "Pileup... [ $(date) ]"
n=$(samtools view -c $temp_d/Aligned.Primary.Unique.Read2.sorted.bam)
echo "Found $n Five Prime Reads in BOTH strands"

scale=$(bc <<< "scale=6;1000000/$n")
echo "Five_prime_reads_BOTH_strands $n\n" > $cwd_path/Read_numbers.txt

bedtools genomecov -bga -5 -scale $scale -ibam $temp_d/Aligned.Primary.Unique.Read2.sorted.bam > $temp_d/primary.unique.bdg
bedSort $temp_d/primary.unique.bdg $temp_d/primary.unique.sort.bdg
bedGraphToBigWig $temp_d/primary.unique.sort.bdg $chrom_sizes $temp_d/primary.unique.bw
mv $temp_d/primary.unique.bdg $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.Read2.sorted.rpm.bdg
mv $temp_d/primary.unique.bw $cwd_path/$id.FivePrime_BothStrands.Primary.Unique.Read2.sorted.rpm.bw

echo "=========================================================="
echo "Splitting into PLUS strand... [ $(date) ]"
samtools view -F 16 -@ $(($cores - 1)) -b -o $temp_d/plus.unique.bam $temp_d/Aligned.Primary.Unique.Read2.sorted.bam
samtools index $temp_d/plus.unique.bam

n=$(samtools view -c $temp_d/plus.unique.bam)
echo "Found $n Five Prime Reads in PLUS strand"

scale=$(bc <<< "scale=6;1000000/$n")
echo "Five_prime_reads_PLUS_strand $n\n" >> $cwd_path/Read_numbers.txt

bedtools genomecov -bga -5 -scale $scale -ibam $temp_d/plus.unique.bam > $temp_d/plus.unique.bdg
bedSort $temp_d/plus.unique.bdg $temp_d/plus.unique.sort.bdg
bedGraphToBigWig $temp_d/plus.unique.sort.bdg $chrom_sizes $temp_d/plus.unique.bw

mv $temp_d/plus.unique.bam $cwd_path/$id.FivePrime_PlusStrand.Primary.Unique.Read2.sorted.bam
mv $temp_d/plus.unique.bam.bai $cwd_path/$id.FivePrime_PlusStrand.Primary.Unique.Read2.sorted.bam.bai
mv $temp_d/plus.unique.bw $cwd_path/$id.FivePrime_PlusStrand.Primary.Unique.Read2.sorted.rpm.bw

echo "=========================================================="
echo "Splitting into MINUS strand... [ $(date) ]"
samtools view -f 16 -@ $(($cores - 1)) -b -o $temp_d/minus.unique.bam $temp_d/Aligned.Primary.Unique.Read2.sorted.bam
samtools index $temp_d/minus.unique.bam

n=$(samtools view -c $temp_d/minus.unique.bam)
echo "Found $n Five Prime Reads in MINUS strand"

scale=$(bc <<< "scale=6;1000000/$n")
echo "Five_prime_reads_MINUS_strand $n\n" >> $cwd_path/Read_numbers.txt

bedtools genomecov -bga -5 -scale $scale -ibam $temp_d/minus.unique.bam > $temp_d/minus.unique.bdg
bedSort $temp_d/minus.unique.bdg $temp_d/minus.unique.sort.bdg
bedGraphToBigWig $temp_d/minus.unique.sort.bdg $chrom_sizes $temp_d/minus.unique.bw

mv $temp_d/minus.unique.bam $cwd_path/$id.FivePrime_MinusStrand.Primary.Unique.Read2.sorted.bam
mv $temp_d/minus.unique.bam.bai $cwd_path/$id.FivePrime_MinusStrand.Primary.Unique.Read2.sorted.bam.bai
mv $temp_d/minus.unique.bw $cwd_path/$id.FivePrime_MinusStrand.Primary.Unique.Read2.sorted.rpm.bw


#Clip adaptors and low quality ends
echo "=========================================================="
echo "=========================================================="
echo "Starting 'else' processing... [ $(date) ]"
echo "=========================================================="
echo "Trimming adaptors... [ $(date) ]"
trimmomatic PE $temp_d/$id\_else_1.fastq.gz $temp_d/$id\_else_2.fastq.gz -quiet -baseout $temp_d/$id.else_clean.fastq.gz -threads $cores TRAILING:10 MINLEN:30 &>> $cwd_path/QC_trimmomatic_else.txt
rm $temp_d/$id\_else_1.fastq.gz $temp_d/$id\_else_2.fastq.gz

#align using STAR
echo "=========================================================="
echo "Aligning... [ $(date) ]"
STAR --runMode alignReads --runThreadN $cores --genomeDir $genome_dir --readFilesIn $temp_d/$id.else_clean_1P.fastq.gz $temp_d/$id.else_clean_2P.fastq.gz --readFilesCommand zcat --outFilterType BySJout --outFilterMultimapNmax 100 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix $temp_d/ --outMultimapperOrder Random --outSAMprimaryFlag AllBestScore --outSAMstrandField intronMotif --runRNGseed 23 --outSAMtype BAM Unsorted --quantMode GeneCounts --twopassMode Basic
mv $temp_d/Log.final.out $cwd_path/$id.Else_Log.final.out
mv $temp_d/ReadsPerGene.out.tab $cwd_path/$id.Else_ReadsPerGene.out.tab
mv $temp_d/SJ.out.tab $cwd_path/$id.Else.SJ.out.tab
mv $temp_d/Log.out $cwd_path/$id.Else.Log.out
rm -r $temp_d/_STARgenome $temp_d/_STARpass1 $temp_d/$id.else_clean_1P.fastq.gz $temp_d/$id.else_clean_2P.fastq.gz



echo "=========================================================="
echo "Removing secondary alignments, singletons and unmapped reads... [ $(date) ]"
samtools view -F 268 -@ $(($cores - 1)) $temp_d/Aligned.out.bam -o $temp_d/Aligned.Primary.Multimapping.bam

#########################################
#######starting unique analysis##########
echo "=========================================================="
echo "Splitting into uniquely-mapping and multimapping... [ $(date) ]"
#Multimapping bam includes all alignments, not only primary ones
samtools view -q 255 -@ $(($cores - 1)) -b -o $temp_d/primary.unique.bam $temp_d/Aligned.Primary.Multimapping.bam
mv $temp_d/Aligned.Primary.Multimapping.bam $temp_d/all.multimapping.bam

echo "=========================================================="
echo "Sorting... [ $(date) ]"
samtools sort -T $temp_d -@ $(($cores - 1)) -o $temp_d/primary.unique.sorted.bam $temp_d/primary.unique.bam
rm $temp_d/primary.unique.bam
cp $temp_d/primary.unique.sorted.bam $cwd_path/$id.Else.unique.sorted.bam

#Indexing
echo "=========================================================="
echo "Indexing... [ $(date) ]"
samtools index $temp_d/primary.unique.sorted.bam
cp $temp_d/primary.unique.sorted.bam.bai $cwd_path/$id.Else.unique.sorted.bam.bai

#Pileup
echo "=========================================================="
echo "Pileup... [ $(date) ]"

n=$(samtools view -c $temp_d/primary.unique.sorted.bam)
echo "Found $n Internal Reads in BOTH strands"

scale=$(bc <<< "scale=6;1000000/$n")
echo "Internal_Reads_BOTH_strands $n\n" >> $cwd_path/Read_numbers.txt

bedtools genomecov -bga -split -scale $scale -ibam $temp_d/primary.unique.sorted.bam > $temp_d/primary.unique.bdg
bedSort $temp_d/primary.unique.bdg $temp_d/primary.unique.sort.bdg
bedGraphToBigWig $temp_d/primary.unique.sort.bdg $chrom_sizes $temp_d/primary.unique.bw
mv $temp_d/primary.unique.bdg $cwd_path/$id.Else.unique.rpm.bdg
mv $temp_d/primary.unique.bw $cwd_path/$id.Else.unique.rpm.bw


#########################################
######starting multimapping analysis#####
echo "=========================================================="
echo "Sorting... [ $(date) ]"
samtools sort -T $temp_d -@ $(($cores - 1)) -o $temp_d/all.multimapping.sorted.bam $temp_d/all.multimapping.bam
rm $temp_d/all.multimapping.bam
cp $temp_d/all.multimapping.sorted.bam $cwd_path/$id.Else.all.multimapping.sorted.bam

#Indexing
echo "=========================================================="
echo "Indexing.. [ $(date) ]"
samtools index $temp_d/all.multimapping.sorted.bam
cp $temp_d/all.multimapping.sorted.bam.bai $cwd_path/$id.Else.all.multimapping.sorted.bam.bai


#Clip adaptors and low quality ends
echo "=========================================================="
echo "=========================================================="
echo "Starting 'All' processing... [ $(date) ]"
echo "=========================================================="
echo "Trimming adaptors... [ $(date) ]"
trimmomatic PE $temp_d/$id\_all_1.fastq.gz $temp_d/$id\_all_2.fastq.gz -quiet -baseout $temp_d/$id.all_clean.fastq.gz -threads $cores TRAILING:10 MINLEN:30 &>> $cwd_path/QC_trimmomatic_all.txt
rm $temp_d/$id\_all_1.fastq.gz $temp_d/$id\_all_2.fastq.gz

#align using STAR
echo "=========================================================="
echo "Aligning... [ $(date) ]"
STAR --runMode alignReads --runThreadN $cores --genomeDir $genome_dir --readFilesIn $temp_d/$id.all_clean_1P.fastq.gz $temp_d/$id.all_clean_2P.fastq.gz --readFilesCommand zcat --outFilterType BySJout --outFilterMultimapNmax 100 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix $temp_d/ --outMultimapperOrder Random --outSAMprimaryFlag AllBestScore --outSAMstrandField intronMotif --runRNGseed 23 --outSAMtype BAM Unsorted --quantMode GeneCounts --twopassMode Basic
mv $temp_d/Log.final.out $cwd_path/$id.All_Log.final.out
mv $temp_d/ReadsPerGene.out.tab $cwd_path/$id.All_ReadsPerGene.out.tab
mv $temp_d/SJ.out.tab $cwd_path/$id.All.SJ.out.tab
mv $temp_d/Log.out $cwd_path/$id.All.Log.out
rm -r $temp_d/_STARgenome $temp_d/_STARpass1 $temp_d/$id.all_clean_1P.fastq.gz $temp_d/$id.all_clean_2P.fastq.gz

#echo "=========================================================="
#echo "Sorting and indexing... [ $(date) ]"
#samtools sort -n -@ $(($cores - 1)) -o $temp_d/Aligned.out.sorted.bam $temp_d/Aligned.out.bam

echo "=========================================================="
echo "Removing secondary alignments, singletons and unmapped reads... [ $(date) ]"
samtools view -F 268 -@ $(($cores - 1)) $temp_d/Aligned.out.bam -o $temp_d/Aligned.Primary.Multimapping.bam

#########################################
#######starting unique analysis##########
echo "=========================================================="
echo "Splitting into uniquely-mapping and multimapping... [ $(date) ]"
#Multimapping bam includes all alignments, not only primary ones
samtools view -q 255 -@ $(($cores - 1)) -b -o $temp_d/primary.unique.bam $temp_d/Aligned.Primary.Multimapping.bam
mv $temp_d/Aligned.Primary.Multimapping.bam $temp_d/all.multimapping.bam

echo "=========================================================="
echo "Sorting... [ $(date) ]"
samtools sort -T $temp_d -@ $(($cores - 1)) -o $temp_d/primary.unique.sorted.bam $temp_d/primary.unique.bam
rm $temp_d/primary.unique.bam
cp $temp_d/primary.unique.sorted.bam $cwd_path/$id.All.unique.sorted.bam

#Indexing
echo "=========================================================="
echo "Indexing... [ $(date) ]"
samtools index $temp_d/primary.unique.sorted.bam
cp $temp_d/primary.unique.sorted.bam.bai $cwd_path/$id.All.unique.sorted.bam.bai

#Pileup
echo "=========================================================="
echo "Pileup... [ $(date) ]"
n=$(samtools view -c $temp_d/primary.unique.sorted.bam)
echo "Found $n Total Reads in BOTH strands"

scale=$(bc <<< "scale=6;1000000/$n")
echo "All_Reads_BOTH_strands $n\n" >> $cwd_path/Read_numbers.txt

bedtools genomecov -bga -split -scale $scale -ibam $temp_d/primary.unique.sorted.bam > $temp_d/primary.unique.bdg
bedSort $temp_d/primary.unique.bdg $temp_d/primary.unique.sort.bdg
bedGraphToBigWig $temp_d/primary.unique.sort.bdg $chrom_sizes $temp_d/primary.unique.bw
mv $temp_d/primary.unique.bdg $cwd_path/$id.All.unique.rpm.bdg
mv $temp_d/primary.unique.bw $cwd_path/$id.All.unique.rpm.bw


#########################################
######starting multimapping analysis#####
echo "=========================================================="
echo "Sorting... [ $(date) ]"
samtools sort -T $temp_d -@ $(($cores - 1)) -o $temp_d/all.multimapping.sorted.bam $temp_d/all.multimapping.bam
rm $temp_d/all.multimapping.bam
cp $temp_d/all.multimapping.sorted.bam $cwd_path/$id.All.all.multimapping.sorted.bam

#Indexing
echo "=========================================================="
echo "Indexing... [ $(date) ]"
samtools index $temp_d/all.multimapping.sorted.bam
cp $temp_d/all.multimapping.sorted.bam.bai $cwd_path/$id.All.all.multimapping.sorted.bam.bai

#Pileup
echo "=========================================================="
echo "Pileup... [ $(date) ]"
n=$(samtools view -c $temp_d/all.multimapping.sorted.bam)

scale=$(bc <<< "scale=6;1000000/$n")

bedtools genomecov -bga -split -scale $scale -ibam $temp_d/all.multimapping.sorted.bam > $temp_d/primary.multimapping.bdg
bedSort $temp_d/primary.multimapping.bdg $temp_d/primary.multimapping.sort.bdg 
bedGraphToBigWig $temp_d/primary.multimapping.sort.bdg $chrom_sizes $temp_d/primary.multimapping.bw
mv $temp_d/primary.multimapping.bdg $cwd_path/$id.All.multimapping.rpm.bdg
mv $temp_d/primary.multimapping.bw $cwd_path/$id.All.multimapping.rpm.bw

#########
echo "=========================================================="
echo "=========================================================="
echo "Running statistics on uniquely mapping reads with 'All' reads... [ $(date) ]"
cp $temp_d/primary.unique.sorted.bam $temp_d/$id.bam
cp $temp_d/primary.unique.sorted.bam.bai $temp_d/$id.bam.bai

samtools view -F 8 -o $temp_d/$id.fixed.bam $temp_d/$id.bam
rm $temp_d/$id.bam
mv $temp_d/$id.fixed.bam $temp_d/$id.bam
samtools index $temp_d/$id.bam

echo "=========================================================="
echo "Fixing mates for Picard... [ $(date) ]"
picard FixMateInformation --QUIET true -I $temp_d/$id.bam

echo "=========================================================="
echo "Counting duplicates, they will not be removed... NOTE: Update optical values for Hiseq4000 generated data !!!  [ $(date) ]"
picard MarkDuplicates --QUIET true -I $temp_d/$id.bam -O $temp_d/primary.unique.sorted.counted.bam -M $cwd_path/QC_Picard_duplicates.txt -REMOVE_DUPLICATES false -OPTICAL_DUPLICATE_PIXEL_DISTANCE $duplicate_pixel_dist

echo "=========================================================="
echo "Counting GC bias...  [ $(date) ]"
picard CollectGcBiasMetrics --QUIET true -I $temp_d/$id.bam -O $cwd_path/QC_Picard_GCbias_metrics.txt -S $cwd_path/QC_Picard_GCbias_summary.txt -CHART_OUTPUT $cwd_path/QC_Picard_GCchart.pdf -R $ref_genome

echo "=========================================================="
echo "Calculating insert sizes... [ $(date) ]"
picard CollectInsertSizeMetrics --QUIET true -I $temp_d/$id.bam -O $cwd_path/QC_Picard_insert_metrics.txt -H $cwd_path/QC_Picard_insert_chart.pdf


echo "=========================================================="
echo "Calculating complexity... [ $(date) ]"
preseq c_curve -B -v -o $cwd_path/$id.preseq.txt $temp_d/$id.bam


echo "=========================================================="
echo "Calculating chromosome stats... [ $(date) ]"
samtools idxstats $temp_d/$id.bam > $cwd_path/$id.idxstat

rm $temp_d/primary.unique.sorted.bam $temp_d/primary.unique.sorted.bam.bai $temp_d/$id.bam $temp_d/$id.bam.bai

rm -rf $temp_d;

echo "done";






