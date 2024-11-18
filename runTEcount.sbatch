#!/bin/bash

#SBATCH -o ~/slurm_jobs/slurm_%j.log
#SBATCH -e ~/slurm_jobs/slurm_%j.err
#SBATCH -J TEcount
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH -c 4
#SBATCH -t 48:00:00
#SBATCH --mem=100G

echo "SLURM Job ID        : ${SLURM_JOB_ID}";
echo "SLURM Job name      : ${SLURM_JOB_NAME}";
echo "SLURM Node list     : ${SLURM_JOB_NODELIST}";
echo "SLURM partition     : ${SLURM_JOB_PARTITION}";

echo "activate conda";
source activate TEtranscripts;
echo "activated conda";

cwd=$1
gtf=$2
te_gtf=$3


id=$(basename $cwd)
temp_d=$(mktemp -d --tmpdir=$HOME/tmp)
echo "temporary files go to $temp_d"
cd $temp_d

echo "TEcount for 5' reads";
echo "convert Read2 to proper single end file";
perl convert_pebamRead2_single_end.pl $cwd/$id.FivePrime_BothStrands.Primary.Multimapping.Read2.sorted.norRNA.bam;
samtools view -b $cwd/$id.FivePrime_BothStrands.Primary.Multimapping.Read2.sorted.norRNA.bam.single_end.sam.gz -o $cwd/$id.FivePrime_BothStrands.Primary.Multimapping.Read2.sorted.norRNA.bam.single_end.bam;
rm $cwd/$id.FivePrime_BothStrands.Primary.Multimapping.Read2.sorted.norRNA.bam.single_end.sam.gz;
TEcount -b $cwd/$id.FivePrime_BothStrands.Primary.Multimapping.Read2.sorted.norRNA.bam.single_end.bam -L 115 --GTF $gtf --TE $te_gtf --sortByPos --mode multi --project $cwd/$id.TEcount.dfam.FivePrime.Read2 --verbose 0

rm -r $temp_d

echo "done";


