file1=$1
file2=$2
adaptor_root=$3
adaptor_five=$4
adaptor_three=$5
final_output_five=$6
final_output_three=$7
final_output_else=$8
cores=$9

#note that fastqs should already be in a temp folder
temp_d=$(mktemp -d --tmpdir=$HOME/tmp)
echo "temporary files go to $temp_d";

#count number of reads in one fastq
n=`zcat $file1 | wc -l`
n=${n[0]}
n=$((n / 4))

#one core less than assigned will be used for fastq processing
usable_cores=$(($cores - 1))

#numer of reads processed per core rounded up
x=$(($(($n + $usable_cores - 1))/$usable_cores))

#arrays harboring temp file names
declare -a list_five_files
declare -a list_three_files
declare -a list_else_files

#decompress files once for speed
file1_decomp=$(mktemp --tmpdir=$temp_d)
file2_decomp=$(mktemp --tmpdir=$temp_d)
zcat $file1 > $file1_decomp
zcat $file2 > $file2_decomp

#set pids list
pids=""
RESULT=0

echo "Parameters loaded. Launching $usable_cores processes, processing $x reads each... [ $(date) ]"

for i in `seq 0 $(($usable_cores-1))`; do
  initial=$(($(($i*$x))+1))
  if [ "$i" -eq $(($usable_cores-1)) ]; then
    last=$n
  else
    last=$(($(($i+1))*$x))
  fi

  temp1=$(mktemp --tmpdir=$temp_d)
  cat $file1_decomp | head -n $(($last*4)) | tail -n $(($(($last-$initial+1))*4)) > ${temp1}

  temp2=$(mktemp --tmpdir=$temp_d)
  cat $file2_decomp | head -n $(($last*4)) | tail -n $(($(($last-$initial+1))*4)) > ${temp2}

  out_five=$(mktemp --tmpdir=$temp_d)
  out_three=$(mktemp --tmpdir=$temp_d)
  out_else=$(mktemp --tmpdir=$temp_d)

  list_five_files[i]=$out_five
  list_three_files[i]=$out_three
  list_else_files[i]=$out_else

  echo "Values are n = $n x = $x cores = $usable_cores initial = $initial final = $last"

  python processTailedRNAseq.py $temp1 $temp2 $adaptor_root $adaptor_five $adaptor_three $out_five $out_three $out_else &

  pids="$pids $!"
done

echo "=========================================================="
echo "Waiting for processes to finish... [ $(date) ]"

#wait for processes to finish
for pid in $pids; do
   wait $pid || let "RESULT=1"
done

if [ "$RESULT" == "1" ];
   then
     echo "=========================================================="
     echo "ERROR!!"
     echo "Some instance failed to execute!"
     exit 1
fi

echo "=========================================================="
echo "Merging files... [ $(date) ]"

#merge files
for i in `seq 0 $(($usable_cores-1))`; do
  cat "${list_five_files[$i]}_1.fastq" | pigz -p $cores >> "${final_output_five}_1.fastq.gz"
  cat "${list_five_files[$i]}_2.fastq" | pigz -p $cores >> "${final_output_five}_2.fastq.gz"
  cat "${list_three_files[$i]}_1.fastq" | pigz -p $cores >> "${final_output_three}_1.fastq.gz"
  cat "${list_three_files[$i]}_2.fastq" | pigz -p $cores >> "${final_output_three}_2.fastq.gz"
  cat "${list_else_files[$i]}_1.fastq" | pigz -p $cores >> "${final_output_else}_1.fastq.gz"
  cat "${list_else_files[$i]}_2.fastq" | pigz -p $cores >> "${final_output_else}_2.fastq.gz"

rm -rf $temp_d
done;
