#! /usr/bin/env bash

# stop when fail
set -e
set -o pipefail 

export PYTHONPATH=/home/gjr/data/SSUsearch/compareTools/curtis/lib:$PYTHONPATH
python -c "import khmer; print khmer.__version__"

#Seq=/mnt/research/tiedjelab/temp/samples/SSUsearch/data/MSB1_LN1.fastq.bz2.kmerout/MSB1_LN1.afterMerge.fa
Seq=/home/gjr/data/SSUsearch/compareTools/data/test.5M.fa
#Seq=/home/gjr/data/SSUsearch/compareTools/data/test.5M.fa.ribopickerout/ribopick_rrna.fa
Ref=/home/gjr/software/RNA/ribopicker-standalone-0.4.3/db/ssr108_qi_rep.fa
Outdir=check_ribopicker
Prefix=$(basename $Seq .fasta).curtis_nosweep
Log="$Prefix".log

#module load bwa
Scriptpath=$( cd "$(dirname "$0")" && pwd -P ) 
Binpath=$( cd "$(dirname "$0")" && cd ../bin && pwd -P ) 
Dbpath=$( cd "$(dirname "$0")" && cd ../db && pwd -P ) 

mkdir -p $Outdir

date|tee -a $Log
echo 'start curtis_nosweep'|tee -a $Log
#python sweep_read.py -k 21 -N 4 -x 1e9 -P 10 $Ref $Seq $Outdir/$Prefix.sweep.nonrrna - | \
time(
cat $Seq | \
$Binpath/bwa mem -t 1 -c 10000000 $Ref - | \
python $Scriptpath/parse_sam.py -o $Outdir -p $Prefix.parsesam -c 75 -i 50 -

echo "                               " | tee -a $Log  # donot konw why have to
) 2>&1 | tee -a $Log
echo 'curtis_nosweep finished'|tee -a $Log
echo |tee -a $Log
echo |tee -a $Log

cat $Outdir/$Prefix.sweep.nonrrna $Outdir/$Prefix.parsesam.nonrrna > $Outdir/$Prefix.nonrrna
#mv $Outdir/$Prefix.parsesam.rrna > $Outdir/$Prefix.rrna

