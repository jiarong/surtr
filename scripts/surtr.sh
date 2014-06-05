#! /usr/bin/env bash

# stop when fail

set -e
set -o pipefail 

if [ $# -ne 3 ]
then
  echo "Usage: $(basename $0) <ref> <seqfile> <outdir>"
  exit 1
fi

Ref=$1
Seq=$2
Outdir=$3


#Prefix=5M
Prefix=$(basename $Seq .fasta).surtr
Log="$Outdir"/"$Prefix".log


#module load bwa
Scriptpath=$( cd "$(dirname "$0")" && pwd -P ) 
Binpath=$( cd "$(dirname "$0")" && cd ../bin && pwd -P ) 
Dbpath=$( cd "$(dirname "$0")" && cd ../db && pwd -P ) 
Libpath=$( cd "$(dirname "$0")" && cd ../lib && pwd -P ) 

#Ref=$Dbpath/rrna.fasta

export PYTHONPATH=$Libpath:$PYTHONPATH

#python -c "import khmer; print 'khmer version: {}'.format(khmer.__version__)"

if [ ! -f $Ref.bwt ];
then
  $Binpath/bwa index $Ref
fi

mkdir -p $Outdir
date|tee -a $Log
echo 'start surtr'|tee -a $Log
time(
python $Scriptpath/sweep_read.py -k 21 -N 4 -x 1e8 -P 10 $Ref $Seq $Outdir/$Prefix.sweep.nonrrna - | \
$Binpath/bwa mem -t 1 -c 10000000 $Ref - | \
python $Scriptpath/parse_sam.py -o $Outdir -p $Prefix.parsesam -c 75 -i 50 -

echo "                               " | tee -a $Log  # donot konw why have to
) 2>&1 | tee -a $Log
echo 'surtr finished'|tee -a $Log
echo |tee -a $Log
echo |tee -a $Log

cat $Outdir/$Prefix.sweep.nonrrna $Outdir/$Prefix.parsesam.nonrrna > $Outdir/"$Prefix"_nonrrna.fasta

cat $Outdir/$Prefix.parsesam.rrna $Outdir/$Prefix.parsesam.nonrrna > $Outdir/"$Prefix"_sweep.rrna

mv $Outdir/$Prefix.sweep.nonrrna $Outdir/"$Prefix"_sweep.nonrrna
mv $Outdir/$Prefix.parsesam.rrna  $Outdir/"$Prefix"_rrna.fasta
mv $Outdir/$Prefix.parsesam.csv  $Outdir/"$Prefix"_rrna.csv

rm -f $Outdir/$Prefix.parsesam.rrna $Outdir/$Prefix.parsesam.nonrrna $Outdir/$Prefix.parsesam.temp  $Outdir/$Prefix.sweep.nonrrna
rm -f $Outdir/"$Prefix"_sweep.nonrrna $Outdir/"$Prefix"_sweep.rrna
