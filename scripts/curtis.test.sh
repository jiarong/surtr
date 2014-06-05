#! /usr/bin/env bash

# stop when fail
set -e
set -o pipefail 

#Seq=/home/gjr/data/SSUsearch/compareTools/data/simulated_data/mg_sim_l200_e0_1.fasta
#Seq=/home/gjr/data/SSUsearch/compareTools/data/simulated_data/mg_sim_l200_e2_1.fasta
#Seq=/home/gjr/data/SSUsearch/compareTools/data/simulated_data/mg_sim_l200_e5_1.fasta
Seq=/home/gjr/data/SSUsearch/compareTools/data/test.5M.fa

Outdir=test_curtis
#Outdir=simu_out_rrna
Ref=/home/gjr/software/RNA/ribopicker-standalone-0.4.3/db/ssr108_qi_rep.fa
#Ref=/home/gjr/data/SSUsearch/compareTools/curtis/db/rrna.fasta

#Prefix=5M
Prefix=$(basename $Seq .fasta).curtis
Log="$Prefix".log

export PYTHONPATH=/home/gjr/data/SSUsearch/compareTools/curtis/lib:$PYTHONPATH
python -c "import khmer; print 'khmer version: {}'.format(khmer.__version__)"

#module load bwa
Scriptpath=$( cd "$(dirname "$0")" && pwd -P ) 
Binpath=$( cd "$(dirname "$0")" && cd ../bin && pwd -P ) 
Dbpath=$( cd "$(dirname "$0")" && cd ../db && pwd -P ) 

mkdir -p $Outdir
python $Scriptpath/sweep_read.py -k 21 -N 4 -x 1e8 -P 10 $Ref $Seq $Outdir/$Prefix.sweep.nonrrna - | \
$Binpath/bwa mem -t 1 -c 10000000 $Ref - | \
python $Scriptpath/parse_sam.py -o $Outdir -p $Prefix.parsesam -c 75 -i 50 -

cat $Outdir/$Prefix.sweep.nonrrna $Outdir/$Prefix.parsesam.nonrrna > $Outdir/"$Prefix"_nonrrna.fasta

cat $Outdir/$Prefix.parsesam.rrna $Outdir/$Prefix.parsesam.nonrrna > $Outdir/"$Prefix"_sweep.rrna

mv $Outdir/$Prefix.sweep.nonrrna $Outdir/"$Prefix"_sweep.nonrrna
mv $Outdir/$Prefix.parsesam.rrna  $Outdir/"$Prefix"_rrna.fasta
mv $Outdir/$Prefix.parsesam.csv  $Outdir/"$Prefix"_rrna.csv

rm -f $Outdir/$Prefix.parsesam.rrna $Outdir/$Prefix.parsesam.nonrrna $Outdir/$Prefix.parsesam.temp  $Outdir/$Prefix.sweep.nonrrna
