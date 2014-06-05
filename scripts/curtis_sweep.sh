#! /usr/bin/env bash

# stop when fail
set -e
set -o pipefail 

export PYTHONPATH=/home/gjr/data/SSUsearch/compareTools/curtis/lib:$PYTHONPATH
python -c "import khmer; print 'khmer version: {}'.format(khmer.__version__)"

Seq=/home/gjr/data/SSUsearch/compareTools/data/simulated_data/mg_sim_l200_e0_1.fasta
#Seq=/home/gjr/data/SSUsearch/compareTools/data/simulated_data/mg_sim_l200_e2_1.fasta
#Seq=/home/gjr/data/SSUsearch/compareTools/data/simulated_data/mg_sim_l200_e5_1.fasta
#Seq=/home/gjr/data/SSUsearch/compareTools/data/test.5M.fa
#Ref=/home/gjr/software/RNA/ribopicker-standalone-0.4.3/db/ssr108_qi_rep.fa
Ref=/home/gjr/data/SSUsearch/compareTools/curtis/db/rrna.fasta
#Outdir=test_curtis
Outdir=simu_out_rrna
#Prefix=5M
Prefix=$(basename $Seq .fasta).sweep

#module load bwa
Scriptpath=$( cd "$(dirname "$0")" && pwd -P ) 
Binpath=$( cd "$(dirname "$0")" && cd ../bin && pwd -P ) 
Dbpath=$( cd "$(dirname "$0")" && cd ../db && pwd -P ) 

mkdir -p $Outdir
time(
python $Scriptpath/sweep_read.py -k 21 -N 4 -x 1e8 -P 10 $Ref $Seq $Outdir/$Prefix.sweep.nonrrna $Outdir/$Prefix.sweep.rrna
)

