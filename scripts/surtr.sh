#! /usr/bin/env bash

# stop when fail

set -e
set -o pipefail 

module load bwa
module load screed
module load khmer/1.0.1-rc2
#export PYTHONPATH=/mnt/home/guojiaro/Documents/lib/git/khmer/python
Scriptpath=/mnt/home/guojiaro/Documents/software/RNA/surtr/scripts


if [ $# -ne 3 ]
then
  echo "Usage: $(basename $0) <ref> <seqfile> <outdir>"
  exit 1
fi

show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-c coverage] [-i identity_alignment] [-I identity_read] ref_file read_file out_dir
    
    -h          display this help and exit
    -f OUTFILE  write the result to OUTFILE instead of standard output.
    -v          verbose mode. Can be used multiple times for increased
                verbosity.
EOF
}                

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

while getopts "h:c:i:I:" opt;
do
  case "$opt" in
    h)
      show_help
      exit 0
      ;;
    c)
      coverage=$OPTARG
      ;;
    i)
      iden_ali=$OPTARG
      ;;
    I)
      iden_read=$OPTARG
      ;;
    '?')
      show_help >&2
      exit 1
      ;;
  esac
done

shift "$((OPTIND-1))" # Shift off the options and optional --.
echo "coverage=$coverage, iden_ali=$iden_ali, iden_read=$iden_read" 

Ref=$1
Seq=$2
Outdir=$3

#Prefix=5M
Prefix=$(basename $Seq .fasta).surtr
Log="$Outdir"/"$Prefix".log



###Scriptpath=$( cd "$(dirname "$0")" && pwd -P ) 
###Binpath=$( cd "$(dirname "$0")" && cd ../bin && pwd -P ) 
###Dbpath=$( cd "$(dirname "$0")" && cd ../db && pwd -P ) 
###Libpath=$( cd "$(dirname "$0")" && cd ../lib && pwd -P ) 
###Ref=$Dbpath/rrna.fasta


python -c "import khmer; print 'khmer version: {}'.format(khmer.__version__)"

if [ ! -f $Ref.bwt ];
then
  bwa index $Ref
fi

mkdir -p $Outdir
date
echo 'start surtr'

python $Scriptpath/filter_low_comp.py $Seq - | \
python $Scriptpath/sweep_read.py -k 21 -N 4 -x 1e8 -P 1 $Ref $Seq $Outdir/$Prefix.sweep.nonrrna - | \
bwa mem -t 1 -c 10000000 $Ref - | \
python $Scriptpath/parse_sam.py -o $Outdir -p $Prefix.parsesam -I 50 -

#python $Scriptpath/parse_sam.py -o $Outdir -p $Prefix.parsesam -c 75 -i 50 -

echo 'surtr finished'

cat $Outdir/$Prefix.sweep.nonrrna $Outdir/$Prefix.parsesam.nonrrna > $Outdir/"$Prefix"_nonrrna.fasta

cat $Outdir/$Prefix.parsesam.rrna $Outdir/$Prefix.parsesam.nonrrna > $Outdir/"$Prefix"_sweep.rrna

mv $Outdir/$Prefix.sweep.nonrrna $Outdir/"$Prefix"_sweep.nonrrna
mv $Outdir/$Prefix.parsesam.rrna  $Outdir/"$Prefix"_rrna.fasta
mv $Outdir/$Prefix.parsesam.csv  $Outdir/"$Prefix"_rrna.csv

mkdir -p $Outdir/temp_files
#rm -f $Outdir/$Prefix.parsesam.rrna $Outdir/$Prefix.parsesam.nonrrna $Outdir/$Prefix.parsesam.temp  $Outdir/$Prefix.sweep.nonrrna
#rm -f $Outdir/"$Prefix"_sweep.nonrrna $Outdir/"$Prefix"_sweep.rrna

mv  $Outdir/$Prefix.parsesam.nonrrna $Outdir/$Prefix.parsesam.temp  $Outdir/temp_files
mv $Outdir/"$Prefix"_sweep.nonrrna $Outdir/"$Prefix"_sweep.rrna $Outdir/temp_files
