surtr
=====

rrna search in shotgun metagenome

Brief instructions:
---------------------
git clone git@github.com:jiarong/surtr.git
cd surtr
ls

#copy reference sequence file to db directory
cd db
wget http://lyorn.idyll.org/~gjr/public2/misc/surtr_db/rrna.fasta
cd ..

#scripts/surtr.sh is the pipeline. For help, type:
bash scripts/surtr.sh

#run with test data
bash scripts/surtr.sh db/rrna.fasta test/2k.fa test/2k.test.out/

#to see the output files
ls test/2k.test.out/

Existing problems
-----------------
1) false positives from low complexity reads

