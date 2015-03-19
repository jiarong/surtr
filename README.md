<<<<<<< HEAD
surtr
=====

rrna search in shotgun metagenome

Brief instructions:
---------------------
###Dependencies:
- khmer
- screed
- bwa

```bash
git clone git@github.com:jiarong/surtr.git
cd surtr
ls

####copy reference sequence file to db directory
cd db
#refs for all rRNA (SSU, LSU, 5s)
wget http://lyorn.idyll.org/~gjr/public2/misc/surtr_db/rrna.fasta

#refs for SSU rRNA (16s + 18s), *use this one for 16s*
wget http://lyorn.idyll.org/~gjr/public2/misc/surtr_db/ssu.fasta
cd ..

####scripts/surtr.sh is the pipeline. For help, type:
bash scripts/surtr.sh

####run with test data
####it takes several mins for bwa to index refs if ran the first time
bash scripts/surtr.sh db/ssu.fasta test/2k.fa test/2k.test.out/

####to see the output files
ls test/2k.test.out/
```

Existing problems
-----------------
1) false positives from low complexity reads

=======
# surtr
Screen for or out target sequences in metagenome
>>>>>>> b592544344e9bfb60d1a77ee56c87c8f65fe9ecd
