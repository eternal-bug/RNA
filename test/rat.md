
+ 质控
```bash
cd ~/data/rat/sequence
mkdir -p output/fastqc

fastqc *.gz -t 20 -o output/fastqc

```

+ 比对

```
cd ~/database/genome/rn6/

~/Applications/biosoft/hisat2-2.1.0/hisat2-build  -p 18 rn6.fa rn6



parallel -k -j 6 "
    ~/Applications/biosoft/hisat2-2.1.0/hisat2 \
      -t -x ~/database/genome/rn6/rn6 \
      -1 {1}_1.fq.gz -2 {1}_2.fq.gz -S ../output/align/{1}.sam \
      2>../output/align/log/{1}.log
" ::: $(ls *.gz | perl -n -e 'print $1."\n" if m/(\w+?)_/' | uniq)



parallel -k -j 6 "
    ~/Applications/biosoft/hisat2-2.1.0/hisat2 \
      -t -x ~/database/genome/rn6/rn6 \
      -1 {1}_1.fq.gz -2 {1}_2.fq.gz -S ../output/align/{1}.sam \
      2>../output/align/log/{1}.log
" ::: LHA{1..3} LHB{1..3}
```

| sample | ratio | time |
| ---- | ----- | ----- |
| LHA1 | 96.09 | 106.10 |
| LHA2 | 96.20 | 101.70 |
| LHA3 | 96.34 | 86.15 |
| LHB1 | 96.29 | 80.02 |
| LHB2 | 96.31 | 87.13 |
| LHB3 | 0.00  | 0.00 |
| LHC1 | 96.21 | 59.00 |
| LHC2 | 96.10 | 56.18 |
| LHC3 | 96.18 | 57.20 |
| LMA1 | 96.75 | 64.57 |
| LMA2 | 96.68 | 71.27 |
| LMA3 | 96.74 | 87.63 |
| LMB1 | 96.53 | 107.83 |
| LMB2 | 96.45 | 88.12 |
| LMB3 | 96.46 | 82.23 |
| LMC1 | 96.57 | 78.52 |
| LMC2 | 96.39 | 91.32 |
| LMC3 | 0.00  | 0.00 |
| NH   | 96.35 | 92.30 |
| NM   | 96.41 | 87.52 |
