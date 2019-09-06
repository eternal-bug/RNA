
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

~/Applications/biosoft/hisat2-2.1.0/hisat2
```
