
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
