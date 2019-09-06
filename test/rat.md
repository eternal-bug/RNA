
+ 质控
```bash
cd ~/data/rat/sequence
mkdir -p output/fastqc

parallel -j 10 "
    fastqc {1} -t 4 -o output/fastqc
" ::: $(ls *.gz )
```

+ 比对

```

```
