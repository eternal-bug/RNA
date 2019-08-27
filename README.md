# RNA-seq

## 0. 介绍


```bash
      database                   Workflow                        tools
======================================================================================
  
+=================+     +-------------------------+                               
|     database    |     |      Quality Analysis   |---------------> fastqc 
+=================+     +-------------------------+                                
|+------+         |                 v                                      
|| rRNA |---------|--+  +-------------------------+
|+------+         |  |  | Base Quality Filtering  |------------> TrimGalore
|  +------+       |  |  +-------------------------+
|  |genome|-------|-+|              v
|  +------+       | ||  +-------------------------+
|     +----------+| |+->| rRNA Sequence Filtering |------------> SortMeRNA
|     |  Genome  || |   +-------------------------+
|     |Annotation|| |               v
|     +----------+| |   +-------------------------+
|          |      | +-->|   Genome Alignment      |------------> hisat2
+----------|------+     +-------------------------+
           |                        v
           |            +-------------------------+
           +----------->|  Count Mapped Reads     |------------> HTseq
                        +-------------------------+
                                    v
                        +-------------------------+
                        | Differential Expression |------------> DESeq2
                        +-------------------------+
                                    v
                        +-------------------------+
                        |     Pathway analysis    |------------> ClusterProfiler
                        +-------------------------+
```

## 1. 前期准备

在进行数据处理之前，需要将大致的目录生成好，这样一来将数据存放的位置有序一些，便于查看，二来在程序运行的时候能够更加方便的明确填写文件路径。将后面可能会用到的文件存放的文件夹大致拟定一下，便于自己理解和查找文件，这里不需要将后面所有可能用到的文件夹都新建，只是建立一个整体的文件夹的目录结构就行了。这样层级清晰，可以防止出错，在某种程序上提高了效率节省时间。

```bash
# 首先定位到自己这个用户的位置
$ cd ~

# 然后新建项目的文件夹
# 新建一个biosoft文件夹用于存放生物软件工具
$ mkdir biosoft

# 新建一个项目文件夹其中包含大鼠的文件夹
$ mkdir -p project/rat

# 进入
$ cd project/rat

# 将今后各个文件需要存放的文件夹可以事先拟定一下
$ mkdir annotation genome sequence output script
```

| 文件夹名 | 说明 |
| ---| --- |
| annotation | 存放注释文件(.gff .bed .gff3) |
| genome | 存放基因组与索引文件(.fa .bt)|
| sequence | 存放测序数据(.fastq.gz) |
| output | 存放各种处理的输出文件 |
| script | 存放脚本的位置 |

使用`tree`命令看一下我们设置的目录结构

```bash
# 首先进入大鼠的项目文件夹中
$ cd ~/project/rat

$ tree

.
├── annotation  用于存放大鼠的基因组注释信息(.gff/gtf)
├── genome      用于存放大鼠的基因组数据(.fasta)
├── output      用于存放处理后的数据
└── sequence    用于存放测序原始数据
```

## 2. 工具下载

### 2.0 conda

管理生信工具

```bash
# 下载文件
wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh

# 安装
bash miniconda.sh

# 配置
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
```

建立python3.6的环境 

```bash

```

### 2.1 sratoolkit

```bash
$ cd ~/biosoft
# mac
$ wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6-1/sratoolkit.2.9.6-1-mac64.tar.gz
$ tar -xzvf sratoolkit.2.9.6-1-mac64.tar.gz
$ mv sratoolkit.2.9.6-1-mac64 sratoolkit.2.9.6
$ cd sratoolkit.2.9.6/bin

# 导入临时环境变量
$ export PATH="$(pwd):$PATH"
$ prefetch --help
```

### 2.2 fastqc

对测序文件质量控制

|     | 站点 |
| --- | --- |
| 官网 | http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| 手册 | http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/ |
| 中文解释 | https://www.plob.org/article/5987.html |


```bash
cd ~/biosoft
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
unzip fastqc_v0.11.7.zip
cd FastQC

# 导入临时环境变量
export PATH="$(pwd):$PATH"

# 测试是否能运行
fastqc --help
```
### 2.3 multiqc

将fastqc的统计结果汇聚成一个网页可视化文件，便于查看

|      | 站点 |
| --- | --- |
| 官网 |  https://multiqc.info/ |
| 文档 | https://multiqc.info/docs/ |
| 中文解读 | http://fbb84b26.wiz03.com/share/s/3XK4IC0cm4CL22pU-r1HPcQQ1iRTvV2GwkwL2AaxYi2fXHP7 |

```bash
# 使用python的安装器安装
pip install multiqc
```
### 2.4 cutadapt

用于去除测序接头

|      | 站点 |
| --- | --- |
| 手册 | https://cutadapt.readthedocs.io/en/stable/guide.html |

```bash
pip install cutadapt
```
### 2.5 trimmomatic

trimmomatic是一款多线程命令行工具，可以用来修剪Illumina (FASTQ)数据以及删除接头，是目前使用最多的高通量测序数据清洗的工具。

|   | 站点 |
| --- | --- |
| 官网 |  http://www.usadellab.org/cms/index.php?page=trimmomatic |
| 手册 | http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf |
| 中文解读 |  |


```bash
cd ~/biosoft
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
unzip Trimmomatic-0.38.zip

cd Trimmomatic-0.38

# 导入临时环境变量
export PATH="$(pwd):$PATH"
```
### 2.6 hisat2

作为bowtie2和tophat的继任者，它在RNA-seq中使用较多。

|        | 站点 |
| --- | --- |
| 官网 | https://ccb.jhu.edu/software/hisat2/index.shtml |
| 手册 | https://ccb.jhu.edu/software/hisat2/manual.shtml |

+ 首先用浏览器进入`hisat2`[官网](https://ccb.jhu.edu/software/hisat2/index.shtml)

+ 在右侧有对应的不同系统的安装版本，根据自己的系统进行选择

| version 2.1.0 | 6/8/2017 |
| --- | --- |
| [   Source code](http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-source.zip) |
| [   Linux x86_64 binary](http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-Linux_x86_64.zip) |
| [   Mac OS X x86_64 binary](http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-OSX_x86_64.zip) |
| [   Windows binary](http://www.di.fc.ul.pt/~afalcao/hisat2_windows.html) |

另外是参考文献

Kim D, Langmead B and Salzberg SL. [**HISAT: a fast spliced aligner with low memory requirements**](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3317.html). _[Nature Methods](http://www.nature.com/nmeth)_2015


+ 这里我们选择`Mac OS X`的，对着[Mac OS X x86_64 binary](http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-OSX_x86_64.zip)右键，复制链接地址目

+ 回到终端上来：

```bash
# 下载
$ cd ~/biosoft/
$ wget http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-OSX_x86_64.zip

# 解压缩
$ unzip hisat2-2.1.0-OSX_x86_64.zip
$ cd hisat2-2.1.0

# 导入临时环境变量
$ export PATH="~/biosoft/hisat2-2.1.0:$PATH"

# 测试是否可用
$ hisat2 -h
```
### 2.7 samtools

比对得到的sam或者bam文件的各种操作

|   |  站点 |
| --- | --- |
| 官网 | http://www.htslib.org/ |
| 手册 | http://quinlanlab.org/tutorials/samtools/samtools.html |
| 中文解读 | https://www.jianshu.com/p/6b7a442d293f |

samtools是对比对后得到的文件进行格式转化处理合并等操作的工具。

```bash
cd ~/biosoft
wget -c https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar jxvf samtools-1.9.tar.bz2
cd samtools-1.9
./configure --prefix=$(pwd)
make

# 导入临时环境变量
export PATH="$(pwd):$PATH"
```
### 2.8 HTseq

对比对后的文件进行read计数

```
pip install HTseq
```

### 2.9 R

- 官网：https://www.r-project.org

R语言中集合了多种生物信息学的分析工具，其中RNA-seq分析的工具更是丰富，R语言最擅长统计学分析，这里在后续的基因表达量分析，差异分析以及作图时候会用到

- 下载

点击左上角`CRAN`，往下拉找到`china`的站点，就选第一个[清华的站点](https://mirrors.tuna.tsinghua.edu.cn/CRAN/)，点击进去，之后最上面有三个对应系统的安装包，按照自己的系统下载，这里点击[Mac OS](https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/macosx/)，点击[R-3.6.1.pkg](https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/macosx/R-3.6.1.pkg)就开始下载了，下载之后双击安装包安装。

### 2.10 Rstudio

因为R自身带的界面使用起来不是特别方便和美观，这里使用Rstudio来对R的显示效果进行提升，除此之外还有其他功能。

进入网站：https://www.rstudio.com/products/rstudio/download/

点击`RStudio Desktop Open Source License` 下面的`DOWNLOAD`，之后选择对应的版本的，这里选择`MAC OS`版本[RStudio 1.2.1335 - Mac OS X 10.12+ (64-bit)](https://download1.rstudio.org/desktop/macos/RStudio-1.2.1335.dmg) 下载完成之后双击安装。

**注意**：安装完成`R`之后再安装`Rstudio`

### 2.11 parallel

parallel是进行多线程运行的工具，并行运行可以提升效率，节省时间

```
brew install parallel
```

### === StringTie - Ballgown流程还未走通，

### StringTie

能够应用流神经网络算法和可选的de novo组装进行转录本组装并预计表达水平。与Cufflinks等程序相比，StringTie实现了更完整、更准确的基因重建，并更好地预测了表达水平。

|  | 站点 |
| --- | --- |
| 官网 | http://ccb.jhu.edu/software/stringtie/ |
| 手册 | http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual |
| 中文解读 | https://www.plob.org/article/12865.html |

+ 安装

```bash
cd ~/biosoft

wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.6.OSX_x86_64.tar.gz

tar -xzvf stringtie-1.3.6.OSX_x86_64.tar.gz
mv stringtie-1.3.6.OSX_x86_64 stringtie-1.3.6
cd stringtie-1.3.6

export PATH="$(pwd):$PATH"

stringtie --help
```
### Ballgown

是R语言中基因差异表达分析的工具，能利用RNA-Seq实验的数据(StringTie, RSEM, Cufflinks)的结果预测基因、转录本的差异表达。

+ 安装



