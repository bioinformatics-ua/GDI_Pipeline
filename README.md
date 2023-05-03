# GDI_Pipeline


This Snakemake pipeline implements the GATK best-practices workflow for calling small germline variants with possibility of using GATK Hard/VQSR filter or using deepvariant variant caller.


Improvements from pipeline (https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling.git)

1) Update software versions;
2) Use GATK SPARK tools, when possible;
3) Set an human genome version more robust;
4) Add GATK VDSQ filter version;
5) Add deep variant software call available;

## How to install:

```
$ mamba create -c conda-forge -c bioconda --name call_variants snakemake snakedeploy
$ source activate call_variants
$ mkdir -p path/to/project-workdir
$ cd path/to/project-workdir
$ snakedeploy deploy-workflow https://github.com/bioinformatics-ua/GDI_Pipeline.git . --tag main
```

:warning: Important,  If intend to use DeepVariant check if your server/cluster has glibc >= 2.18 

```
$ ldd --version
```

## How to use

```
$ snakemake --jobs 20 --use-conda
```

## using SGE or Slurm

Activate your model

```
$ snakemake --profile sge --jobs 20 --use-conda
$ snakemake --profile slurm --jobs 20 --use-conda
```

