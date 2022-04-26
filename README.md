# Structural Variant Discovery
## Overview
These pipelines were developed to call large structural variants (SVs) from reads mapped to [repeat-masked](https://www.repeatmasker.org/) reference genomes. The two callers used here are [Delly](https://github.com/dellytools/delly) and [Lumpy](https://github.com/arq5x/lumpy-sv), and the overlap among calls is analysed using the [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) consensus framework.

I've developed the pipeline using [Nextflow](https://www.nextflow.io/) version 20.10.0 (the most recent version on Imperial College London's HPC; date: 26/04/22). To map reads to the repeat-masked reference genome, please see a [previous nextflow workflow](https://github.com/dthorburn/Genomic_Read_Processing). 

## Pipeline framework
**BP**: Split and discordant reads are extracted from indexed bam files.

**LC**: Call SVs using LumpyExpress. 

**LT**: Genotype Lumpy SVs. 

**DC**: Call SVs using Delly. 

**DF**: Remove failed SVs using GATK SelectVariants. 

**SV**: Merge overlapping SVs to create a consensus call. 

## Steps to run this pipeline:
  1.  Download a copy of the appropraite reference genome (i.e., *Anopheles gambiae* genome from [vectorbase](https://vectorbase.org/vectorbase/app/record/dataset/DS_2251b21396#description)). 
  2. Clone this repository. 
  3. Update the project directory path and add required (and optional) arguments to the `Structural_Variants.sh` PBS job submission script. 
  4. Initialise the conda environment using the command ```conda env create --file Structural_Variants.yaml```. You can check the installation by comparing the versions using ```nextflow run Structural_Variants.nf -c Structural_Variants.config --profile imperial --help``` with those installed using ```module load anaconda3/personal; source activate SVs; conda list```. 
  5. Run the pipeline using `qsub Structural_Variants.sh`. 
  
***NB*** The pipeline and example code will only work on Imperial's HPC. If you are not using a PBS job submission system, you'll need to update all the scripts to reflect the job submission system. See [Nextflow doucmentation](https://www.nextflow.io/docs/latest/executor.html) for help. 

The help message from the nextflow script is below:
```
Usage:
  These pipelines were developed to call large structural variants from reads mapped to repeatmaked reference genomes. The two
  callers used here are Delly and Lumpy, and the overlap among calls is analysed using the SURVIVOR consensus framework.

  Pipeline framework:
        LC → LT
     ⭧          ⭨
  BP               SV
     ⭨          ⭧
        DC → DF

  To use, there are 3 steps:
  1. Update project directory path in Structural_Vars.sh
  2. Add required arguments listed below
  3. Submit pipeline coordinator using qsub Structural_Vars.sh

  If you require available HPC jobs for alternative scripts lower job concurrency options.

  Required arguments:
    --RefGen                                        Path to reference fasta. Usage '--RefGen /path/to/genome.fasta'
    --InDir                                         Path to input bam directory. Required even if skipping BP step.

  Optional arguments:
    -w                                              Path to nextflow working directory. (Default: ./work)
    --help                                          Show this message
    --version                                       Show versions used to develop pipeline
    --Chroms                                        User defined chromosome selection (Default: all major LGs in AgamP4).
                                                    Usage '--Chroms "AgamP4_2R,AgamP4_3R"'. Selection must be comma
                                                    delimited in quotes and match the names of the contigs in the
                                                    fasta index file.
    --LC_args                                       Optional arguments for LumpyExpress
    --DC_args                                       Optional arguments for Delly
    --LT_args                                       Optional arguments for svtyper
    --DF_args                                       Optional arguments for GATK SelectVariants

  Concurrency options:                              Imperial HPC only permits 50 jobs per user. These options limit the
                                                    number of concurrent processes running per step. NB. Multiple
                                                    processes can be running at the same time.
    --BP_Forks                                      Default: 10
    --LC_Forks                                      Default: 10
    --DC_Forks                                      Default: 10
    --LT_Forks                                      Default: 10
    --DF_Forks                                      Default: 10

  Debugging arguments:
    -resume                                         Resumes pipeline once errors are resolved. Usage: '-resume curious_borg'
                                                    when log file shows "Launching `GATK_Variant_Call.nf` [curious_borg]"
    --Skip_BP                                       Skips processing bams
    --Skip_LC                                       Skips LumpyExpress Call
    --Skip_DC                                       Skips Delly Call
    --Skip_LT                                       Skips SVTyper
    --Skip_DF                                       Skips GATK SelectVariants
    --BP_threads                                    Number of threads for each subprocess - swap BP for process any acronym to
                                                    alter other processes. (i.e., LC_walltime = 24)
    --BP_memory                                     Number of Gb of memory for each subprocess
    --BP_walltime                                   Number of hours for each subprocess (72 is maximum)
```
