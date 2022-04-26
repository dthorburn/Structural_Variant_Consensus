#!/usr/bin/env nextflow

/*
 * Pipeline developed for consensus (SURVIVOR) structural variant discovery using Lumpy and Delly.
 *
 * Author: Miles Thorburn <d.thorburn@imperial.ac.uk>
 *
 * Date last modified: 26/04/2022
 */

def helpMessage() {
  log.info """
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
"""
}

def versionMessage() {
  log.info """
==============================================================================================================================
_libgcc_mutex             0.1                        main
_openmp_mutex             4.5                       1_gnu
bamkit                    16.07.26                   py_0    bioconda
bcftools                  1.9                  h68d8f2e_9    bioconda
blas                      1.0                         mkl
boost-cpp                 1.70.0               ha2d47e9_1    conda-forge
bzip2                     1.0.8                h7b6447c_0
ca-certificates           2022.3.18            h06a4308_0
certifi                   2020.6.20          pyhd3eb1b0_3
curl                      7.71.1               h8f29fe8_2
cytoolz                   0.10.1           py27h7b6447c_0
delly                     0.8.3                h43566fd_0    bioconda
gawk                      5.1.0                h7b6447c_0
gsl                       2.5                  h294904e_1    conda-forge
htslib                    1.9                  h4da6232_3    bioconda
icu                       58.2                 he6710b0_3
intel-openmp              2022.0.1          h06a4308_3633
krb5                      1.19.1               h3535a68_0
ld_impl_linux-64          2.35.1               h7274673_9
ldc                       1.20.0               h9a1ace1_1    conda-forge
libblas                   3.9.0           1_h6e990d7_netlib    conda-forge
libcblas                  3.9.0           3_h893e4fe_netlib    conda-forge
libcurl                   7.71.1               h303737a_2
libdeflate                1.6                  h516909a_0    conda-forge
libedit                   3.1.20191231         h46ee950_2    conda-forge
libffi                    3.2.1             hf484d3e_1007
libgcc                    7.2.0                h69d50b8_2
libgcc-ng                 9.3.0               h5101ec6_17
libgfortran-ng            7.5.0               ha8ba4b0_17
libgfortran4              7.5.0               ha8ba4b0_17
libgomp                   9.3.0               h5101ec6_17
libssh2                   1.9.0                h1ba5d50_1
libstdcxx-ng              9.3.0               hd4cf53a_17
lumpy-sv                  0.2.14a              hdfb72b2_2    bioconda
mkl                       2020.2                      256
mkl-service               2.3.0            py27he904b0f_0
mkl_fft                   1.0.15           py27ha843d7b_0
mkl_random                1.1.0            py27hd6b4f25_0
ncurses                   6.1                  he6710b0_1
numpy                     1.16.6           py27hbc911f0_0
numpy-base                1.16.6           py27hde5b4d6_0
openssl                   1.1.1n               h7f8727e_0
perl                      5.26.2               h14c3975_0
pip                       19.3.1                   py27_0
pysam                     0.16.0.1         py27ha863e18_1    bioconda
python                    2.7.15          h5a48372_1011_cpython    conda-forge
python_abi                2.7                    1_cp27mu    conda-forge
readline                  8.0                  h7b6447c_0
samblaster                0.1.26               h7d875b9_1    bioconda
samtools                  1.9                 h10a08f8_12    bioconda
scipy                     1.2.1            py27h7c811a0_0
setuptools                44.0.0                   py27_0
six                       1.16.0             pyhd3eb1b0_1
sqlite                    3.38.0               hc218d9a_0
survivor                  1.0.7                h9a82719_1    bioconda
svtyper                   0.7.1                      py_0    bioconda
tk                        8.6.11               h1ccaba5_0
toolz                     0.10.0             pyhd3eb1b0_0
util-linux                2.35             py27hfbfebf4_1    conda-forge
wheel                     0.37.1             pyhd3eb1b0_0
xz                        5.2.5                h7b6447c_0
zlib                      1.2.11               h7f8727e_4
==============================================================================================================================
"""
}

log.info """
==============================================================================================================================
                                        Structual Variant Discovery Pipeline v1
==============================================================================================================================

Reference       : ${params.RefGen}
Input           : ${params.InDir}
Processed Bams  : ${PWD}/03_Processed_Bams/
Raw Lumpy VCFs  : ${PWD}/04_Lumpy_Raw/
Raw Delly VCFs  : ${PWD}/05_Delly_Raw/
Lumpy VCFs      : ${PWD}/06_Lumpy_Typed/
Delly VCFs      : ${PWD}/07_Delly_Filtered/
Merged VCFs     : ${PWD}/08_Surviror_Output/

==============================================================================================================================
"""

if (params.help) {
    helpMessage()
    exit 0
}

if (params.version) {
    versionMessage()
    exit 0
}

// Adding error messages for missing required parameters
if(!params.RefGen) {
  log.info"""
ERROR: No reference genome path provided! --RefGen /path/to/genome.fasta
==============================================================================================================================
  """
  helpMessage()
  exit 0
}
if(!params.InDir) {
  log.info"""
ERROR: No input bam directory path provided! --InDir /path/to/bams/
==============================================================================================================================
  """
  helpMessage()
  exit 0
}
                                                            // =========================================================
                                                            // Setting the value channels (can be read unlimited times)
                                                            // =========================================================
ref_genome = file( params.RefGen, checkIfExists: true )
ref_dir    = ref_genome.getParent()
ref_name   = ref_genome.getBaseName()
ref_index  = file( "${ref_dir}/${ref_name}*.fai", checkIfExists: true ).first()
split_sct  = file( params.BWA_split_script, checkIfExists: true )
                                                            // =========================================================
                                                            // Step 1: Processing BAM files
                                                            // =========================================================
if( params.Skip_BP == false ) {
  Channel
    .fromFilePairs("${params.InDir}/*{bam,bai}") { file -> file.name.replaceAll(/.bam|.bai$/,'') }
    .ifEmpty { error "No bams found in ${params.InDir}" }
    .map { ID, files -> tuple(ID, files[0], files[1]) }
    .set { raw_bams }

  // Setting up the chromosome channel
  if( params.Chroms == "" ){
    // Defaulting to using all major An. coluzzii chromosomes
    chromosomes_ch = Channel
                        .from("AgamP4_2L", "AgamP4_2R", "AgamP4_3L", "AgamP4_3R", "AgamP4_X", "AgamP4_Y_unplaced", "AgamP4_UNKN") //AgamP4
//                        .from("CM029350.1", "CM029348.1", "CM029349.1") //AcolMOP1
//                        .from("CM029352.1","CM029353.1","CM029354.1") //AaraD3
    println "No chromosomes specified, using all major chromosomes: AgamP4_2L, AgamP4_2R, AgamP4_3L, AgamP4_3R, AgamP4_X, AgamP4_Y_unplaced, AgamP4_UNKN"
  } else {
    // User option to choose which chromosome will be used
    // This worked with the following syntax nextflow run testing.nf --profile imperial --Chroms "AgamP4_3R,AgamP4_2L"
    chrs = params.Chroms.split(",")
    chromosomes_ch = Channel
                      .from( chrs )
    println "User defined chromosomes set: ${params.Chroms}"
  }

  process Process_Bams {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
    maxRetries 3
    maxForks params.BP_Forks

    tag { SampleID+"-"+chrom }

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.BP_threads}:mem=${params.BP_memory}gb -lwalltime=${params.BP_walltime}:00:00"

    publishDir(
      path: "${params.BP_Dir}",
      mode: 'copy',
    )

    input:
    each chrom from chromosomes_ch
    set SampleID, path(bam), path(bai) from raw_bams
    path split_sct

    output:
    tuple val(SampleID), val(chrom), path("${SampleID}.all.${chrom}.bam"), path("${SampleID}.discordant.${chrom}.bam"), path("${SampleID}.splitters.${chrom}.bam") into processed_bams_lumpy, processed_bams_delly

    beforeScript 'module load python/3.5.4; module load samtools/1.2'

    script:
    """
    ## Discordant Reads
    samtools view -b -F 1294 ${bam} ${chrom} \\
      | samtools sort > ${SampleID}.discordant.${chrom}.bam

    ## Split Reads
    samtools view -h ${bam} ${chrom} \\
      | python ${split_sct} -i stdin \\
      | samtools view -Sb - \\
      | samtools sort > ${SampleID}.splitters.${chrom}.bam

    ## Chrom Bam
    samtools view -h ${bam} ${chrom} \\
      | samtools sort > ${SampleID}.all.${chrom}.bam
    """
  }
}
                                                            // =========================================================
                                                            // Step 2a: LUMPY Call
                                                            // =========================================================
if( params.Skip_LC == false ){
  if( params.Skip_BP ){
    // Both of these solutions should have the same outcome, but let's see if it works with the simpler version first.
    Channel
      .fromFilePairs("${params.BP_Dir}/*{all,discordant,splitters}*.bam")
      .ifEmpty { error "no bams found in ${params.BP_Dir}" }
      .map { files -> tuple(files[0].baseName.replaceAll(".all.*", ""), files[0].baseName.replaceAll("*.all.", "").replaceAll(".bam", ""), files[0], files[1], files[2]) }
      .groupTuple( by: [0,1] )
      set { processed_bams_lumpy }

//    Channel
//      .fromPath("${params.BP_Dir}/*all*.bam")
//      .ifEmpty { error "no full bams found in ${params.BP_Dir}" }
//      .map { file -> tuple(file.baseName.replaceAll(".all.*", ""), file.baseName.replaceAll("*.all.", "").replaceAll(".bam", ""), file) }
//      .set { all_bams }

//    Channel
//      .fromPath("${params.BP_Dir}/*discordant*.bam")
//      .ifEmpty { error "no discordant bams found in ${params.BP_Dir}" }
//      .map { file -> tuple(file.baseName.replaceAll(".discordant.*", ""), file.baseName.replaceAll("*.discordant.", "").replaceAll(".bam", ""), file) }
//      .set { dis_bams }

//    Channel
//      .fromPath("${params.BP_Dir}/*splitters*.bam")
//      .ifEmpty { error "no splitters bams found in ${params.BP_Dir}" }
//      .map { file -> tuple(file.baseName.replaceAll(".splitters.*", ""), file.baseName.replaceAll("*.splitters.", "").replaceAll(".bam", ""), file) }
//      .set { spl_bams }

//    all_bams
//      .combine( dis_bams, by = [0,1] )
//      .combine( spl_bams, by = [0,1] )
//      .set { processed_bams_lumpy }
}
  process Lumpy_Call {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
    maxRetries 3
    maxForks params.LC_Forks

    tag { SampleID+"-"+chrom }

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.LC_threads}:mem=${params.LC_memory}gb -lwalltime=${params.LC_walltime}:00:00"

    publishDir(
      path: "${params.LC_Dir}",
      mode: 'copy',
    )

    input:
    set SampleID, chrom, all_bam, dis_bam, spl_bam from processed_bams_lumpy
    path ref_genome
    path ref_index

    output:
    tuple SampleID, chrom, path("${SampleID}_${chrom}_lumpy.bam.vcf") into raw_lumpy_vcfs

    beforeScript 'module load anaconda3/personal; source activate SVs'

    script:
    """
    lumpyexpress -k -P ${params.LC_args} \\
      -R ${ref_genome} \\
      -B ${all_bam} \
      -S ${spl_bam} \\
      -D ${dis_bam} \\
      -o ${SampleID}_${chrom}_lumpy.bam.vcf
    """
  }
}
                                                            // =========================================================
                                                            // Step 2b: DELLY Call
                                                            // =========================================================
if( params.Skip_DC == false ) {
  if( params.Skip_BP ){
    // See solution to Step 2a if this doesn't work.
    Channel
      .fromFilePairs("${params.BP_Dir}/*{all,discordant,splitters}*.bam")
      .ifEmpty { error "no bams found in ${params.BP_Dir}" }
      .map { files -> tuple(files[0].baseName.replaceAll(".all.*", ""), files[0].baseName.replaceAll("*.all.", "").replaceAll(".bam", ""), files[0], files[1], files[2]) }
      .groupTuple( by: [0,1] )
      set { processed_bams_lumpy }
  }
  process Delly_Call {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
    maxRetries 3
    maxForks params.DC_Forks

    tag { SampleID+"-"+chrom }

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.DC_threads}:mem=${params.DC_memory}gb -lwalltime=${params.DC_walltime}:00:00"

    publishDir(
      path: "${params.DC_Dir}",
      mode: 'copy',
    )

    input:
    set SampleID, chrom, all_bam, dis_bam, spl_bam from processed_bams_delly
    path ref_genome
    path ref_index

    output:
    tuple SampleID, chrom, path("${SampleID}_${chrom}_delly.bam.vcf") into raw_delly_vcfs

    beforeScript 'module load anaconda3/personal; source activate SVs'

    script:
    """
    delly call ${params.DC_args} \\
      -g ${ref_genome} \\
      -o ${SampleID}_${chrom}_delly.bam.bcf \\
      ${all_bam}
    bcftools view ${SampleID}_${chrom}_delly.bam.bcf > ${SampleID}_${chrom}_delly.bam.vcf
    """
  }
}
                                                            // =========================================================
                                                            // Step 3a: LUMPY Typing
                                                            // =========================================================
if( params.Skip_LT == false ) {
  if( params.Skip_LC ){
    Channel
      .fromPath("${params.LC_Dir}/*_lumpy.bam.vcf")
      .map { file -> tuple(file.baseName.replaceAll("_*", ""), file.baseName.replaceAll("_lumpy.bam.vcf", "").replaceAll("*_", ""), file) }
      .set { raw_lumpy_vcfs }
  }

  // This step requires the original bam too, so a new tuple needs to be made by joining.                                                CHECK THIS WORKS!
  Channel
    .fromFilePairs("${params.InDir}/*{bam,bai}") { file -> file.name.replaceAll(/.bam|.bai$/,'') }
    .ifEmpty { error "No bams found in ${params.InDir}" }
    .map { ID, files -> tuple(ID, files[0], files[1]) }
    .set { og_bams }

  raw_lumpy_vcfs
    .combine( og_bams, by: 0 )
    .set { raw_lumpy_vcfs_bams }

  process Lumpy_Type {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
    maxRetries 3
    maxForks params.LT_Forks

    tag { SampleID+"-"+chrom }

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.LT_threads}:mem=${params.LT_memory}gb -lwalltime=${params.LT_walltime}:00:00"

    publishDir(
      path: "${params.LT_Dir}",
      mode: 'copy',
    )

    input:
    set SampleID, chrom, vcf, bam from raw_lumpy_vcfs_bams

    output:
    tuple SampleID, chrom, path("${SampleID}_${chrom}_lumpy.typed.NoMissing.vcf") into lumpy_out

    beforeScript 'module load anaconda3/personal; source activate SVs'

    script:
    """
    svtyper \
      --verbose \
      -B ${bam} \\
      -i ${vcf} ${params.LT_args} \\
      > ${SampleID}_${chrom}_lumpy.typed.vcf

    grep -vw "0/0" ${SampleID}_${chrom}_lumpy.typed.vcf > ${SampleID}_${chrom}_lumpy.typed.NoMissing.vcf
    """
  }
}
                                                            // =========================================================
                                                            // Step 3b: DELLY Filter LowQuals
                                                            // =========================================================
if( params.Skip_DF == false ) {
  if( params.Skip_DC ) {
    Channel
      .fromPath("${params.DC_Dir}/*_delly.bam.vcf")
      .ifEmpty { error "No VCFs found in ${params.DC_Dir}" }
      .map { file -> tuple(file.baseName, file.baseName.replaceAll("_delly.bam.vcf", "").replaceAll("*_", ""), file) }
      .set { raw_delly_vcfs }
  }

  process Delly_Filter {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
    maxRetries 3
    maxForks params.DF_Forks

    tag { SampleID+"-"+chrom }

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.DF_threads}:mem=${params.DF_memory}gb -lwalltime=${params.DF_walltime}:00:00"

    publishDir(
      path: "${params.DF_Dir}",
      mode: 'copy',
    )

    input:
    set SampleID, chrom, vcf from raw_delly_vcfs
    path ref_genome
    path ref_index

    output:
    tuple SampleID, chrom, path("${SampleID}_${chrom}_delly.filtered.NoMissing.vcf") into delly_out

    beforeScript 'module load gatk/4.0'

    script:
    """
    gatk SelectVariants \\
      -R ${ref_genome} \\
      -V ${vcf} \\
      -O ${SampleID}_${chrom}_delly.filtered.vcf \\
      -select 'vc.isNotFiltered()' ${params.DF_args}

    grep -vw "0/0" ${SampleID}_${chrom}_delly.filtered.vcf > ${SampleID}_${chrom}_delly.filtered.NoMissing.vcf
    """
  }
}
                                                            // =========================================================
                                                            // Step 4: SURVIVOR Merging SV Calls
                                                            // =========================================================
// TBC - Testing the first part of this scripts first.
