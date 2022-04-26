#!/bin/sh
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=4:mem=12gb
#PBS -N NF_SV_Coordinator
#PBS -j oe

Project_Dir=/rds/general/user/dthorbur/home/ephemeral/00_Work/04_Collarless_Raw/05_SVs/01_AgamP4
cd $Project_Dir

module load nextflow/20.10.0

## To see help message paste the following command (without ##) into the console and press enter:
## module load nextflow/20.10.0; nextflow run GATK_Variant_Call.nf -c GATK_Variant_Call.config --profile imperial --help

echo "Starting: `date`"
nextflow run Structural_Vars.nf -c Structural_Vars.config \
	--profile imperial \
	--RefGen "/rds/general/user/dthorbur/projects/tmstorage/live/acolmop1_ref_genome/GCA_016920705.1_AcolMOP1_genomic.fna" \
	--InDir  "/rds/general/user/dthorbur/home/tmstorage/live/Collarless/01_Reads/02_AcolMOP1" \
	--Chroms "AgamP4_2L,AgamP4_2R,AgamP4_3L,AgamP4_3R,AgamP4_X"
	
echo "Finished: `date`"

