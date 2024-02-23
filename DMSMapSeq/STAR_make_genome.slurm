#!/bin/bash
#SBATCH --partition=ycga
#SBATCH --job-name=make_genome
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=100GB
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2

module load STAR


STAR --runThreadN 4 --runMode genomeGenerate --genomeDir DMS_genome --genomeFastaFiles short_long_MS2.fa --genomeSAindexNbases 9 --limitGenomeGenerateRAM 37000000000