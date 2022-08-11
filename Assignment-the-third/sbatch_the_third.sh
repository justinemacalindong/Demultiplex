#!/bin/bash

#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=demux!     ### Job Name
#SBATCH --output=results-%j.out         ### File in which to store job output
#SBATCH --error=results-%j.err          ### File in which to store job error messages
#SBATCH --time=0-08:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --account=bgmp      ### Account used for job submission
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1


dir="/projects/bgmp/justinem/bioinfo/Bi622/Demultiplex/Assignment-the-third"
sdir="/projects/bgmp/shared/2017_sequencing"

/usr/bin/time -v $dir/demultiplex.py -r1 $sdir/1294_S1_L008_R1_001.fastq.gz -r2 $sdir/1294_S1_L008_R4_001.fastq.gz -i1 $sdir/1294_S1_L008_R2_001.fastq.gz -i2 $sdir/1294_S1_L008_R3_001.fastq.gz -k $sdir/indexes.txt
/usr/bin/time -v pigz $dir/read*