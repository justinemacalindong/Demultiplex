#!/bin/bash

#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=R1_1      ### Job Name
#SBATCH --output=results-%j.out         ### File in which to store job output
#SBATCH --error=results-%j.err          ### File in which to store job error messages
#SBATCH --time=0-08:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --account=bgmp      ### Account used for job submission
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --nodelist=n278   

dir="/projects/bgmp/justinem/bioinfo/Bi622/Demultiplex/Assignment-the-first"
sdir="/projects/bgmp/shared/2017_sequencing"

/usr/bin/time -v $dir/second_test.py -f $sdir/1294_S1_L008_R1_001.fastq.gz -pn R1_1.png -r 101 -pt 1294_S1_L008_R1_001.fastq.gz