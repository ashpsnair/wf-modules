#!/bin/bash

#PBS -l select=1:ngpus=1
#PBS -l walltime=10:00:00
#PBS -P 11003581
#PBS -N fast5-to-pod5
#PBS -j oe

# Change to the directory where the job was submitted 
cd $PBS_O_WORKDIR

#Codes for converting fast5 to pod5
wt_pod5_dir="/home/users/nus/ash.ps/scratch/MCF10A-lab/pod5/WT/"
p3l_pod5_dir="/home/users/nus/ash.ps/scratch/MCF10A-lab/pod5/P3L/"
wt_fast5_dir="/home/users/nus/ash.ps/scratch/MCF10A-lab/WT/MCF10A_WT_fast5_pass"
p3l_fast5_dir="/home/users/nus/ash.ps/scratch/MCF10A-lab/P3L/MCF10A_P3L_2_fast5_pass"

#Convert each fast5 to its relative converted output. 
pod5 convert fast5 $wt_fast5_dir/*.fast5 --output $wt_pod5_dir/ --one-to-one $wt_fast5_dir/
pod5 convert fast5 $p3l_fast5_dir/*.fast5 --output $p3l_pod5_dir/ --one-to-one $p3l_fast5_dir/