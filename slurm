#!/bin/bash
#SBATCH -J GaAs_gpu
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --gres=gpu:NVIDIAA800-SXM4-80GB:1
#SBATCH -o %j.log
#SBATCH -e %j.err
echo "SLURM_JOB_PARTITION=$SLURM_JOB_PARTITION"
echo "SLURM_JOB_NODELIST=$SLURM_JOB_NODELIST"
echo Time is `date`
echo Directory is $PWD
echo This job runs on the following nodes:
echo $SLURM_JOB_NODELIST
echo This job has allocated $SLURM_JOB_CPUS_PER_NODE cpu cores.

ulimit -s unlimited
module purge
module load apps/matlab/R2024a
matlab -nodisplay -nodesktop -nosplash -r test_GaAs_fullfreq > log 

YY=`date +%Y`
MM=`date +%m`
DD=`date +%d`
echo '  '$PWD'  Finished at  ' `date` '  ' >> /public/home/mxy/slurm_log/$YY-$MM-$DD.log
