#!/bin/bash
#SBATCH -J itensor10    
##SBATCH --partition=intel
#SBATCH --partition=amd
##SBATCH -N 1
#SBATCH --nodes=1         
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=400GB
##SBATCH --ntasks-per-node=16  
##SBATCH --nodelist=XHPC-6150-20
##SBATCH --exclude=XHPC-9654-[01-03] 
##SBATCH --time=dd-hh:mm:ss    
##SBATCH --output=file_name   
##SBATCH --error=file_name    




#For module 
#module purge
#module load compiler/intel/intel-8458p apps/gsl/gsl-2.8 
#echo $PATH > path.log

#source /public/software/mkl/setvars.sh intel64
#make slurm_cpus_per_node=$SLURM_JOB_CPUS_PER_NODE > Make_out.log
export OMP_NUM_THREADS=10
make clean
make

##echo “SLURM_JOB_PARTITION=$SLURM_JOB_PARTITION”
##echo “SLURM_JOB_NODELIST=$SLURM_JOB_NODELIST”

./t-J 

