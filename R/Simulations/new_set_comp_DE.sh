#!/bin/bash
#############################
# les directives Slurm vont ici:

# Your job name (displayed by the queue)
#SBATCH -J DE_m100

# walltime (hh:mm::ss)
#SBATCH -t 2:30:00

# total memory per node
#SBATCH --mem=46GB

# Specify the number of nodes(nodes=) and the number of cores per nodes(tasks-pernode=) to be used
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# use preemptible to speed-up launches
#SBATCH -p preemptible

# change working directory
#SBATCH --chdir=.

##SBATCH --array 1-500

##SBATCH --output=/dev/null
##SBATCH --error=/dev/null


# fin des directives PBS
#############################

# useful informations to print
echo "#############################"
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "SLURM_TASK_ID:" $SLURM_ARRAY_TASK_ID
echo "#############################"

#############################

P1=$1
P2=$2
P3=$3
P4=$4



# on charge R
module load R/4.5.1


# on lance le script exemple_lm.R
Rscript citcdf_gsa/comp_sea_citcdf/new_all_ok/code/comp_methods_DE.R "$P1" "$P2" "$P3" "$P4"
#R CMD BATCH 