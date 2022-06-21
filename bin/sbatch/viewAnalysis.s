#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00
#SBATCH --mem=10GB
#SBATCH --job-name=viewAnalysis
#SBATCH --mail-type=END
#SBATCH --output=slurm_out/slurm_%x.out
#SBATCH --mail-user=ntf229@nyu.edu

module purge

cd /home/ntf229/containers

singularity exec --overlay overlay-15GB-500K.ext3:ro \
	    /scratch/work/public/singularity/cuda11.0-cudnn8-devel-ubuntu18.04.sif \
/bin/bash -c "source /ext3/env.sh; 
python /home/ntf229/sphRad/bin/viewAnalysis.py \
--ageSmooth=True --SF=False --tauClear=2 --clumps=False --numCells=8 --numClumps=34 \
--dust=False --numPhotons=1e8 --pixels=500 \
--dustFraction=0.1 --SSP=FSPS_Chabrier --maxTemp=16000"



