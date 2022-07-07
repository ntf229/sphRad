#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00
#SBATCH --mem=40GB
#SBATCH --job-name=analysis
#SBATCH --mail-type=END
#SBATCH --output=slurm_out/slurm_%x.out
#SBATCH --mail-user=ntf229@nyu.edu

module purge

cd /home/ntf229/containers

singularity exec --overlay overlay-15GB-500K.ext3:ro \
	    /scratch/work/public/singularity/cuda11.0-cudnn8-devel-ubuntu18.04.sif \
/bin/bash -c "source /ext3/env.sh; 
python /home/ntf229/sphRad/bin/analysis.py \
--ageSmooth=True --SF=True --tauClear=2 --clumps=False --numCells=8 --numClumps=34 \
--inc=0 --az=0 --numPhotons=1e9 --pixels=2000 \
--dustFraction=0.1 --SSP=FSPS_Chabrier --maxTemp=16000"



