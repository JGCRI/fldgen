#!/bin/zsh

#SBATCH -p short
#SBATCH -t 180
#SBATCH -A IHESD


module purge
module load gcc/8.1.0
module load netcdf
module load R/3.4.3

## 
echo Rscript -e \"source('train-emulators.R'); train_models('$1')\"

Rscript -e "source('train-emulators.R'); train_models('$1')"
