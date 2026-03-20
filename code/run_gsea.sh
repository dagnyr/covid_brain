#!/bin/bash
#SBATCH --job-name=gsea
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=gsea.log

source ~/miniconda3/bin/activate covid_brain

pip install gseapy --prefer-binary --quiet

python /scratch/users/dagnyr/brain/code/gsea.py

