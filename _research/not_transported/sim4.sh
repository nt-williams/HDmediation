#!/bin/sh
#$ -l mem=10G,time=2:00:00
cd HDmediation
Rscript=/nfs/apps/R/4.0.3/bin/Rscript
export R_LIBS_USER=/ifs/home/msph/epi/ntw2117/R_4.0
${Rscript} _research/not_transported/sim4.R
