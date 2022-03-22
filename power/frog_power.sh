#!/bin/bash
# --------------------------------
#$ -S /bin/bash
# frog_power
#$ -N frog_power
#
# output file
#$ -o $HOME/projects/frog_calls/power/logs/log_o.$TASK_ID.out
# error file
#$ -e $HOME/projects/frog_calls/power/logs/log_e.$TASK_ID.out
#$ -M timothy.farkas@gmail.com
#$ -m aes 
cd $HOME/projects/frog_calls/power

#
#$ -t 1-56
#$ -cwd

# load R
module load R/3.1.1

ncalls=$(head -$SGE_TASK_ID ncalls_list.txt | tail -1)
sims=10000

Rscript --no-save frog_power.R $ncalls $sims > output.txt



