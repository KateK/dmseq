#!/bin/bash

#PBS -N aggregateOutput
#PBS -q short
#PBS -l nodes=1:ppn=12
#PBS -o /home3/struck/scripts/log/
#PBS -e /home3/struck/scripts/log/
#PBS -j oe

/home3/struck/Analysis_Projects/DMseq/bin/run_summarizeMISOoutput.R