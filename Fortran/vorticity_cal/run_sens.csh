#!/bin/csh
#$ -V
#$ -cwd
#$ -S /bin/csh
#$ -N made_file
#$ -o out.runsens
#$ -j y
#$ -q normal
#$ -P hrothgar
#$ -pe fill 12

date
./made_file
date
