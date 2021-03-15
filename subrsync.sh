#!/bin/bash
## user:=mrajeev
## cluster:=della.princeton.edu
## dirput:=analysis/slurm
## -r recurse through sub-directories
## -L transform symlink into reference file/dir
## -v be verbose about printing messages
## -z compresses data before transfer and decompresses after transfer
## -t pass the timestamp when syncing

rsync -rLvzt --update --exclude '*.git' --exclude '.Rproj*' --exclude '*archive*' --exclude analysis/out --exclude 'data-raw*'  --exclude "*.DS_Store*" ~/Documents/Projects/dynamicSD mrajeev@della.princeton.edu:~/


# Push up the outputs to scratch (so always working from there for outputs)
# Should already be there 
rsync -rLvzt --update --exclude "*.DS_Store*" --exclude "*archive*" ~/Documents/Projects/dynamicSD/analysis/out mrajeev@della.princeton.edu:/scratch/gpfs/mrajeev/dynamicSD/analysis/
