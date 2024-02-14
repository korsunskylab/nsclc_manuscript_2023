#!/bin/bash
#SBATCH -c 15                               # Request one core
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=30G                        # Memory total in MiB (for all cores)
#SBATCH -o /n/data1/bwh/medicine/korsunsky/lab/ik97/NirMerfish/Revisions/commot/hostname_%j.out 
#SBATCH -e /n/data1/bwh/medicine/korsunsky/lab/ik97/NirMerfish/Revisions/commot/hostname_%j.err 


/n/data1/bwh/medicine/korsunsky/lab/ik97/NirMerfish/Revisions/do_commot.R $@

