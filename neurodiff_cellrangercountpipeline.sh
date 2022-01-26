#!/bin/bash
# Job name:
#SBATCH --job-name=neurodiff-ctrl
#
# Project:
#SBATCH --account=projectname
#
# Wall clock limit:
#SBATCH --time=28:00:00
#
# Max memory usage: Size is a number plus M (megabyte) or G (gigabyte), e.g., 3M or 5G
#SBATCH --mem-per-cpu=4G
#
# Number of tasks (cores): this is added to make it easier for you to do exercises
#SBATCH --ntasks=12
#TASK STATUS EMAILS
#SBATCH --mail-user=email@email.com
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## Set up job environment: (this is done automatically behind the scenes)
## (make sure to comment '#' or remove the following line 'source ...')
# source /cluster/bin/jobsetup

module purge   # clear any inherited modules
#module load SoftWare/Version #nb: 'Version' is mandatory! There are no default versions of modules on Saga!

# It is also recommended to to list loaded modules, for easier debugging:
module list

## Copy input files to the work directory:

## Make sure the results are copied back to the submit directory:
# chkfile MyResultFile
# chkfile is replaced by 'savefile' on Saga

#paths to folder
HHL7VBGXB=~/Raw-Data_HHL7VBGXB
H72LWBGXB=~/Raw-Data_H72LWBGXB
## Do some work:
export PATH=~/software/cellranger-3.1.0:$PATH
cellranger count --id=sc1d0-ctrl \
                 --transcriptome=/REFERENCE_GENOMES/SCRNA_10X/refdata-cellranger-GRCh38-3.0.0 \
                 --fastqs=$HHL7VBGXB/sc1d0-ctrl,$H72LWBGXB/sc1d0-ctrl\
                 --sample=sc1d0-ctrl \
                 --expect-cells=1400
echo done

cellranger count --id=sc2d0-ctrl  \
                 --transcriptome=/REFERENCE_GENOMES/SCRNA_10X/refdata-cellranger-GRCh38-3.0.0 \
                 --fastqs=$HHL7VBGXB/sc2d0-ctrl,$H72LWBGXB/sc2d0-ctrl \
                 --sample=sc2d0-ctrl \
                 --expect-cells=1700
echo done
cellranger aggr --id aggr_d0_ctrl --csv aggregationd0.csv > aggr-log.txt
echo done

cellranger count --id=sc1d7-ctrl  \
                 --transcriptome=/REFERENCE_GENOMES/SCRNA_10X/refdata-cellranger-GRCh38-3.0.0 \
                 --fastqs=$HHL7VBGXB/sc1d7-ctrl,$H72LWBGXB/sc1d7-ctrl \
                 --sample=sc1d7-ctrl  \
                 --expect-cells=1800

echo done

cellranger count --id=sc2d7-ctrl \
                 --transcriptome=/REFERENCE_GENOMES/SCRNA_10X/refdata-cellranger-GRCh38-3.0.0 \
                 --fastqs=$HHL7VBGXB/sc2d7-ctrl,$H72LWBGXB/sc2d7-ctrl \
                 --sample=sc2d7-ctrl\
                 --expect-cells=1700
echo done
cellranger aggr --id aggr_d7_ctrl --csv aggregationd7.csv > aggr-log.txt
echo done 

cellranger count --id=sc1d13-ctrl \
                 --transcriptome=/REFERENCE_GENOMES/SCRNA_10X/refdata-cellranger-GRCh38-3.0.0 \
                 --fastqs=$HHL7VBGXB/sc1d13-ctrl,$H72LWBGXB/sc1d13-ctrl \
                 --sample=sc1d13-ctrl \
                 --expect-cells=1200
echo done


cellranger count --id=sc2d13-ctrl  \
                 --transcriptome=/REFERENCE_GENOMES/SCRNA_10X/refdata-cellranger-GRCh38-3.0.0 \
                 --fastqs=$HHL7VBGXB/sc2d13-ctrl,$H72LWBGXB/sc2d13-ctrl \
                 --sample=sc2d13-ctrl  \
                 --expect-cells=1700
echo done
cellranger aggr --id aggr_d13_ctrl --csv aggregationd13.csv > aggr-log.txt
echo done 

cellranger count --id=sc1d20-ctrl \
                 --transcriptome=/REFERENCE_GENOMES/SCRNA_10X/refdata-cellranger-GRCh38-3.0.0 \
                 --fastqs=$HHL7VBGXB/sc1d20-ctrl,$H72LWBGXB/sc1d20-ctrl \
                 --sample=sc1d20-ctrl \
                 --expect-cells=2100

echo done

cellranger count --id=sc2d20-ctrl  \
                 --transcriptome=/REFERENCE_GENOMES/SCRNA_10X/refdata-cellranger-GRCh38-3.0.0 \
                 --fastqs=$HHL7VBGXB/sc2d20-ctrl,$H72LWBGXB/sc2d20-ctrl \
                 --sample=sc2d20-ctrl \
                 --expect-cells=1900
echo done 

cellranger aggr --id aggr_d20_ctrl --csv aggregationd20.csv > aggr-log.txt

echo done 
