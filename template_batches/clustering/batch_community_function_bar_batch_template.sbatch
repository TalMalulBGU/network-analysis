#!/bin/bash

################################################################################################
### sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --partition main			### specify partition name where to run a job. short: 7 days limit; gtx1080: 7 days; debug: 2 hours limit and 1 job at a time
#SBATCH --time 3-00:00:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name tm_ddg_met			### name of the job
#SBATCH --output job-%J_%a.out			### output log for running job - %J for job number
##SBATCH --mail-user=talmalu@post.bgu.ac.il	### user's email for sending job status messages
##SBATCH --mail-type=NONE			### conditions for sending the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
#SBATCH --mem=80G				### ammount of RAM memory
#SBATCH --cpus-per-task=6			### number of CPU cores

### Print some data to output file ###
echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

### Start your code below ####
module load anaconda				### load anaconda module (must be present when working with conda environments)
source activate NewNetworkAnalysis			### activate a conda environment, replace my_env with your conda environment
FILES=(###JSON_FILES###)
python /home/talmalu/thesis/projects/python/NetworkAnalysis/Scripts/python/clustering/community_function_bar.py \
--protein ###PROTEIN### \
--xgmml_file ###XGMML_FILE### \
--alignment_weight ###ALIGNMENT_WEIGHT### \
--algorithm_file ${FILES[$SLURM_ARRAY_TASK_ID]} \
--function_column ###FUNCTION_COLUMN### \
--relevance_function_xlsx ###RELEVANCY_XLSX### \
--image_title ###IMAGE_TITLE### \
--min_in_cluster ###MIN_IN_CLUSTER### \
--database ###DATABASE### \
--out_file ${FILES[$SLURM_ARRAY_TASK_ID]}.function.png
