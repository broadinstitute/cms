#$ -P sabeti_lab
#$ -wd /idi/sabeti-scratch/jvitti/cms/cms/
#$ -t 1-3360

export PATH=/idi/sabeti-scratch/jvitti/cms/cms:$PATH
export PATH=/home/unix/vitti/miniconda3/bin:$PATH
source activate /home/unix/vitti/miniconda3/envs/cms-env3
source /broad/software/scripts/useuse
MYFILE=/idi/sabeti-scratch/jvitti/cms/cms/fpr_debug_composite_nsl_noihs_taskargs.txt
SAMPLE=$(awk "NR==$SGE_TASK_ID" $MYFILE)
python /idi/sabeti-scratch/jvitti/power.py fpr $SAMPLE
