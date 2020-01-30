#$ -P sabeti_lab
#$ -wd /idi/sabeti-scratch/jvitti/cms/cms/
#$ -t 1-4

export PATH=/idi/sabeti-scratch/jvitti/cms/cms:$PATH
export PATH=/home/unix/vitti/miniconda3/bin:$PATH
source activate /home/unix/vitti/miniconda3/envs/cms-env3
source /broad/software/scripts/useuse
MYFILE=/idi/sabeti-scratch/jvitti/cms/cms/redo_delihh_nd_1_bins_taskargs.txt
SAMPLE=$(awk "NR==$SGE_TASK_ID" $MYFILE)
python /idi/sabeti-scratch/jvitti/power.py composite_sims $SAMPLE
