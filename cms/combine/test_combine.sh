#$ -l h_vmem=3g
#$ -cwd
#$ -P sabeti_lab
#$ -l h_rt=48:00:00
#$ -M vitti@broadinstitute.org -m eas
export PATH=/idi/sabeti-scratch/jvitti/miniconda3/bin:$PATH
source activate /idi/sabeti-scratch/jvitti/miniconda3/envs/cms-env3
source /broad/software/scripts/useuse

/idi/sabeti-scratch/jvitti/cms/cms/combine/collate_scores_v2  /idi/sabeti-scratch/jvitti/1kg_composite/gw/unmasked_aa_v2/chr19_CHS.cms.gw.out.122217_b  /idi/sabeti-scratch/jvitti/1kg_composite/runparams/params_3.122217  /idi/sabeti-scratch/jvitti/1kg_components/pairs/unmasked_aa_v2/chr19_CHS_YRI.pair.122217  /idi/sabeti-scratch/jvitti/1kg_components/pairs/unmasked_aa_v2/chr19_CHS_CEU.pair.122217  /idi/sabeti-scratch/jvitti/1kg_components/pairs/unmasked_aa_v2/chr19_CHS_BEB.pair.122217
