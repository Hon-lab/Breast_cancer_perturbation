#!/bin/tcsh

#SBATCH --job-name=pySpade_DEobs_MB231_GWAS               # job name
#SBATCH --partition=256GB,256GBv1,384GB                   # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=60-00:00:00                                # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                         # standard output file name
#SBATCH --error=serialJob.%j.time                         # standard error output file name
#SBATCH --mail-user=Yihan.Wang@utsouthwestern.edu         # specify an email address
#SBATCH --mail-type=end                                   # send email when job status change (start, end, abortion and etc.)


setenv PATH /project/GCRB/Hon_lab/s426305/Conda/py37-cluster2/bin:$PATH
echo 'Program is running with the current python version:'
which python
python --version

set TRANS_FILE=/project/GCRB/Hon_lab/s426305/Analysis/Spade_test/MB231/MB231_GWAS/Singlet_sub_df.h5
set SG_FILE=/project/GCRB/Hon_lab/s426305/Analysis/Spade_test/MB231/MB231_GWAS/Singlet_sgRNA_df.h5
set ANNOT_FILE=/project/GCRB/Hon_lab/s426305/Analysis/Mosaic-seq/MB231-dCas9-KRAB+YWsg1P3/enhancer_regions_combined.txt

mkdir MB231_GWAS_DEobs
pySpade DEobs\
	-t $TRANS_FILE\
	-s $SG_FILE\
	-d $ANNOT_FILE\
	-r $SLURM_CPUS_ON_NODE\
	-n 'cpm'\
	-o ./MB231_GWAS_DEobs/
