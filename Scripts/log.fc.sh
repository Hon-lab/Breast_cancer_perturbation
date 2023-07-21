#!/bin/tcsh

#SBATCH --job-name=pySpade-fc-MB231                       # job name
#SBATCH --partition=256GB,256GBv1                         # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=60-00:00:00                                # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                         # standard output file name
#SBATCH --error=serialJob.%j.time                         # standard error output file name
#SBATCH --mail-user=Yihan.Wang@utsouthwestern.edu         # specify an email address
#SBATCH --mail-type=end                                   # send email when job status change (start, end, abortion and etc.)

setenv PATH /project/GCRB/Hon_lab/s426305/Conda/py37-cluster2/bin:$PATH

set GWAS_DICT=/project/GCRB/Hon_lab/s426305/Analysis/Mosaic-seq/MB231-dCas9-KRAB+YWsg1P3/enhancer_regions_combined.txt
set MB231_GWAS_DIR=/project/GCRB/Hon_lab/s426305/Analysis/Spade_test/MB231/MB231_GWAS/

#GWAS
pySpade fc \
	-t $MB231_GWAS_DIR\
	-d $GWAS_DICT\
	-r ./GWAS_repression_query.txt\
	-o ./Jupyter_plots/
