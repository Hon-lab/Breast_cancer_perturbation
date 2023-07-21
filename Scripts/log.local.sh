#!/bin/tcsh

#SBATCH --job-name=pySpade_local_MB231_GWAS	          # job name
#SBATCH --partition=256GB,256GBv1                         # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=60-00:00:00                                # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                         # standard output file name
#SBATCH --error=serialJob.%j.time                         # standard error output file name
#SBATCH --mail-user=Yihan.Wang@utsouthwestern.edu         # specify an email address
#SBATCH --mail-type=end                                   # send email when job status change (start, end, abortion and etc.)

setenv PATH /project/GCRB/Hon_lab/s426305/Conda/py37-cluster2/bin:$PATH
set DIR=/project/GCRB/Hon_lab/s426305/Analysis/Spade_test/MB231/
set OBS_DIR=$DIR/MB231_GWAS_DEobs/
set RAND_DIR=$DIR/MB231_GWAS_DErand/
set SGRNA_DICT=/project/GCRB/Hon_lab/s426305/Analysis/Mosaic-seq/MB231-dCas9-KRAB+YWsg2P1/annotation/enhancer_regions_hg38.txt

    
 pySpade local \
 	-f $DIR \
	-d $OBS_DIR\
	-t $RAND_DIR\
	-s $SGRNA_DICT\
	-o ./MB231_GWAS_local_hit.csv
