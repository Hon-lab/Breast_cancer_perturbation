#!/bin/tcsh

#SBATCH --job-name=pySpade_local_CMpilot2_cluster0        # job name
#SBATCH --partition=256GB,256GBv1                         # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=60-00:00:00                                # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                         # standard output file name
#SBATCH --error=serialJob.%j.time                         # standard error output file name
#SBATCH --mail-user=Yihan.Wang@utsouthwestern.edu         # specify an email address
#SBATCH --mail-type=end                                   # send email when job status change (start, end, abortion and etc.)

setenv PATH /project/GCRB/Hon_lab/s426305/Conda/py37-cluster2/bin:$PATH
set DIR=/project/GCRB/Hon_lab/s426305/Analysis/IGVF/20230612_CMPilot2_cluster0_ver0/pySpade/
set OBS_DIR=$DIR/CMPilot2_cluster0_ver0_DEobs/
set RAND_DIR=$DIR/CMPilot2_cluster0_ver0_DErand/
set SGRNA_DICT=$DIR/sgRNA_dict_hg38.txt

    
 pySpade local \
 	-f $DIR \
	-d $OBS_DIR\
	-t $RAND_DIR\
	-s $SGRNA_DICT\
	-o ./CMpilot2_cluster0_local_hit.csv
