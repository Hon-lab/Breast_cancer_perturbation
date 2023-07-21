#!/bin/tcsh

#SBATCH --job-name=Spade_process_MB231_GWAS               # job name
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

set TRANS_DIR=/project/GCRB/Hon_lab/s426305/Sequencing_data_analysis/10X/YW104-YW109/10x/MB231-YWsg1_combine_nova_new/outs/
set SG_FILE=/project/GCRB/Hon_lab/s426305/Sequencing_data_analysis/10X/MB231-GWAS_SM-combine/fba_sgRNA/curve/MB231_GWAS.pkl.gz
set MULTIPLEX_FILE=/project/GCRB/Hon_lab/s426305/Analysis/Spade_test/MB231/HTO.txt

mkdir MB231_GWAS
pySpade process\
	-f $TRANS_DIR\
	-s $SG_FILE\
	-m $MULTIPLEX_FILE\
	-o ./MB231_GWAS/
