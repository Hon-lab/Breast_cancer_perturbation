#!/bin/tcsh

#SBATCH --job-name=pySpade_MB231_explevel                   # job name
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

#set DIR=/project/GCRB/Hon_lab/s426305/Sequencing_data_analysis/10X/DAA271-DAA286-ER+screens/process_singlets/MB361/
set DIR=/project/GCRB/Hon_lab/s426305/Analysis/Spade_test/MB231/MB231_SM/
#set DIR=/project/GCRB/Hon_lab/s426305/Analysis/Spade_test/MB231/MB231_GWAS/
pySpade explevel \
      -t $DIR\
      -g ./query_gene.txt\
      -o ./MB231.gene_expression.txt



