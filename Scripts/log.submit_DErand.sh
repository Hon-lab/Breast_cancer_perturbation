#!/bin/tcsh

foreach NUM (600\
	     700\
	     800\
	     900\
	     1000\
	     1100\
	     1200\
	     1300\
	     1400\
	     1500\
	     1600\
	     1700\
	     1800\
	     1900\
	     2000\
	    )

    sed -e "s/NUM/$NUM/g" log.DErand.sh\
    > log.DErand_$NUM.sh

    sbatch log.DErand_$NUM.sh
end
