# Sample synthetic dataset

This direcotry contains synthetic dataset for testing RF_INV. Note that the dataset here is paricuraly 'easy' dataset that are tuned such that one can obtain resonable solutions quickly. For application to real dataset, you will need ~100 times more iterations. DON'T copy & paste the parameter file here for your dataset without noticing this. 

## Usage
1.  `mkdir rslt`
1.  `mpirun -np 20 ../bin/rf_inv params.in`
1.  `python ../util/plot.py params.in`
