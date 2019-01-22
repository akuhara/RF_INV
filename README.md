# RF_INV

Receiver function inversion by reversible-jump Markov-chain Monte Carlo

(c) 2018 Takeshi Akuhara

IMPORTANT NOTE (1/18/2019) : THIS PROGRAM IS UNDER DEVELOPMENT. ANY BUG REPORT IS WELCOME, BUT I DON'T GURANTEE THAT THE PROGRAM WORKS CORRECTLY AT THIS STAGE. A COMPLETE VERSION WILL BE RELEASED IN THE NEAR FUTURE.

# Features

There are many literature exists for trans-dimensional inversion of receiver functions. Nevertheless I would like to make this code open becuase of the following originalities.

* Applicable to OBS data
    * Model can include sea water on its top. A station is assumed to locate on the seafloor.  
* Multiple input trace
    * You can asign different ray parameter and Gaussian-filter for each trace.    
* Parallel tempering
    * More efficient than conventional MCMC.

# Requirement

The following libraries must be installed. 

* [FFTW library](http://fftw.org/)
* [LAPACK library](http://www.netlib.org/lapack/)

The appropriate locations of the above libraries must be specified in `Makefile`. 

* [Open MPI](https://www.open-mpi.org/)

It is recommended to complie Open MPI with GNU fortran compilar (i.e, `gfortran`). In that case, you can use Makefile without changing compiler options.

# How to install 

`make`
