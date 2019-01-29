# RF_INV

Receiver function inversion by reversible-jump Markov-chain Monte Carlo

(c) 2018 Takeshi Akuhara (Email: akuhara @ eri.u-tokyo.ac.jp)

IMPORTANT NOTE (1/18/2019) : THIS PROGRAM IS UNDER DEVELOPMENT. ANY BUG REPORT IS WELCOME, BUT I DON'T GURANTEE THAT THE PROGRAM WORKS CORRECTLY AT THIS STAGE. A COMPLETE VERSION WILL BE RELEASED IN THE NEAR FUTURE.

# Features

There are many softwares and literatures for trans-dimensional inversion of receiver functions. This inversion code has following originalities.

* Applicable to OBS data
    * Model can include sea water on its top. A station is assumed to locate on the seafloor.  
* Multiple input traces
    * Can asign different ray parameter and Gaussian-filter for each trace.    
* Parallel tempering
    * More efficient than conventional MCMC.
    
See [Wiki](https://github.com/akuhara/RF_INV/wiki) for more details.

# Terms of Use
* Please clarify the URL of the GitHub repository (https://github.com/akuhara/RF_INV) and developer's name (Takeshi Akuhara) when you make any presentation or publish articles using this program.
* This program is licensed under the GNU General Public License v3.0.

# Requirement

The following libraries must be installed. 

* [FFTW library](http://fftw.org/)
* [LAPACK library](http://www.netlib.org/lapack/)

The appropriate locations of the above libraries must be specified in `Makefile`. 

* [Open MPI](https://www.open-mpi.org/)

It is recommended to compile Open MPI with GNU Fortran compiler (i.e, `gfortran`). In that case, you can use Makefile without changing compiler options.

---

# Manual 

Check [Wiki](https://github.com/akuhara/RF_INV/wiki).
