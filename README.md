# RF_INV

Transdimensional inversion of receiver function waveforms by reversible-jump Markov-chain Monte Carlo

(c) 2018-2019 Takeshi Akuhara (Email: akuhara @ eri.u-tokyo.ac.jp)

Any bug report and suggestions are welcome!

# Features
 
* Applicable to OBS data
    * Model can include sea water on its top. A station is assumed to locate on the seafloor.  
* Multiple input traces
    * Can asign different ray parameter and Gaussian-filter for each trace.
* Can use S receiver functions
    * Joint inversion of P and S receiver functions is also possible.
* Parallel tempering
    * More efficient than conventional MCMC.
* Invert for velocity perturbation
    * Non-uniqueness of inversion can be mitigated by constraint from a reference velocity model.
    
See [Wiki](https://github.com/akuhara/RF_INV/wiki) for more details.

# Terms of Use
* Please clarify the URL of the GitHub repository (https://github.com/akuhara/RF_INV) and developer's name (Takeshi Akuhara) when you make any presentation or publish articles using this program.
* This program is licensed under the GNU General Public License v3.0.

# Requirements
* [FFTW library](http://fftw.org/)
* [LAPACK library](http://www.netlib.org/lapack/)
* [Open MPI](https://www.open-mpi.org/)

---

# Manual 

A manual is available [here](https://github.com/akuhara/RF_INV/wiki).
