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

# How to install 

Use `make` in the root directory of this package. An executable files `rf_inv` will be created under `bin` directory. 

# How to run

Use `mpirun` because this program requires parallel computation.

#### Example
`mpirun -np <N_proc> rf_inv`
* N_proc = Number of processes
* Each process has multiple McMC chains that are interacts each other through parallel tempering scheme. 

---

# Instruction

## Theory
See our [Wiki](https://github.com/akuhara/RF_INV/wiki) for details. 

## Input files

### Parameter file (params.in)
A parameter file, in which other input files and tuning parameters are specified, must exist in the current directory with the name 'params.in'. In the parameter file, all parameters must appear in the exact order as shown in the example below. Note that lines beginning with '#' are ignored.  

#### Example

```
# Example of parameter file for RF_INV (params.in)
#
# Output directory (must exist)
Out
#
#---------------------------------------------------------
# Iteration
#---------------------------------------------------------
# Iteration number in burn-in period
10000
# Iteration number after burn-in
50000
# Iteration number per saving one sample
10
# Number of chains per processor
5
# Number of non-tempered chains per processor (>= 1)
1
# Highest temperature (>= 1.0)
50.0
# Random number seed
56763298
#
#---------------------------------------------------------
# Input Data
#---------------------------------------------------------
# Number of input waveforms
3
# Ray parameter (s/km) (* Need to repeat as many time as the number of input waveforms) 
#     Note: If ray parameters are common for all input waveforms, then forward computation occurs only once.
#           Otherwise, the computation is repeated as many times as the number of input. 
0.06
0.06
0.06
# Gaussian parameter (i.e., a in G(omega) = exp(-(omega^2/(4a^2))) 
#                    (* Need to repeat as many time as the number of input waveforms)
8.0
4.0
1.0
# Sample number used in FFT (higher number takes longer time for forward computation, 
#                            but may offer better resolution in the frequency domain)
2048
# Observed RF file (* Need to repeat as many time as the number of input waveforms)
#    Note: Observed RF file must be SAC format. 
#          A header 'delta' is required to be common for all input waveforms. 
KMB06_120_150_a8.0_MC3.sac 
# Start/End time
-1.0 10.0
# Station depth (km below the sea surface)
#    Note: A station is assumed to be located on the seafloor. 
2.499
# Reference velocity file
#    Format: z(km)  Vp(km/s)  Vs(km/s)
#    Note: Need evenly spaced grids in depth
KMB06.vel.formatted
#
#---------------------------------------------------------
# Inversion mode
#---------------------------------------------------------
# Vp mode (0: Fixed at reference model, 1: Solved)
0
#
#---------------------------------------------------------
# Prior probability 
#---------------------------------------------------------
# Min./Max. # of interfaces [min., max.)
1 31
# Min./Max. of interface depth (km; below sea surface)
2.499 20.0
# Min./Max. of dVs (km/s)
-1.5 1.5
# Min./Max. of dVp (km/s)
-2. 2.0
# Min./Max. of noise sigma
0.005 0.08
#
#---------------------------------------------------------
# Proposal
#---------------------------------------------------------
# standard deviation for depth proposal
0.03
# standard deviation for dVs proposal
0.02
# standard deviation for dVp proposal
0.02
# standard deviation for noise sigma
0.005
#
#---------------------------------------------------------
# Figure
#---------------------------------------------------------
# Number of bins for depth
100
# Number of bins for Vs
50
# Number of bins for Vp
50
# Number of bins for noise sigma
50
# Number of bins for amplitudes
100
# Min./Max. amplitudes to be displayed
-0.6 0.6
# Min./Max. Vp to be displayed
0.1 8.6
# Min./Max. Vs to be displayed
0.001 5.0
```

## Output files

### num_interface.ppd
The marginal posterior probability distribution of the number of layer interfaces. 

#### Format 
|1st column| 2nd column|
|:--:|:--:|
|# of layer interfaces|Probability|

### interface_depth.ppd
The marginal posterior probability distribution of interface depths.  

#### Format
|1st column| 2nd column|
|:--:|:--:|
|Depth (km)|Probability|

### vs_z.ppd
The marginal posterior probability distribution of S-wave velocity profile.

#### Format 
|1st column| 2nd column| 3rd column|
|:--:|:--:|:--:|
|Vs (km/s)|Depth (km)|Probability|

### vp_z.ppd
The marginal posterior probability distribution of P-wave velocity profile.

#### Format 
|1st column| 2nd column| 3rd column|
|:--:|:--:|:--:|
|Vp (km/s)|Depth (km)|Probability|

### sigma.ppd
The marginal posterior probability distribution of noise standard deviation.

#### Format
|1st column| 2nd column|
|:--:|:--:|
|Standard deviation|Probability|

