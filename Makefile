MF90 = mpifort

#--------------
# For fast computation
FFLAGS = -ffast-math -march=native -mtune=native -O3 -fno-range-check
#--------------

#--------------
# For debug

#FFLAGS = -pg -Wall -pedantic -std=f2003 -fbounds-check -O0 -Wuninitialized \
            -ffpe-trap=invalid,zero,overflow -fbacktrace \
            -fno-range-check 
#--------------

#--------------
# Locations of library
FFTW   = -I/usr/local/include -lfftw3
LAPACK = -llapack -lblas
#--------------


#
# Do not change below
#

BINDIR  = ./bin
TARGET1 = $(BINDIR)/rf_inv
OBJS1   = src/rf_inv.o src/params.o src/mt19937.o \
         src/fftw.o src/model.o src/sort.o src/likelihood.o src/forward.o \
         src/pt_mcmc.o src/mcmc_out.o src/math.o src/prior.o
TARGET2 = $(BINDIR)/make_syn
OBJS2   = src/make_syn.o src/model.o src/params.o src/sort.o src/mt19937.o \
	  src/likelihood.o src/forward.o src/fftw.o \
	  src/pt_mcmc.o src/math.o src/prior.o




all: $(TARGET1) $(TARGET2)

$(TARGET1): $(OBJS1)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(MF90) $(FFLAGS) $(FFTW) $(LAPACK) $^ -o $@

$(TARGET2): $(OBJS2)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(MF90) $(FFLAGS) $(FFTW) $(LAPACK) $^ -o $@

src/rf_inv.o: params.mod mt19937.mod fftw.mod model.mod likelihood.mod \
              forward.mod pt_mcmc.mod mcmc_out.mod
src/fftw.o: params.mod
src/model.o: params.mod mt19937.mod sort.mod math.mod prior.mod
src/likelihood.o: params.mod model.mod forward.mod
src/forward.o: params.mod fftw.mod model.mod
src/pt_mcmc.o: params.mod mt19937.mod model.mod likelihood.mod math.mod \
	       prior.mod
src/mcmc_out.o: params.mod
src/make_syn.o: params.mod model.mod pt_mcmc.mod likelihood.mod forward.mod \
	        fftw.mod mt19937.mod math.mod
src/math.o: mt19937.mod

clean:
	rm -f *.mod bin/inv_PT_RF src/*.o *.o


#------------------------------------------------------------
# Pattern rule
#------------------------------------------------------------
%.o: %.f90
	$(MF90) $(FFLAGS) -c $< $(FFTW) $(LAPACK) -o $*.o 
%.mod: src/%.f90 src/%.o
	@:
