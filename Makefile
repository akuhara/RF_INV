MF90 = mpif90
#FFLAGS = -ffast-math -march=native -mtune=native -O3 -fno-range-check
FFLAGS = -pg -Wall -pedantic -std=f95 -fbounds-check -O -Wuninitialized \
            -ffpe-trap=invalid,zero,overflow -fbacktrace \
            -fno-range-check 

FFTW   = -I/usr/local/include -lfftw3
LAPACK = -llapack -lblas


BINDIR  = ./bin
TARGET = $(BINDIR)/rf_inv
OBJS   = src/rf_inv.o src/params.o src/mt19937.o src/read_obs.o \
         src/fftw.o src/model.o src/sort.o src/likelihood.o src/forward.o \
         src/pt_mcmc.o


all: $(TARGET)

$(TARGET): $(OBJS)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(MF90) $(FFLAGS) $(FFTW) $(LAPACK) $^ -o $@

src/rf_inv.o: params.mod mt19937.mod fftw.mod model.mod likelihood.mod \
              forward.mod pt_mcmc.mod
src/read_obs.o: params.mod
src/fftw.o: params.mod
src/model.o: params.mod mt19937.mod sort.mod forward.mod
src/likelihood.o: params.mod model.mod
src/forward.o: params.mod fftw.mod
src/pt_mcmc.o: params.mod mt19937.mod model.mod likelihood.mod

clean:
	rm -f *.mod bin/inv_PT_RF src/*.o *.o


#------------------------------------------------------------
# Pattern rule
#------------------------------------------------------------
$(OBJS): %.o: %.f90
	$(MF90) $(FFLAGS) -c $< $(FFTW) $(LAPACK) -o $*.o 
%.mod: src/%.f90 src/%.o
	@:
