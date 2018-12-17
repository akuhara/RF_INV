MF90 = mpif90
FFLAGS = -ffast-math -march=native -mtune=native -O3 -fno-range-check
#FFLAGS = -pg -Wall -pedantic -std=f95 -fbounds-check -O -Wuninitialized \
            -ffpe-trap=invalid,zero,overflow -fbacktrace \
            -fno-range-check 

FFTW  = -I/usr/local/include -lfftw3



BINDIR  = ./bin
TARGET  = $(BINDIR)/rf_inv
OBJ_MPI = src/rjmcmc_inv_rf.o
OBJ_F90 = src/fwd_rf.o src/fwd_seis.o src/pt.o src/params.o src/mt19937.o
OBJ_F   = src/SVDlibNR.o
OBJ_ALL = $(OBJ_MPI) $(OBJ_F90) $(OBJ_F)


all: $(TARGET)

$(TARGET): $(OBJ_ALL) 
	$(MF90) $(FFLAGS)  -o $(TARGET) $(OBJ_ALL) $(LIBRARY) $(INCLUDE)


src/rjmcmc_inv_rf.o: params.mod mt19937.mod
src/fwd_rf.o: params.mod
src/fwd_seis.o: params.mod

clean:
	rm -f *.mod bin/inv_PT_RF src/*.o *.o


#------------------------------------------------------------
# Pattern rule
#------------------------------------------------------------
$(OBJ_F90): %.o: %.f90
	$(MF90) $(FFLAGS) -c $< $(LIBRARY) $(INCLUDE) -o $*.o 
$(OBJ_F): %.o: %.f
	$(MF90) $(FFLAGS) -c $< $(LIBRARY) $(INCLUDE) -o $*.o 
$(OBJ_MPI): %.o: %.f90
	cp $*.f90 $*.F90
	$(MF90) $(FFLAGS) -c $*.F90 $(LIBRARY) $(INCLUDE) -o $*.o
%.mod: src/%.f90 src/%.o
	@:
