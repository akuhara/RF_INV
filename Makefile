MF90 = mpif90
#FFLAGS = -ffast-math -march=native -mtune=native -O3 -fno-range-check
FFLAGS = -pg -Wall -pedantic -std=f95 -fbounds-check -O -Wuninitialized \
            -ffpe-trap=invalid,zero,overflow -fbacktrace \
            -fno-range-check 

FFTW  = -I/usr/local/include -lfftw3



BINDIR  = ./bin
TARGET  = $(BINDIR)/rf_inv
OBJS = src/rf_inv.o src/params.o src/mt19937.o



all: $(TARGET)

$(TARGET): $(OBJS)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(MF90) $(FFLAGS) $(FFTW) $^ -o $@


src/rf_inv.o: params.mod

clean:
	rm -f *.mod bin/inv_PT_RF src/*.o *.o


#------------------------------------------------------------
# Pattern rule
#------------------------------------------------------------
$(OBJS): %.o: %.f90
	$(MF90) $(FFLAGS) -c $< $(LIBRARY) $(INCLUDE) -o $*.o 
%.mod: src/%.f90 src/%.o
	@:
