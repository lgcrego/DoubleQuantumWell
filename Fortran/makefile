#

SUFFIX=.f90

#FC=gfortran
#FFLAGS = -ffree-form -ffree-line-length-none 

FC=/opt/intel/bin/ifort
FFLAGS = -O3 -free #(Kewin sua bichona, eu comentei aqui!)

LIB_BLAS    = -L/opt/intel/mkl/lib/intel64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
LIB_LAPACK  = -L/opt/intel/mkl/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

INCS_MKL   = -I/opt/intel/mkl/include/intel64/lp64

LIB    =  $(LIB_BLAS) $(LIB_LAPACK)
INCS   =  $(INCS_MKL)

#-----------------------------------------------------------------------
# general rules
#-----------------------------------------------------------------------

SOURCE =constants.o \
	grade.o \
        roots.o \
        ground_state.o \
        wf.o \
	wfcont.o \
	orthogonalbasis.o \
	width.o \
	dirs.o \
	main.o \

a: $(SOURCE)  
	rm -f a-high2
	$(FC) $(INCS) -o a-high2 $(SOURCE) $(LIB) 
	-rm -f *.log
.f.o:
	$(FC) $(FFLAGS) $(INCS) -c $*$(SUFFIX)
clean: 
	-rm -f *.o *.mod *.dat ; touch *.f




