SUFFIX=.f

##########################
# FORTRAN CONFIGURATIONS #
##########################
# Compiler
FC=ifort

# Applied to all fortran files
FC_ALL = -xHost -ip -align -no-wrap-margin

# Parallelization flags
FC_PARALLEL = -qopenmp -parallel
F_FLAGS =  -O3 $(FC_PARALLEL) $(FC_ALL)

############################
# LIBARRIES CONFIGURATIONS #
############################
LIB_BLAS = -lmkl_blas95_lp64
LIB_LAPACK = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
LIB_OMP = -liomp5 -lpthread

INCLUDES_MKL = -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include/fftw

LIBS  = $(LIB_GPU) $(LIB_BLAS) $(LIB_LAPACK) $(LIB_OMP) -lrt
INCLUDES = $(INCLUDES_MKL)

#-----------------------------------------------------------------------
# general rules
#-----------------------------------------------------------------------

SOURCE = constants.o \
	     grade.o \
         IO_stuff.o \
         roots.o \
         initial_state.o \
         wf.o \
	     wfcont.o \
         ensemble.o \
	     time_propagator.o \
         QW_dynamics.o \
	     width.o \
	     main.o \

a: $(SOURCE)  
	rm -f a
	$(FC) $(INCLUDES) -o a $(SOURCE) $(LIBS) 
	-rm -f *.log

# Program runs very slowly with this                                                                                                                                    
safe: F_FLAGS += -check all -traceback -fstack-protector -assume protect_parens -implicitnone -warn all,noexternal -fpe-all=0
safe: a

.f.o:
	$(FC) -free $(F_FLAGS) $(INCLUDES) -c $*$(SUFFIX)
clean: 
	-rm -f *.o *.mod *.dat ; touch *.f




