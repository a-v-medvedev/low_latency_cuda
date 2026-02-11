FC = nvfortran
CXX = nvcc

# NOTE: tuned for MareNostrum-5 Accelerated partition
# module load nvidia-hpc-sdk/25.7

F90ARCHFLAGS=-acc=gpu -gpu=cc90
CXXARCHFLAGS=-gencode arch=compute_90,code=sm_90
LINKARCHFLAGS=-acc=gpu -gpu=cc90

F90OPTFLAGS=-O3     #  -O0 -g 
CXXOPTFLAGS=-O3     #  -O0
LINKOPTFLAGS=       #  -g 

F90GENERICFLAGS=-Mpreprocess -Minfo=accel -module .
CXXGENERICFLAGS=--x cu
LINKGENERICFLAGS=

LINKERLIBS=-lstdc++ -L$(NVHPC_ROOT)/cuda/lib64 -lcudart

# NOTE: comment this out to build a binary that is free from CUTENSOREX
F90EXTRAFLAGS= -DWITH_CUTENSOREX
CXXEXTRAFLAGS=
LINKEXTRALIBS= -cudalib=cutensor


all: sc_test

sc_test.o: sc_test.f90
	$(FC) -c -o $@ $(F90GENERICFLAGS) $(F90EXTRAFLAGS) $(F90ARCHFLAGS) $(F90OPTFLAGS) $^

print_norms.o: print_norms.f90
	$(FC) -c -o $@ $(F90GENERICFLAGS) $(F90EXTRAFLAGS) $(F90ARCHFLAGS) $(F90OPTFLAGS) $^

sum_prefix_custom.o: sum_prefix_custom.f90
	$(FC) -c -o $@ -cuda $(F90GENERICFLAGS) $(F90EXTRAFLAGS) $(F90ARCHFLAGS) $(F90OPTFLAGS) $^ 
# NOTE: -cuda only because of the "device" attribute

sum_prefix_custom_cxx.o: sum_prefix_custom_cxx.cpp
	$(CXX) -c $(CXXARCHFLAGS) $(CXXGENERICFLAGS) $(CXXEXTRAFLAGS) $(CXXOPTFLAGS) $^

sc_test: sc_test.o sum_prefix_custom.o sum_prefix_custom_cxx.o print_norms.o
	$(FC) $(LINKARCHFLAGS) $(LINKOPTFLAGS) -o $@ $^ $(LINKERLIBS) $(LINKEXTRALIBS)

sum_prefix_custom.o: sum_prefix_custom_cxx.o
sc_test.o: sum_prefix_custom.o print_norms.o

clean:
	rm -f *.o *.mod sc_test

