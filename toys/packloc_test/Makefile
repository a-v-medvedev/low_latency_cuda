BINARY=packloc_test

FC = nvfortran
CXX = nvcc

# NOTE: tuned for MareNostrum-5 Accelerated partition
# module load nvidia-hpc-sdk/25.7

F90ARCHFLAGS=-acc=gpu -gpu=cc90
CXXARCHFLAGS=-gencode arch=compute_90,code=sm_90
LINKARCHFLAGS=-acc=gpu -gpu=cc90

F90OPTFLAGS=-O3    
CXXOPTFLAGS=-O3    
LINKOPTFLAGS=       

#F90OPTFLAGS=-O0 -g 
#CXXOPTFLAGS=-O0 -g -G
#LINKOPTFLAGS=-g 

F90GENERICFLAGS=-Mpreprocess -Minfo=accel -cuda -module .
CXXGENERICFLAGS=--x cu --std=c++20
LINKGENERICFLAGS= -cuda

LINKERLIBS=-lstdc++ -L$(NVHPC_ROOT)/cuda/lib64 -lcudart

# NOTE: comment this out to build a binary that is free from CUTENSOREX
F90EXTRAFLAGS= 
CXXEXTRAFLAGS=
LINKEXTRALIBS=

ifeq ($(WITH_PACKLOC_CUTENSOREX),TRUE)
override F90EXTRAFLAGS+= -DWITH_PACKLOC_CUTENSOREX=TRUE
override LINKEXTRALIBS+= -cudalib=cutensor
endif
ifeq ($(WITH_PACKLOC_CUSTOM_OLD),TRUE)
override F90EXTRAFLAGS+= -DWITH_PACKLOC_CUSTOM_OLD=TRUE
endif
ifeq ($(WITH_PACKLOC_CUSTOM_NEW),TRUE)
override F90EXTRAFLAGS+= -DWITH_PACKLOC_CUSTOM_NEW=TRUE
endif

all: $(BINARY)

packloc_test_fortran.o: packloc_test_fortran.f90
	$(FC) -c -o $@ $(F90GENERICFLAGS) $(F90EXTRAFLAGS) $(F90ARCHFLAGS) $(F90OPTFLAGS) $<

sum_prefix_custom.o: sum_prefix_custom.f90
	$(FC) -c -o $@ -cuda $(F90GENERICFLAGS) $(F90EXTRAFLAGS) $(F90ARCHFLAGS) $(F90OPTFLAGS) $<
# NOTE: -cuda only because of the "device" attribute

sum_prefix_custom_cxx.o: sum_prefix_custom_cxx.cpp
	$(CXX) -c $(CXXARCHFLAGS) $(CXXGENERICFLAGS) $(CXXEXTRAFLAGS) $(CXXOPTFLAGS) $<

packloc_custom.o: packloc_custom.f90
	$(FC) -c -o $@ -cuda $(F90GENERICFLAGS) $(F90EXTRAFLAGS) $(F90ARCHFLAGS) $(F90OPTFLAGS) $<
# NOTE: -cuda only because of the "device" attribute

packloc_custom_cxx.o: packloc_custom_cxx.cpp
	$(CXX) -c $(CXXARCHFLAGS) $(CXXGENERICFLAGS) $(CXXEXTRAFLAGS) $(CXXOPTFLAGS) $<

$(BINARY): packloc_test_fortran.o sum_prefix_custom.o sum_prefix_custom_cxx.o packloc_custom.o packloc_custom_cxx.o

	$(FC) $(LINKARCHFLAGS) $(LINKGENERICFLAGS) $(LINKOPTFLAGS) -o $@ $^ $(LINKERLIBS) $(LINKEXTRALIBS)

sum_prefix_custom.o: sum_prefix_custom_cxx.o Makefile psum.inl.cu
packloc_custom.o: packloc_custom_cxx.o Makefile psum.inl.cu
packloc_test_fortran.o: sum_prefix_custom.o packloc_custom.o Makefile

clean:
	rm -f *.o *.mod $(BINARY)

