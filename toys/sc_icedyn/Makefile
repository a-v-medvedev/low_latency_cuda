FC = nvfortran
CXX = nvcc

all: sc_test

sc_test.o: sc_test.f90
	$(FC) -c -o $@ -Mpreprocess -acc=gpu -gpu=cc90 -cudalib=cutensor -Minfo=accel -O3 $^

# FIXME thrust_wrapper.o and sum_prefix_custom.o are temporarely switched off
#sum_prefix_custom.o: sum_prefix_custom.f90
#	$(FC) -c -o $@ -Mpreprocess -acc=gpu -gpu=cc90 -Minfo=accel -O3 $^
#
#thrust_wrapper.o: thrust_wrapper.cpp
#	$(CXX) -c -gencode arch=compute_80,code=sm_80 --x cu -O3 $^

# FIXME thrust_wrapper.o and sum_prefix_custom.o are temporarely switched off
#sc_test:  thrust_wrapper.o sc_test.o sum_prefix_custom.o 
sc_test:  sc_test.o #sum_prefix_custom.o thrust_wrapper.o
	$(FC) -acc=gpu -gpu=cc90 -o $@ $^ -lstdc++ -L$(NVHPC_ROOT)/cuda/lib -lcudart -cudalib=cutensor

## does this really work on mn5-acc? we tend to have reports on device code incompatibilities
## nvfortran -O3 -cuda -Mpreprocess -acc=gpu -gpu=cc90 -o sc_test sc_test.f90 sum_prefix_custom.f90 thrust_wrapper.o -lstdc++ -L/apps/ACC/NVIDIA-HPC-SDK/25.7/Linux_x86_64/25.7/cuda/lib -lcudart -cudalib=cutensor

clean:
	rm -f *.o *.mod sc_test

