FC = nvfortran
CXX = nvcc

all: sc_test

sc_test.o: sc_test.f90
	$(FC) -c -o $@ -Mpreprocess -acc=gpu -gpu=cc90 -cuda -Minfo=accel -O3 $^

sum_prefix_custom.o: sum_prefix_custom.f90
	$(FC) -c -o $@ -Mpreprocess -acc=gpu -gpu=cc90 -cuda -Minfo=accel -O3 $^

thrust_wrapper.o: thrust_wrapper.cpp
	$(CXX) -c -gencode arch=compute_90,code=sm_90 --x cu -O3 $^

#sc_test:  thrust_wrapper.o sc_test.o sum_prefix_custom.o 
sc_test:  sc_test.o sum_prefix_custom.o thrust_wrapper.o
	$(FC) -acc=gpu -gpu=cc90 -o $@ $^ -lstdc++ -L$(NVHPC_ROOT)/cuda/lib -lcudart -cudalib=cutensor

## what really works:
## nvfortran -O3 -cuda -Mpreprocess -acc=gpu -gpu=cc90 -o sc_test sc_test.f90 sum_prefix_custom.f90 thrust_wrapper.o -lstdc++ -L/apps/ACC/NVIDIA-HPC-SDK/25.7/Linux_x86_64/25.7/cuda/lib -lcudart -cudalib=cutensor

clean:
	rm -f *.o *.mod sc_test

