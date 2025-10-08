FC = nvfortran
CXX = nvcc

all: sc_test

sc_test.o: sc_test.f90
	$(FC) -c -Mpreprocess -acc=gpu -gpu=cc90 -cuda -Minfo=accel -O3 $^

sum_prefix_custom.o: sum_prefix_custom.f90
	$(FC) -c -Mpreprocess -acc=gpu -gpu=cc90 -cuda -Minfo=accel -O3 $^

thrust_wrapper.o: thrust_wrapper.cpp
	$(CXX) -c -gencode arch=compute_90,code=sm_90 --x cu -O3 $^

sc_test:  thrust_wrapper.o sc_test.o sum_prefix_custom.o 
	$(FC) -acc=gpu -gpu=cc90 -o $@ $^ -cudalib=cutensor -lstdc++ -L$(NVHPC_ROOT)/cuda/lib64 -lcudart

clean:
	rm -f *.o *.mod sc_test

