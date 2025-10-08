FC = nvfortran

all: tabtest

iceglobals.o: iceglobals.f90
	$(FC) -c -acc=gpu -Minfo=accel -O3 iceglobals.f90

icetab.o: icetab.f90 iceglobals.o
	$(FC) -c -acc=gpu -Minfo=accel -O3 icetab.f90

tabtest.o: tabtest.f90
	$(FC) -c -acc=gpu -Minfo=accel -O3 tabtest.f90

tabtest: icetab.o tabtest.o 
	$(FC) -acc=gpu -Minfo=accel -O3 iceglobals.o icetab.o tabtest.o -o tabtest

run: tabtest
	NVCOMPILER_ACC_NOTIFY=3 ./tabtest

profile: tabtest
	TMPDIR=${PWD} nsys profile --stats=true ./tabtest
	#TMPDIR=${PWD} ncu ./tabtest

clean:
	rm -f *.o *.mod tabtest


