#    Basic Hartree-Fock Proxy Application 

 FCSEQ=gfortran
#FCSEQ=ifort
 FCMPI=mpifort
#FCMPI=mpiifort
#FLAGS=
 FLAGS=-O3 -w
#OMP=
 OMP=-fopenmp
#OMP=-qopenmp

all: bhfp-seq bhfp-mpi bhfp-omp bhfp-gpu bhfp-mpi+omp bhfp-mpi+gpu 

bhfp-seq:  src/basic-hf-proxy-seq.f90 
	$(FCSEQ) $(FLAGS) -o bin/bhfp-seq src/basic-hf-proxy-seq.f90 

bhfp-mpi:  src/basic-hf-proxy-mpi.f90 
	$(FCMPI) $(FLAGS) -o bin/bhfp-mpi src/basic-hf-proxy-mpi.f90 

bhfp-omp:  src/basic-hf-proxy-omp.f90 
#	$(FCSEQ) $(OMP) $(FLAGS) -i8 -o bin/bhfp-omp src/basic-hf-proxy-omp.f90 
	$(FCSEQ) $(OMP) $(FLAGS) -o bin/bhfp-omp src/basic-hf-proxy-omp.f90 

bhfp-gpu:  src/basic-hf-proxy-gpu.f90 
	$(FCSEQ) $(OMP) $(FLAGS) -o bin/bhfp-gpu src/basic-hf-proxy-gpu.f90 

bhfp-mpi+omp:  src/basic-hf-proxy-mpi+omp.f90 
	$(FCMPI) $(OMP) $(FLAGS) -o bin/bhfp-mpi+omp src/basic-hf-proxy-mpi+omp.f90 

bhfp-mpi+gpu:  src/basic-hf-proxy-mpi+gpu.f90 
	$(FCMPI) $(OMP) $(FLAGS) -o bin/bhfp-mpi+gpu src/basic-hf-proxy-mpi+gpu.f90 

clean:
	rm params.mod

