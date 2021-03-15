# Basic Hartree-Fock Proxy Application 

Software language: Fortran (90+)

### Overview 

Hartree-Fock is considered the 'entry-level' method of quantum chemistry. 
Its computational characteristics are similar to those of Density Functional Theory. 
As with many quantum chemistry methods, the implementation typically depends on a choice for the integral computation method. 
There are many integral methods and each has its particular performance characteristics. 
This proxy application aims to simulate the compute load and data-movement of the kernel of Hartree-Fock in a way that is independent of the particular integral method. 
GPU offloading is also a key aim.

### Objective

The kernel of Hartree-Fock for a given molecular system computes the electron-repulsion term of the Fock operator. 
Electron repulsion integrals over gaussian-type functions located on different atomic centers (four per integral) are computed and contracted with elements of the (square) density matrix before being summed to elements of the (square) Fock matrix. 
The cost increases as the fourth power of the problem size. 
However, the magnitude of an integral declines with distance between the centers and screening (based on the Schwarz Inequality) can be used to avoid small terms. 
Furthermore, within the integral calculations the gaussian-product factors can also be tested to avoid small terms. 
While the number of atoms is a primary cost factor, it is through such screening that the atomic positions also influence the overall cost.

### Model

To simulate the compute load, a cluster of helium-like atoms each with a single contracted (s-type) gaussian type function (orbital) occupied by a pair of opposite-spin electrons, is used. 
The function cost can be varied by changing the expansion length while the number and positions of the atoms can also be varied (see below).
As mentioned above, the compute load increases as N<sup>4</sup>, where N is the number of helium-like atoms. 
Screening typically reduces the cost to O(N<sup>3</sup>). 
The individual integral cost is, likewise, fourth-order in the orbital expansion length(s) but very long expansions (>>30) are not typical of routine calculations (see below). 
Memory requirements and data-movement scale as O(N<sup>2</sup>). 


### Variables

In order to simulate different compute loads, the following are variable- 
* Number of helium-like atoms ('natom' in the code) 
Determines the memory requirements (order of the Fock and Density matrices), the data-movement and the primary compute loads. 
* Contraction length of atomic orbital ('ngauss' in the code)
Determines the per-atom compute load. 
The gaussian exponents and contraction coefficients can be entirely fabricated, but taking values from a real basis set yields a more realistic calculation. 
* Arrangement of atoms (the input geometry)
This varies the compute load via integral screening with closer spacing giving higher complexity. 


### Current Test Set Input Choices
Helium atoms are arranged on a regular 3D grid.
The He...He distance is chosen to be 1.4 A, approximating the C-C bond length, and yielding a cost complexity of around 2.6. 
A separation of 1.0 A would be closer to bond-to-hydrogen (OH,CH,NH...).
1.2A represents a compromise, yielding a complexity of around 2.9. 


### Software Overview 

| Proxy App Versions         | Description  | 
|  :-:        |     :-:                |
| basic-hf-proxy-seq.f90     | Plain sequential version  | 
| basic-hf-proxy-mpi.f90     | MPI parallel version  | 
| basic-hf-proxy-omp.f90     | Multi-threaded parallel, using OpenMP  | 
| basic-hf-proxy-gpu.f90     | Single GPU offload, using OMP5.0  | 
| basic-hf-proxy-mpi+omp.f90 | Combining MPI for outer loops, OMP for inner loops  | 
| basic-hf-proxy-mpi+gpu.f90 | MPI, with GPU offload for inner loops  | 



## Build instructions

| Code                       | Dependencies  | 
|  :-:        |     :-:                |
| basic-hf-proxy-seq.f90     | Fortran 90+ compiler  | 
| basic-hf-proxy-mpi.f90     | Fortran 90+ compiler, MPI library  | 
| basic-hf-proxy-omp.f90     | Fortran 90+ compiler supporting OpenMP  | 
| basic-hf-proxy-gpu.f90     | Fortran 90+ compiler supporting OpenMP 5.0 or higher  | 
| basic-hf-proxy-mpi+omp.f90 | Fortran 90+ compiler supporting OpenMP, MPI library  | 
| basic-hf-proxy-mpi+gpu.f90 | Fortran 90+ compiler supporting OpenMP 5.0 or higher, MPI library  | 
 
Note that the 'mpi+gpu' version has not been tested yet (TBD).  

In the root dir, 

$ mkdir bin 

$ make 


### Specific build examples 
 
Sequential version, with GNU fortran compiler
 
$ gfortran basic-hf-proxy-seq.f90

MPI parallel version, with compiler-wrapper

$ mpifort basic-hf-proxy-mpi.f90

MPI parallel version, with MPI library

$ gfortran basic-hf-proxy-mpi.f90 -L[path to MPI lib] -I[path to mpif.h] 

OpenMP parallel version

$ gfortran -fopenmp basic-hf-proxy-omp.f90

For GPU offload, eg. with pgf90 (supporting OMP5.0+) 

$ pgf90 -mp=gpu -Mcuda=cc80,cuda11.2 basic-hf-proxy-gpu.f90

See also Nsight Compute profile, ./docs/ncu.profile.bhfp-gpu.he64 



## Run Examples (see ./tests) 

Sequential version

$ ./a.out < he4

MPI version

$ mpirun -n 2 ./a.out < he4

OMP version

$ export OMP_NUM_THREADS=2

$ ./a.out < he8

#### Small regression tests  (see ./tests) 
|no. He atoms | 2-el. E/AU to 6 places |
|  :-:        |     :-:                |
|   4 |    4.050176 | 
|   8 |   11.585213 | 
|  16 |   36.902417 | 
|  32 |  122.248336 | 

#### Large regression tests  

|no. He atoms | 2-el. E/AU    | 
|  :-:        |     :-:                |
|  64 |  419.672121   | 
| 128 | 1300.093626   | 
| 256 | 4204.804478   | 
| 512 |14131.752017   | 

Agreement across different hardware, compilers, etc, to 6 places is typical.


## Figure-of-Merit (FOM)

Sample performance data with various hardware and execution modes (default optimization levels). 
The purpose here is not to present optimal performance figures but rather to suggest how improved might be improved. 


#### GPU Speedup on Nvidia A100 

<img src="./docs/gpuSpeedup.png" height="350"/> 

See also Nsight Compute profile, ./docs/ncu.profile.bhfp-gpu.he64 

#### MPI with OpenMP multi-threading (on Skylake node) 

<img src="./docs/mpi+ompSkylake.png" height="350"/> 

#### MPI parallelism on KNL cluster (ALCF "Theta" platform) 

<img src="./docs/mpiThetaKNL.png" height="350"/> 


### How to Validate
#### Example using GAMESS

* 1. See comments in the proxy app code about reading a density matrix from a file, rebuild the proxy executable as appropriate.

* 2. Convert the input to a GAMESS .inp file, see examples (he4.inp, he8.inp), be sure to include  NPRINT=5  so that GAMESS prints the density matrix

* 3. Run GAMESS.

* 4. Extract the final density matrix from the GAMESS output and convert it to a simple triangular format, appending it to the input file for the proxy app.

* 5. Run the new proxy executable with the new input, and compare output to the line matching-
  'ELECTRON-ELECTRON POTENTIAL ENERGY'
   in the GAMESS output.


