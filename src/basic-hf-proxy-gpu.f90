! Copyright © 2021, UChicago Argonne, LLC, OPEN SOURCE LICENSE
! See (root dir/)LICENSE.TXT 

!  Basic Hartree-Fock Proxy Application,  single-offload GPU version 
!  Written by Graham D. Fletcher, Computational Science Division, 
!  Argonne National Laboratory 
!  Computes the 2-electron (repulsion) term of the Fock operator  
!  for a cluster of helium-like atoms 


 module params  
 real(8)      pi  
 parameter  ( pi = 3.1415926535897931d0 )
 real(8)      sqrpi2
 parameter  ( sqrpi2 = ( pi**( -0.5d0 ) )*2.0d0 )
 real(8)      dtol,          rcut 
!parameter  ( dtol= 1.0d-10, rcut= 1.0d-12 ) 
 parameter  ( dtol= 1.0d-12, rcut= 1.0d-12 ) 
 real(8)      tobohrs 
 parameter  ( tobohrs = 1.889725987722d0 )
 end module params 

 program   basic_hf_proxy   
 use       params 
 implicit  none 
 integer   natom,ngauss 
 real(8),  allocatable ::  xpnt( : )
 real(8),  allocatable ::  coef( : )
 real(8),  allocatable ::  geom( : , : )
 real(8),  allocatable ::  fock( : , : )
 real(8),  allocatable ::  dens( : , : )
 real(8),  allocatable ::  schwarz( : )
 real(8)   erep

 integer    ib,jb,kb,lb,  ij,i,j,kl,k,l,nn 
 real(8)    aij,dij,xij,yij,zij, akl,dkl, aijkl,tt,f0t, eri 
 integer    nnnn,ijkl,n    !  need i8?


! ngauss is the number of gaussian-type functions (GTF) per atom.
! The GTF are specified by their exponents (xpnt) and contracted 
! using the expansion coefficients (coef) to make an orbital which 
! is occupied by a pair of opposite-spin electrons.
! Use ngauss (,xpnt,coef) to adjust the per-atom workload.
! natom is the number of helium-like atoms.
 read *,  ngauss, natom 

 allocate( xpnt( ngauss ) )  
 allocate( coef( ngauss ) )  
 allocate( geom( 3, natom ) )  

 read *,  ( xpnt( i ), coef( i ), i = 1, ngauss ) 
 do  i  =  1,  natom       !  read the geometry 
 read *,  ( geom( j, i ), j = 1, 3 )  
 end do

 allocate( dens( natom, natom ) )  

! Typical Hartree-Fock calculations build the density matrix from 
! a guess wave function (not a significant computational effort). 
! To simulate this just make a fake density here.
 do  i  =  1,  natom  
 do  j  =  1,  natom  
 dens( i, j )  =  0.1d0   
 end do
 dens( i, i )  =  1.0d0   
 end do

! Uncomment the lines below to read a density matrix from a file 
! (then broadcast it) to enable a direct comparison with a code 
! such as GAMESS. 
!do  i  =  1,  natom       !  GAMESS outputs a triangle 
!read *,  ( dens( i, j ), j = 1, i )  
!do  j  =  1,  i-1  
!dens( j, i )  =  dens( i, j )  
!end do
!end do

 
! Along with the density matrix, above, 3*(N^2) storage is typical 
! of real Hartree-Fock calculations. 
 allocate( fock( natom, natom ) )  

! Normalize the primitive GTO weights.
 do i = 1, ngauss 
 coef( i ) = coef( i )*( ( 2.0d0*xpnt( i ) )**0.75d0 )
 end do

!  Scale the geometry to Bohrs for energy calculations in AU.
 do  i  =  1,  natom 
 geom( 1, i )  =  geom( 1, i )*tobohrs   
 geom( 2, i )  =  geom( 2, i )*tobohrs   
 geom( 3, i )  =  geom( 3, i )*tobohrs   
 end do

 do  i  =  1,  natom  
 do  j  =  1,  natom  
 fock( i, j )  =  0.0d0  
 end do
 end do


! compute Schwarz Inequality factors for integral screening 

 nn = ((natom**2)+natom)/2  
 allocate( schwarz( nn ) )  

 ij = 0 
 do  i  =  1,  natom  
 do  j  =  1,  i  
 ij = ij + 1 
 call ssss( i, j, i, j, ngauss, xpnt, coef, geom, eri )
 schwarz( ij )  =  sqrt( abs( eri ) )  
 end do
 end do


! The following loop (expanded to four indices, with permutational
! symmetry) represents the kernel of Hartree-Fock calculations.
! Integrals are screened to avoid small terms.
 nnnn = ((nn**2)+nn)/2

!$OMP TARGET MAP(TO:nnnn,ngauss,geom,xpnt,coef,dens,schwarz) MAP(TOFROM:fock)
!$OMP TEAMS DISTRIBUTE 

 do  ijkl  =  1,  nnnn
! decompose triangular ijkl index into ij>=kl
 ij = sqrt( dble( 2*ijkl ) )
 n = ( ij*ij + ij )/2
 do  while ( n .lt. ijkl )
 ij = ij + 1
 n = ( ij*ij + ij )/2
 end do
 kl  =  ijkl - ( ij*ij - ij )/2
 if ( schwarz( ij )*schwarz( kl ) > dtol ) then      

! decompose triangular ij index into i>=j
 i = sqrt( dble( 2*ij ) )
 n = ( i*i + i )/2
 do  while ( n .lt. ij )
 i = i + 1
 n = ( i*i + i )/2
 end do
 j  =  ij - ( i*i - i )/2
! decompose triangular kl index into k>=l
 k = sqrt( dble( 2*kl ) )
 n = ( k*k + k )/2
 do  while ( n .lt. kl )
 k = k + 1
 n = ( k*k + k )/2
 end do
 l  =  kl - ( k*k - k )/2

 eri  =  0.0d0
 do  ib  =  1,  ngauss
 do  jb  =  1,  ngauss
 aij = 1.0d0/( xpnt( ib ) + xpnt( jb ) )
 dij = coef( ib )*coef( jb )*exp( -xpnt( ib )*xpnt( jb )*aij*  &
     ( ( geom( 1, i ) - geom( 1, j ) )**2   &
     + ( geom( 2, i ) - geom( 2, j ) )**2   &
     + ( geom( 3, i ) - geom( 3, j ) )**2  )  )*( aij**1.5d0 )
 if ( abs( dij ) > dtol ) then
 xij = aij*( xpnt( ib )*geom( 1, i ) + xpnt( jb )*geom( 1, j ) )
 yij = aij*( xpnt( ib )*geom( 2, i ) + xpnt( jb )*geom( 2, j ) )
 zij = aij*( xpnt( ib )*geom( 3, i ) + xpnt( jb )*geom( 3, j ) )
 do  kb  =  1,  ngauss
 do  lb  =  1,  ngauss
 akl = 1.0d0/( xpnt( kb ) + xpnt( lb ) )
 dkl = dij*coef( kb )*coef( lb )*exp( -xpnt( kb )*xpnt( lb )*akl*  &
     ( ( geom( 1, k ) - geom( 1, l ) )**2   &
     + ( geom( 2, k ) - geom( 2, l ) )**2   &
     + ( geom( 3, k ) - geom( 3, l ) )**2  )  )*( akl**1.5d0 )
 if ( abs( dkl ) > dtol ) then
 aijkl = ( xpnt( ib ) + xpnt( jb ) )*( xpnt( kb ) + xpnt( lb ) )  &
       / ( xpnt( ib ) + xpnt( jb )  +  xpnt( kb ) + xpnt( lb ) )
 tt = aijkl*( ( xij -akl*( xpnt( kb )*geom( 1, k ) + xpnt( lb )*geom( 1, l ) ) )**2  &
            + ( yij -akl*( xpnt( kb )*geom( 2, k ) + xpnt( lb )*geom( 2, l ) ) )**2  &
            + ( zij -akl*( xpnt( kb )*geom( 3, k ) + xpnt( lb )*geom( 3, l ) ) )**2  )
 f0t  =  sqrpi2
 if ( tt > rcut )  f0t  =  ( tt**( -0.5d0 ) )*erf( sqrt(tt) )
 eri  =  eri  +  dkl*f0t*sqrt(aijkl)
 end if
 end do  ;  end do
 end if
 end do  ;  end do

 if ( i == j ) eri = eri*0.5d0 
 if ( k == l ) eri = eri*0.5d0 
 if ( i == k .and. j == l ) eri = eri*0.5d0
!$OMP ATOMIC
 fock( i, j )  =  fock( i, j )  +dens( k, l )*eri*4.0d0  
!$OMP ATOMIC
 fock( k, l )  =  fock( k, l )  +dens( i, j )*eri*4.0d0  
!$OMP ATOMIC
 fock( i, k )  =  fock( i, k )  -dens( j, l )*eri  
!$OMP ATOMIC
 fock( i, l )  =  fock( i, l )  -dens( j, k )*eri  
!$OMP ATOMIC
 fock( j, k )  =  fock( j, k )  -dens( i, l )*eri  
!$OMP ATOMIC
 fock( j, l )  =  fock( j, l )  -dens( i, k )*eri  

 end if  
 end do  
!$OMP END TEAMS DISTRIBUTE
!$OMP END TARGET 


! It is customary to symmetrize the Fock matrix,
! but this is not necessary just to compute the 2e- energy. 
!do  i  =  1,  natom  
!do  j  =  1,  i - 1
!fock( i, j ) = ( fock( i, j ) + fock( j, i ) )*0.5d0  
!fock( j, i ) = fock( i, j )
!end do
!end do 
!  Print Fock matrix, not needed for timing, validation only.
!do  i  =  1,  natom  
!print *, ( fock( i, j ), j  =  1,  i  ) 
!end do

!  Trace Fock with the density and print the 2e- energy.
 erep  =  0.0d0
 do  i  =  1,  natom  
 do  j  =  1,  natom 
 erep  =  erep  +  fock( i, j )*dens( i, j ) 
 end do
 end do 
 print *, '2e- energy= ',erep*0.5d0  

 deallocate( schwarz )  
 deallocate( xpnt )  
 deallocate( coef )  
 deallocate( dens )  
 deallocate( fock )  
 deallocate( geom )  
 end program basic_hf_proxy   



 subroutine ssss( i, j, k, l, ngauss, xpnt, coef, geom, eri )
 use        params 
 implicit   none 
 integer    i,j,k,l, ngauss 
 real(8)    xpnt(*), coef(*), geom(3,*), eri 
 integer    ib,jb,kb,lb 
 real(8)    aij,dij,xij,yij,zij, akl,dkl, aijkl,tt,f0t  

!  loops over contracted GTOs  

 eri  =  0.0d0 
 do  ib  =  1,  ngauss  
 do  jb  =  1,  ngauss  
 aij = 1.0d0/( xpnt( ib ) + xpnt( jb ) ) 
 dij = coef( ib )*coef( jb )*exp( -xpnt( ib )*xpnt( jb )*aij*  &  
     ( ( geom( 1, i ) - geom( 1, j ) )**2   &
     + ( geom( 2, i ) - geom( 2, j ) )**2   &
     + ( geom( 3, i ) - geom( 3, j ) )**2  )  )*( aij**1.5d0 )  
 if ( abs( dij ) > dtol ) then      
 xij = aij*( xpnt( ib )*geom( 1, i ) + xpnt( jb )*geom( 1, j ) )  
 yij = aij*( xpnt( ib )*geom( 2, i ) + xpnt( jb )*geom( 2, j ) )  
 zij = aij*( xpnt( ib )*geom( 3, i ) + xpnt( jb )*geom( 3, j ) )  
 do  kb  =  1,  ngauss  
 do  lb  =  1,  ngauss 
 akl = 1.0d0/( xpnt( kb ) + xpnt( lb ) ) 
 dkl = dij*coef( kb )*coef( lb )*exp( -xpnt( kb )*xpnt( lb )*akl*  &  
     ( ( geom( 1, k ) - geom( 1, l ) )**2   &
     + ( geom( 2, k ) - geom( 2, l ) )**2   &
     + ( geom( 3, k ) - geom( 3, l ) )**2  )  )*( akl**1.5d0 )  
 if ( abs( dkl ) > dtol ) then      
 aijkl = ( xpnt( ib ) + xpnt( jb ) )*( xpnt( kb ) + xpnt( lb ) )  &  
       / ( xpnt( ib ) + xpnt( jb )  +  xpnt( kb ) + xpnt( lb ) )  
 tt = aijkl*( ( xij -akl*( xpnt( kb )*geom( 1, k ) + xpnt( lb )*geom( 1, l ) ) )**2  & 
            + ( yij -akl*( xpnt( kb )*geom( 2, k ) + xpnt( lb )*geom( 2, l ) ) )**2  & 
            + ( zij -akl*( xpnt( kb )*geom( 3, k ) + xpnt( lb )*geom( 3, l ) ) )**2  ) 
 f0t  =  sqrpi2 
 if ( tt > rcut )  f0t  =  ( tt**( -0.5d0 ) )*erf( sqrt(tt) ) 
 eri  =  eri  +  dkl*f0t*sqrt(aijkl)  
 end if 
 end do  ;  end do  
 end if  
 end do  ;  end do 
 end subroutine ssss


