!==================================================================================================
  module utils
! Date: Nov 15, 2020
! This module contains common routine and datatypes

! This is updated 11/2020 based on MathFunda + CG/biEcoBase series + more for SL blends
!		(1) ConstraintRelease parameter (instead of double negative)
!       (2) variable nu
!       (3) Metropolis
!       (4) SL/LL/SS blend
!       (5) Float Z0
!==================================================================================================
	
	type sl 				! slip link type
		integer :: pch		! partner chain
		integer :: psl		! partner sl
		logical :: org		! is sl original?
		integer	:: olp		! old position number
	end type sl


	type chain
		double precision :: nu		! depends on chain length
		logical 		 :: isLin	! is chain linear (assume star-branch otherwise)
		integer			 :: Z		! current number of entanglements
		double precision :: Z0		! equilibrium number of entanglements  **previously int***
		integer			 :: Zinit	! # entanglements at t = 0 (useful for dynamics only)
		
		type(sl), allocatable :: s(:)		
	end type chain

  contains

!==================================================================================================
!   RANDOM NUMBER GENERATOR
!	Using the numerical recipes code to generate double precision
!	random numbers
!==================================================================================================
	FUNCTION rand(idum)
	IMPLICIT NONE
	INTEGER, PARAMETER :: K4B=selected_int_kind(9)
	INTEGER(K4B), INTENT(INOUT) :: idum
	DOUBLE PRECISION :: rand 
	!
	!  Minimal random number generator of Park and Miller combined with a 
	!  Marsaglia shift sequence. Returns a uniform random deviate between 0.0 
	!  and 1.0 (exclusive of the endpoint values). This fully portable, scalar 
	!  generator has the traditional(not Fortran 90) calling sequence with a
	!  random deviate as the returned function value: call with idum a negative 
	!  integer to initialize; thereafter, do not alter idum except to reinitialize.
	!  The period of this generator is about 3.1x10^18.

	INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
	DOUBLE PRECISION, SAVE :: am
	INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
	if (idum <= 0 .or. iy < 0) then 	! Initialize
	am=nearest(1.0,-1.0)/IM
	iy=ior(ieor(888889999,abs(idum)),1)
	ix=ieor(777755555,abs(idum))
	idum=abs(idum)+1
	end if
	ix=ieor(ix,ishft(ix,13)) 	! Marsaglia shift sequence
	ix=ieor(ix,ishft(ix,-17))
	ix=ieor(ix,ishft(ix,5))
	k=iy/IQ
	iy=IA*(iy-k*IQ)-IR*k
	if (iy < 0) iy=iy+IM
	rand=am*ior(iand(IM,ieor(ix,iy)),1) 

	END FUNCTION rand
	
end module utils
