!==================================================================================================
  program dyna
!
! Date: November 15, 2020
!
! This is updated 11/2020 based on MathFunda + CG/biEcoBase series + more for SL blends
!		(1) ConstraintRelease parameter (instead of double negative)
!       (2) variable nu
!       (3) Metropolis hardcoded
!       (4) SL blend <----***---
!		(5) Float Z1, Z2 [encoded as c%Z0] actual instantaneous #sl is integer
!
!
!--------------------------------------------------------------------------------------------------
! (*) changing dt definition - currently okay for LL - may have to revisit for stars
!
! Date: August 18, 2016
!
! (*) This module reads inputs from "inp.dat" and "inSnap.dat" (from Initialize.f90)
! (*) Both Dynamics.f90 and Initialize.f90 use Utils.f90 for data structures and common functions
! (*) Carries out actual dynamics. Monitors stress and dielectric relaxation functions
! (*) By default simulation carried out until sigma(t) drops to 0.01 of initial value
! (*) writes final configuration in "outSnap.dat"
!==================================================================================================
	use utils
	
	logical :: isLin1, isLin2, isCR

	integer :: Np1, Np2, nseed		! and read using subroutine "readInput()"
	double precision :: Z1, Z2		! ***new 11/2020***
	
	logical :: ConstraintRelease    ! ***new 11/2020***
		
	double precision :: printFreq, gammaLin, gammaStr, tau0Lin, tau0Str
	
	integer :: initSeg1, initSeg2   ! initial number of segments for linear/star fractions
	
	type(chain), allocatable :: c(:)

	!
	! Qualitative study
	!
	gammaLin  = 1.4 				! formerly known as nLin
	gammaStr  = 3.5 				! formerly known as nStr

	tau0Lin   = 0.12 
	tau0Str   = 1.8e-3

!~ 	!
!~ 	! Pure samples
!~ 	!
!~ 	gammaLin  = 0.0 				! formerly known as nLin
!~ 	gammaStr  = 0.0 				! formerly known as nStr
!~ 	tau0Lin   = 1.0 
!~ 	tau0Str   = 1.0


	! Desai 2016 and Shivokhin 2014
!~ 	tau0Lin   = 0.1768e-6 
!~ 	tau0Str   = 2.754e-9

!~ 	! Struglinski 2016
!~ 	tau0Lin   = 0.1456e-6 
!~ 	tau0Str   = 2.268e-9

	printFreq = 1.1
	
	call readInput()		! read in "inp.dat"
	call initSetupPlus()	! allocate arrays, read from "inSnap.dat", setup initSeg property calculations
	call equilibrate()		! writes (t, phi, psi) to file "relax.dat"
	call spitOutput()		! writes "outSnap.dat" - the final configuration
	
!	call debugConsistent()  ! debug file
	
	contains
	
	!---------------------------------------------------
		subroutine equilibrate()
	!     	perform a creation/destruction moves until
	! 		sigma(t) drops to 0.01 of initial value
	!---------------------------------------------------
			implicit none
			
			integer 		 :: ich, Np
			double precision :: f1, phi, psi, t, dt, tgt, dt1, dt2
			
			open(3, file='relax.dat')
			open(4, file='extra.dat')

			!
			! New vocabulory and setting of timestep involving different prefactors
			! not just the gamma exponents. (11/6/2020)
			!
			if(isLin1) then
				dt1 = (tau0Lin * Z1**gammaLin)/dfloat(Np1)
			else
				dt1 = (tau0Str * Z1**gammaStr)/dfloat(Np1)
			endif

			if(isLin2) then
				dt2 = (tau0Lin * Z2**gammaLin)/dfloat(Np2)
			else
				dt2 = (tau0Str * Z2**gammaStr)/dfloat(Np2)
			endif

				
			Np  = Np1 + Np2
			phi = 1.0        	 		! need to set if making loop dependent on it
			t   = 0.

			dt  = 1./(1.0/dt1 + 1.0/dt2) ! proportional to Zbar**n, to equilibrate slack
			f1  = dt/dt1                        ! finding probability of picking each fraction


			write(*, '(e12.4,1x,e12.4,1x,e12.4, 1x,f7.4)') dt1, dt2, dt, f1
!~ 			print*, "dt1, dt2, dt2/dt1", dt1, dt2, dt2/dt1
!~ 			print*, "dt = ", dt
!~ 			print*, "f1 = ", f1
!~ 			stop

			tgt = 5.*dt
			
!~ 			do while (phi > 0.001)
			do while (phi > 0.0001)

				t = t + dt
								
				! pick polymer based on size. Longer polymer picked less frequently (stress equilibrate?)
				if(rand(nseed) < f1) then
					ich = int(rand(nseed) * Np1) + 1
				else
					ich = Np1 + int(rand(nseed) * Np2) + 1
				endif

				if(rand(nseed) < 0.5 .and. c(ich)%Z > 0) then
					call destroy(ich)
				else
					call create(ich)
				endif
				
				if(t >= tgt) then
					 call calcProperties(phi, psi)
					 write(3,'(E12.5,2X,F8.4,2X,F8.4)') t, phi, psi
					 tgt = printFreq * t
!~ 					 call calcProperties1()
!~ 					 call flush(3)
				endif

				
			enddo			
					
			close(4)	
			close(3)
			
		end subroutine equilibrate

	!-------------------------------------------------
		subroutine calcProperties1()
	!     calculate stress and dielectric relaxation
	!-------------------------------------------------
			implicit none

			
			integer  :: i, j, isum, jsum
			integer  :: numSL11, numSL12, numSL22, numSL

			! Start with First Fraction of Blend
			
			numSL11 = 0
			numSL12 = 0
			numSL22 = 0
			numSL   = 0
			
			do i = 1, Np1 + Np2
				
				numSL = numSL + c(i)%Z
				
				if(i <= Np1) then
					do j = 1, c(i)%Z
						if(c(i)%s(j)%pch <= Np1) then
							numSL11 = numSL11 + 1
						else
							numSL12 = numSL12 + 1
						endif
					enddo
				else
					do j = 1, c(i)%Z
						if(c(i)%s(j)%pch <= Np1) then
							numSL12 = numSL12 + 1
						else
							numSL22 = numSL22 + 1
						endif
					enddo
				endif
													
			enddo

			write(*, '(F8.4, 1x, F8.4, 1x, F8.4)') numSL11/dfloat(numSL), numSL12/dfloat(numSL), numSL22/dfloat(numSL)
			
		end subroutine calcProperties1

	!-------------------------------------------------
		subroutine calcProperties(phi, psi)
	!     calculate stress and dielectric relaxation
	!-------------------------------------------------
			implicit none
			double precision, intent(out)  :: phi, psi

			
			integer  :: i, j, isum, jsum, isumd, jsumd, icoll
			integer  :: olpFirst, olpLast, numSL, nowSeg, nowEE

			integer  :: nowSeg1, initSeg

			! Start with First Fraction of Blend
			
			nowSeg   = 0  ! # surviving segments (=nowSL for stars; different for Linears)
			nowEE    = 0  ! end-to-end vector sum
			initSeg  = initSeg1 + initSeg2
			
			nowSeg1  = 0
			
			
			do i = 1, Np1 + Np2
			
				olpFirst = 0	! loc of first surviving slip-link
				olpLast  = 0	! loc of first surviving slip-link
				numSL    = 0	! # slip links on current chain
				
				! find number of surviving SL and original position of first one
				do j = 1, c(i)%Z
					if(c(i)%s(j)%org) then
						numSL = numSL + 1							! add original SL	
						if(olpFirst == 0) olpFirst = c(i)%s(j)%olp	! loc of first original SL
						if(olpFirst > 0)  olpLast  = c(i)%s(j)%olp  ! loc of last  original SL
					endif
				enddo						

				! Now find contribution to EE based on architecture (star/lin)
				if(c(i)%isLin) then  ! linear 

					! no stress/dielectric contribution if 0 or 1 
					! original SL left on linear
					if(numSL > 1) then
						nowSeg = nowSeg + numSL - 1
						nowEE  = nowEE  + (olpLast - olpFirst)
					endif

				else ! star
					nowSeg = nowSeg + numSL
					nowEE  = nowEE  + olpLast		! if no orig SL, then olpLast = 0
				endif

				if(i == Np1) then
					nowSeg1 = nowSeg
				endif
													
			enddo

 			phi = dfloat(nowSeg)/dfloat(initSeg)
 			psi = dfloat(nowEE)/dfloat(initSeg) 
			
			write(4, '(F8.4, 1x, F8.4, 1x, F8.4)') nowSeg1/dfloat(initSeg1), dfloat(nowSeg-nowSeg1)/initSeg2, phi 
			
		end subroutine calcProperties


	!----------------------------------------------
		subroutine create(i)
	!     create new slip-link on chain i
	!----------------------------------------------
			implicit none
			integer, intent(in) :: i
	
			integer :: j, jp, ip, ipp, jpp, nsl, pch, psl, olp
			logical :: org
			double precision :: dU, xdum, pacc
						
			! choose partner, protect against self-entanglement: Sep 2018 - pick according to wt fraction
			if (ConstraintRelease) then
				ip = i
				do while(ip == i)
					xdum = rand(nseed)*(Z1*Np1 + Z2*Np2)
					if (xdum < Z1*Np1) then
						ip = int(xdum/Z1) + 1
					else
						ip = Np1 + int((xdum - Z1*Np1)/Z2) + 1
					endif
				enddo
			endif
	
			
			! choose end for creation if linear; default is far end
			if(c(i)%isLin .and. rand(nseed) < 0.5) then
				nsl = 1
			else
				nsl = c(i)%Z + 1
			endif

			! check energy allowed?
			dU = c(i)%nu * (+2. * (c(i)%Z - c(i)%Z0) + 1)/c(i)%Z0

			if(c(i)%Z >= 3*c(i)%Z0) then
				print*, "max slip-links breached", i, c(i)%Z, c(i)%Z0
				stop
			endif

			if (ConstraintRelease) then
				dU = dU + c(ip)%nu * (2. * (c(ip)%Z - c(ip)%Z0) + 1)/c(ip)%Z0
				if(c(ip)%Z >= 3*c(ip)%Z0) then
					print*, "max slip-links breached", ip, c(ip)%Z, c(ip)%Z0
					stop
				endif
			endif

			pacc = min(1.0, dexp(-dU))
			
			! if true, then enact move
			if(rand(nseed) < pacc) then
			
				if (ConstraintRelease) then							
					jp = int(rand(nseed)*(c(ip)%Z + 1)) + 1 	! find partner SL on chain ip
				endif
				
				! renumber "test" chain
				if(nsl == 1) then

					! save current sl#1
					pch = c(i)%s(1)%pch
					psl = c(i)%s(1)%psl
					org = c(i)%s(1)%org
					olp = c(i)%s(1)%olp


					! shift all sl back by 1
					do j = c(i)%Z, 1, -1
					
						c(i)%s(j+1)%pch = c(i)%s(j)%pch
						c(i)%s(j+1)%psl = c(i)%s(j)%psl
						c(i)%s(j+1)%org = c(i)%s(j)%org
						c(i)%s(j+1)%olp = c(i)%s(j)%olp

						if (ConstraintRelease) then
							ipp = c(i)%s(j+1)%pch
							jpp = c(i)%s(j+1)%psl
							
							c(ipp)%s(jpp)%pch = i
							c(ipp)%s(jpp)%psl = j+1
						endif

					enddo

				endif
				
				! just add to the appropriate end
				if (ConstraintRelease) then
					c(i)%s(nsl)%pch = ip
					c(i)%s(nsl)%psl = jp
					c(i)%s(nsl)%org = .false.
					c(i)%s(nsl)%olp = 0
				else
					c(i)%s(nsl)%pch = 0
					c(i)%s(nsl)%psl = 0
					c(i)%s(nsl)%org = .false.
					c(i)%s(nsl)%olp = 0
				endif

				c(i)%Z  = c(i)%Z + 1					
							
				! renumber entangled chain
				if (ConstraintRelease) then
					! save current sl (ip,jp)
					pch = c(ip)%s(jp)%pch
					psl = c(ip)%s(jp)%psl
					org = c(ip)%s(jp)%org
					olp = c(ip)%s(jp)%olp

					if(jp < c(ip)%Z + 1) then
					
						! shift post jp sl by 1
						do j = c(ip)%Z, jp, -1

							c(ip)%s(j+1)%pch = c(ip)%s(j)%pch
							c(ip)%s(j+1)%psl = c(ip)%s(j)%psl
							c(ip)%s(j+1)%org = c(ip)%s(j)%org
							c(ip)%s(j+1)%olp = c(ip)%s(j)%olp

							ipp = c(ip)%s(j+1)%pch
							jpp = c(ip)%s(j+1)%psl
							
							c(ipp)%s(jpp)%pch = ip
							c(ipp)%s(jpp)%psl = j + 1

						enddo				
						
					endif  
					
					! link sl1 to ip, jp
					c(ip)%s(jp)%pch = i
					c(ip)%s(jp)%psl = nsl
					c(ip)%s(jp)%org = .false.
					c(ip)%s(jp)%olp = 0	

					c(ip)%Z = c(ip)%Z + 1
				endif

			endif

		end subroutine create
		
	!----------------------------------------------
		subroutine destroy(i)
	!     destroy SL on chain i
	!----------------------------------------------
			implicit none
			integer, intent(in) :: i		
			
			integer :: nsl, j, ip, jp, ipp, jpp
			double precision :: dU, pacc
			
			! choose end for destruction if linear; default is far end
						
			if(c(i)%isLin .and. rand(nseed) < 0.5) then
				nsl = 1
			else
				nsl = c(i)%Z
			endif				

			! identify partner chain
			if (ConstraintRelease) then
				ip = c(i)%s(nsl)%pch
				jp = c(i)%s(nsl)%psl
			endif
			
			! check energy change for proposed move: disallow?
			dU = c(i)%nu * (-2. * (c(i)%Z - c(i)%Z0) + 1)/c(i)%Z0

			if (ConstraintRelease) then
				dU = dU + c(ip)%nu * (-2. * (c(ip)%Z - c(ip)%Z0) + 1)/c(ip)%Z0
			endif			

			! acceptance or rejection
			pacc = min(1.0, dexp(-dU))			
				
			! if true, then enact move
			if(rand(nseed) < pacc) then
			
				! destroy partner chain sl: renumber if bead in internal
				if (ConstraintRelease) then
					if(jp /= c(ip)%Z) then
					
						do j = jp + 1, c(ip)%Z
						
							c(ip)%s(j-1)%pch = c(ip)%s(j)%pch
							c(ip)%s(j-1)%psl = c(ip)%s(j)%psl
							c(ip)%s(j-1)%org = c(ip)%s(j)%org
							c(ip)%s(j-1)%olp = c(ip)%s(j)%olp

							if(j /= c(ip)%Z .or. ip /= i) then						
								ipp = c(ip)%s(j-1)%pch
								jpp = c(ip)%s(j-1)%psl
								
								c(ipp)%s(jpp)%pch = ip
								c(ipp)%s(jpp)%psl = j-1
							endif
						enddo
					endif
				endif

				! destroy slip link on current chain; renumber if end "1" chosen
				
				if(nsl == 1) then
				
					do j = 2, c(i)%Z
					
						c(i)%s(j-1)%pch = c(i)%s(j)%pch
						c(i)%s(j-1)%psl = c(i)%s(j)%psl
						c(i)%s(j-1)%org = c(i)%s(j)%org
						c(i)%s(j-1)%olp = c(i)%s(j)%olp

						if (ConstraintRelease) then
						
							ipp = c(i)%s(j-1)%pch
							jpp = c(i)%s(j-1)%psl
							
							c(ipp)%s(jpp)%pch = i
							c(ipp)%s(jpp)%psl = j-1
							
						endif

					enddo
				endif

				c(i)%Z  = c(i)%Z - 1
				if (ConstraintRelease) c(ip)%Z = c(ip)%Z - 1
				
			endif
			
		end subroutine destroy

	!-------------------------------
		subroutine readInput()
	!     read input data "inp.dat"
	!-------------------------------
			implicit none
			character(len=25) :: cdum	
			
			! read the input file
			open(1,file='inp.dat')
			
			read(1,*) cdum
			read(1,*) isLin1, isLin2
			read(1,*) cdum
			read(1,*) Np1, Np2
			read(1,*) cdum
			read(1,*) Z1, Z2
			read(1,*) cdum
			read(1,*) nseed
			read(1,*) cdum
			read(1,*) ConstraintRelease
			
			close(1)

		end subroutine readInput

	!----------------------------------------------------------------------
		subroutine initSetupPlus()
	!     allocate arrays and read from "inSnap.dat"
	!     also set up initSeg for use in stress/dielectric calculation
	!----------------------------------------------------------------------
			implicit none

			integer  :: i, j, maxsl, maxslfac, Np, junk1, junk2
			
			maxslfac = 3
			
			! allocate chain arrays
			Np = Np1 + Np2

			allocate(c(Np))
			
			initSeg1 = 0
			initSeg2 = 0
			
			open(2, file='inSnap.dat')

			! allocate/populate the first species (this number Np1 > 0)			
			do i = 1, Np

				if (i <= Np1) then
					c(i)%isLin = isLin1
					c(i)%nu = 1.5 * Z1/(Z1 + 1)
					c(i)%Z0 = Z1
					maxsl = maxslfac*Z1					! may need to keep track in create
				else
					c(i)%isLin = isLin2
					c(i)%nu = 1.5 * Z2/(Z2 + 1)
					c(i)%Z0 = Z2
					maxsl = maxslfac*Z2					
				endif

				allocate(c(i)%s(maxsl))


				read(2,*) c(i)%Z
				c(i)%Zinit = c(i)%Z
				
				! keep track of initial segments
				if(c(i)%isLin) then

					if(i <= Np1) then
						initSeg1 = initSeg1 + c(i)%Z - 1	! correction for linear, segments = #SL - 1
					else
						initSeg2 = initSeg2 + c(i)%Z - 1
					endif

				else

					if(i <= Np1) then
						initSeg1 = initSeg1 + c(i)%Z	! correction for linear, segments = #SL - 1
					else
						initSeg2 = initSeg2 + c(i)%Z
					endif

				endif
				
				do j = 1, c(i)%Z
				
					read(2,'(i4,3x,i4,3x,i4,3x,i4)') junk1,junk2, c(i)%s(j)%pch, c(i)%s(j)%psl

					c(i)%s(j)%org = .true.
					c(i)%s(j)%olp = j

				enddo

			enddo   

		end subroutine initSetupPlus
		
	!----------------------------------------------
		subroutine spitOutput()
	!     write file 'outSnap.dat'
	!----------------------------------------------
			implicit none
			integer :: i, j
			
			open(2, file='outSnap.dat')

			do i = 1, Np1 + Np2
				write(2,*) c(i)%Z
				do j = 1, c(i)%Z
					write(2,'(i4,3x,i4,3x,i4,3x,i4,2x,L1,2x,i4)') i, j, c(i)%s(j)%pch, &
					 c(i)%s(j)%psl, c(i)%s(j)%org, c(i)%s(j)%olp
				enddo
			enddo
			
			close(2)

		end subroutine spitOutput

	!----------------------------------------------
		subroutine debugConsistent()
	!     check if pairs are correct
	!----------------------------------------------
			implicit none
			integer :: i, j, ip, jp
			

			do i = 1, Np1 + Np2
				do j = 1, c(i)%Z

					ip = c(i)%s(j)%pch
					jp = c(i)%s(j)%psl
					
					if(c(ip)%s(jp)%pch /= i .or. c(ip)%s(jp)%psl /= j) then
						print*, "consistency problem", i, j, ip, jp
						stop
					endif
					
					if(i == ip) then
						print*, "self-entanglement problem"
						stop
					endif
					
				enddo
			enddo
			
			print*,"debug: pairing check successful!"

		end subroutine debugConsistent

  end program dyna
