!==================================================================================================
  program init
! Date:  	Nov 15, 2020
! Features:
!
! This is updated 11/2020 based on MathFunda + CG/biEcoBase series + more for SL blends
!		(1) ConstraintRelease parameter (instead of double negative)
!       (2) variable nu
!       (3) Metropolis
!       (4) SL/LL/SS blend
!       (5) Float Z0
!
!
! (*) This module reads from "inp.dat" and prepares an input file "inSnap.dat" for actual dynamics
! (*) This input file "inSnap.dat" can be read using Dynamics.f90
! (*) Both Dynamics.f90 and Initialize.f90 use Utils.f90 for data structures and common functions
!
!==================================================================================================
	use utils
	
	logical :: isLin1, isLin2				! global variables that are defined in "inp.dat"

	integer :: Np1, Np2, nseed				! and read using subroutine "readInput()"
	double precision :: Z1, Z2				! Z1, Z2 float ***new 11/2020***
	
	logical :: ConstraintRelease            ! ***new 11/2020***
	
	type(chain), allocatable :: c(:)		! the datastructure that holds chain-level information
	
	call readInput()			! read in "inp.dat"
	call initSetup()			! allocate arrays and do initial pairing

		
	! randomize pairing - avoid self-entanglement
	if (ConstraintRelease) call randomizePairing()		
	
	call equilibrate(1000)		! equilibrate for nMCS = 1000
	call spitOutput()			! write inSnap.dat
	
	! debugging routine, useful to check if anything amiss in pairing
	if (ConstraintRelease) call debugConsistent()
	
	contains
	
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
		
	!----------------------------------------------
		subroutine initSetup()
	!     allocate arrays and do initial pairing
	!----------------------------------------------
			implicit none

			integer  :: i, j, maxsl, maxslfac, Np
			
			maxslfac = 3
			
			if(modulo(Np1, 2) /= 0 .or. Np1 == 0 .or. modulo(Np2, 2) /= 0) then
				print*, "Need even Np1 > 0; and even Np2 (=0 okay)"
				stop
			endif

			! allocate chain arrays
			Np = Np1 + Np2

			allocate(c(Np))
			
			! allocate/populate the first species (this number Np1 > 0)			
			do i = 1, Np, 2


				if (i <= Np1) then
				
					c(i)%isLin = isLin1
					c(i)%Z  = nint(Z1); c(i)%Z0 = Z1; c(i)%Zinit = nint(Z1)
					c(i)%nu = 1.5 * Z1/(Z1 + 1)

					c(i+1)%Z  = nint(Z1); c(i+1)%Z0 = Z1; c(i+1)%Zinit = nint(Z1)
					c(i+1)%nu = 1.5 * Z1/(Z1 + 1)
	
					maxsl = maxslfac*Z1					! may need to keep track in create

				else

					c(i)%isLin = isLin2				
					c(i)%Z  = nint(Z2); c(i)%Z0 = Z2; c(i)%Zinit = nint(Z2)
					c(i)%nu = 1.5 * Z2/(Z2 + 1)

					c(i+1)%Z  = nint(Z2);	c(i+1)%Z0 = Z2; c(i+1)%Zinit = nint(Z2)
					c(i+1)%nu = 1.5 * Z2/(Z2 + 1)
	
					maxsl = maxslfac*Z2					! may need to keep track in create
				
				endif


				allocate(c(i)%s(maxsl))
				allocate(c(i+1)%s(maxsl))

				! pairing
				if(ConstraintRelease) then

					do j = 1, maxsl
						c(i)%s(j)%pch   = i + 1
						c(i)%s(j)%psl   = j
						c(i)%s(j)%org   = .True.	! placeholders currently
						c(i)%s(j)%olp   = j			! will reassign during dynamics
						
						c(i+1)%s(j)%pch = i
						c(i+1)%s(j)%psl = j
						c(i+1)%s(j)%org = .True.
						c(i+1)%s(j)%olp = j
						
					enddo
				
				else
				
					do j = 1, maxsl
						c(i)%s(j)%pch   = 0
						c(i)%s(j)%psl   = 0
						c(i)%s(j)%org   = .True.	! placeholders currently
						c(i)%s(j)%olp   = j			! will reassign during dynamics
						
						c(i+1)%s(j)%pch = 0
						c(i+1)%s(j)%psl = 0
						c(i+1)%s(j)%org = .True.
						c(i+1)%s(j)%olp = j

					enddo
				
				endif
				
			enddo

		end subroutine initSetup

	!---------------------------------------------------------
		subroutine randomizePairing()
	!     randomize pairing of the slip links (no potential)
	!---------------------------------------------------------
			implicit none

			double precision :: xdum
						
			integer :: i, j, pch, psl, Np
			integer :: iswap1, iswap2, jswap1, jswap2 
			
			
			Np = Np1 + Np2
			
			do i = 1, Np
			
				inner: do j = 1, c(i)%Z
				
					! partner chain: avoid self-entanglement; 9/2018 by weight fraction
					xdum = rand(nseed)*(Z1*Np1 + Z2*Np2)
					if (xdum < Z1*Np1) then
						pch = int(xdum/Z1) + 1
					else
						pch = Np1 + int((xdum - Z1*Np1)/Z2) + 1
					endif
					if (pch == i) cycle inner
					
					! partner slip-link: avoid linking sl on common partner chain
					psl = int(rand(nseed) * c(pch)%Z) + 1
					if (c(i)%s(j)%pch == c(pch)%s(psl)%pch) cycle inner


					! do the actual swap: (i) store partners
					
					iswap1 = c(i)%s(j)%pch
					jswap1 = c(i)%s(j)%psl
					
					iswap2 = c(pch)%s(psl)%pch
					jswap2 = c(pch)%s(psl)%psl

					! (ii) perform swap
					
					c(i)%s(j)%pch = pch
					c(i)%s(j)%psl = psl
					
					c(pch)%s(psl)%pch = i
					c(pch)%s(psl)%psl = j

					! swap the partners
					c(iswap1)%s(jswap1)%pch = iswap2
					c(iswap1)%s(jswap1)%psl = jswap2

					c(iswap2)%s(jswap2)%pch = iswap1
					c(iswap2)%s(jswap2)%psl = jswap1				
					
				enddo inner
				
			enddo
			
		end subroutine randomizePairing
	
	!-------------------------------------------------------
		subroutine equilibrate(nMCS)
	!     perform nMCS creation/destruction moves per chain
	!     uses routines create and destroy
	!--------------------------------------------------------
			implicit none
			integer, intent(in) :: nMCS
			
			integer :: it, ich, Np
			
			Np = Np1 + Np2
			
			do it = 1, nMCS * Np

				ich = int(rand(nseed) * Np) + 1

				if(rand(nseed) < 0.5 .and. c(ich)%Z > 0) then
					call destroy(ich)
				else
					call create(ich)
				endif
				

			enddo			
			
		end subroutine equilibrate


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
		
	!----------------------------------------------
		subroutine spitOutput()
	!     write file 'inSnap.dat'
	!----------------------------------------------
			implicit none
			integer :: i, j
			
			open(2, file='inSnap.dat')

			do i = 1, Np1 + Np2
				write(2,*) c(i)%Z
				do j = 1, c(i)%Z
					write(2,'(i4,3x,i4,3x,i4,3x,i4)') i, j, c(i)%s(j)%pch, c(i)%s(j)%psl 
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
			
!			print*,"debug: pairing check successful!"

		end subroutine debugConsistent

  end program init
