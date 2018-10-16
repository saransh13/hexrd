! ###################################################################
! Copyright (c) 2013-2014, Marc De Graef/Carnegie Mellon University
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
!     - Redistributions of source code must retain the above copyright notice, this list 
!        of conditions and the following disclaimer.
!     - Redistributions in binary form must reproduce the above copyright notice, this 
!        list of conditions and the following disclaimer in the documentation and/or 
!        other materials provided with the distribution.
!     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names 
!        of its contributors may be used to endorse or promote products derived from 
!        this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
! USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ###################################################################

!--------------------------------------------------------------------------
! EMsoft:so3.f90
!--------------------------------------------------------------------------
!
! MODULE: so3
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief everything that has to do with sampling of rotation space SO(3)
!
!> @todo verify that this is correct for point groups with multiple settings, eg, 3m, 32, ...
!
!> @date 05/29/14 MDG 1.0 original
!> @date 10/02/14 MDG 2.0 removed globals + rewrite
!> @date 01/01/15 MDG 2.1 added IsinsideFZ function, also used in dictionary indexing approach
!> @date 01/17/15 MDG 2.2 added gridtype option to SampleRFZ
!> @date 03/03/16 MDG 2.3 added uniform sampling for constant misorientations
!> @date 04/01/17 MDG 3.0 start to add FZs for two-phase systems (e.g., cubic-hexagonal, etc.)
!--------------------------------------------------------------------------
module so3

use local
use constants

real(kind=dbl)          :: LPsPi=3.141592653589793D0       !  pi
real(kind=dbl)          :: LPsiPi=0.318309886183791D0      !  1/pi
real(kind=dbl)          :: LPssPi=1.772453850905516D0      !  sqrt(pi)
real(kind=dbl)          :: LPssPio2=1.253314137315500D0    !  sqrt(pi/2)
real(kind=dbl)          :: LPssPi2=0.886226925452758D0     !  sqrt(pi)/2
real(kind=dbl)          :: LPssrt=0.86602540378D0      !  sqrt(3)/2
real(kind=dbl)          :: LPsisrt=0.57735026919D0    !  1/sqrt(3)
real(kind=dbl)          :: LPsalpha=1.346773687088598D0   !  sqrt(pi)/3^(1/4)
real(kind=dbl)          :: LPsrtt=1.732050807568877D0      !  sqrt(3)
real(kind=dbl)          :: LPsprea=0.525037567904332D0    !  3^(1/4)/sqrt(2pi)
real(kind=dbl)          :: LPspreb=1.050075135808664D0     !  3^(1/4)sqrt(2/pi)
real(kind=dbl)          :: LPsprec=0.906899682117109D0    !  pi/2sqrt(3)
real(kind=dbl)          :: LPspred=2.094395102393195D0     !  2pi/3
real(kind=dbl)          :: LPspree=0.759835685651593D0     !  3^(-1/4)
real(kind=dbl)          :: LPspref=1.381976597885342D0     !  sqrt(6/pi)
real(kind=dbl)          :: LPspreg=1.5551203015562141D0    ! 2sqrt(pi)/3^(3/4)
! the following constants are used for the cube to quaternion hemisphere mapping
real(kind=dbl)          :: LPsa=1.925749019958253D0        ! pi^(5/6)/6^(1/6)
real(kind=dbl)          :: LPsap=2.145029397111025D0       ! pi^(2/3)
real(kind=dbl)          :: LPssc=0.897772786961286D0       ! a/ap
real(kind=dbl)          :: LPsbeta=0.962874509979126D0     ! pi^(5/6)/6^(1/6)/2
real(kind=dbl)          :: LPsR1=1.330670039491469D0       ! (3pi/4)^(1/3)
real(kind=dbl)          :: LPsr2=1.414213562373095D0       ! sqrt(2)
real(kind=dbl)          :: LPsr22=0.707106781186547D0      ! 1/sqrt(2)
real(kind=dbl)          :: LPspi12=0.261799387799149D0     ! pi/12
real(kind=dbl)          :: LPspi8=0.392699081698724D0      ! pi/8
real(kind=dbl)          :: LPsprek=1.643456402972504D0     ! R1 2^(1/4)/beta
real(kind=dbl)          :: LPsr24=4.898979485566356D0      ! sqrt(24)

real(kind=dbl)          :: BP(6)= (/ 0.D0, 1.D0, 0.577350269189626D0, 0.414213562373095D0, 0.D0,  &
                                             0.267949192431123D0 /)       ! used for Fundamental Zone determination in so3 module

integer(kind=irg), parameter        :: FZtypeTable(32,32) = reshape( (/ &
 0, 0, 0, 0, 0, 9, 0, 9, 0, 0, 0, 7, 0, 9, 7, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 4, 4, 3, 4, 3, &
 0, 0, 0, 0, 0, 9, 0, 9, 0, 0, 0, 7, 0, 9, 7, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 4, 4, 3, 4, 3, &
 0, 0, 0, 0, 0, 9, 9, 9, 7, 9, 7, 7, 7, 7, 7, 8, 8, 6, 8, 6, 6, 8, 6, 6, 6, 6, 6, 4, 4, 3, 4, 3, &
 0, 0, 0, 0, 0, 9, 0, 9, 0, 0, 0, 7, 0, 9, 7, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 4, 4, 3, 4, 3, &
 0, 0, 0, 0, 0, 9, 9, 9, 7, 9, 7, 7, 7, 7, 7, 8, 8, 6, 8, 6, 6, 8, 6, 6, 6, 6, 6, 4, 4, 3, 4, 3, &
 9, 9, 9, 9, 9, 9, 9, 9, 7, 9, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 3, 4, 3, &
 0, 0, 9, 0, 9, 9, 0, 9, 0, 0, 0, 7, 0, 9, 7, 0, 0, 6, 0, 6, 0, 0, 0, 6, 0, 6, 6, 4, 4, 3, 4, 3, &
 9, 9, 9, 9, 9, 9, 9, 9, 7, 9, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 3, 4, 3, &
 0, 0, 7, 0, 7, 7, 0, 7, 0, 0, 0, 7, 0, 7, 7, 0, 0, 5, 0, 5, 0, 0, 0, 5, 0, 5, 5, 3, 3, 3, 3, 3, &
 0, 0, 9, 0, 9, 9, 0, 9, 0, 0, 0, 7, 0, 9, 7, 0, 0, 6, 0, 6, 0, 0, 0, 6, 0, 6, 6, 4, 4, 3, 4, 3, &
 0, 0, 7, 0, 7, 7, 0, 7, 0, 0, 0, 7, 0, 7, 7, 0, 0, 5, 0, 5, 0, 0, 0, 5, 0, 5, 5, 3, 3, 3, 3, 3, &
 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 3, 3, 3, &
 0, 0, 7, 0, 7, 7, 0, 7, 0, 0, 0, 7, 0, 7, 7, 0, 0, 5, 0, 5, 0, 0, 0, 5, 0, 5, 5, 3, 3, 3, 3, 3, &
 9, 9, 7, 9, 7, 7, 9, 7, 7, 9, 7, 7, 7, 9, 7, 6, 6, 5, 6, 5, 6, 6, 6, 5, 6, 5, 5, 3, 3, 3, 3, 3, &
 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 3, 3, 3, &
 0, 0, 8, 0, 8, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 2, 2, 1, 2, 1, &
 0, 0, 8, 0, 8, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 2, 2, 1, 2, 1, &
 8, 8, 6, 8, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 6, 8, 6, 6, 6, 8, 6, 2, 2, 1, 2, 1, &
 0, 0, 8, 0, 8, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 8, 6, 2, 2, 1, 2, 1, &
 8, 8, 6, 8, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 6, 8, 6, 6, 6, 8, 6, 2, 2, 1, 2, 1, &
 0, 0, 6, 0, 6, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 6, 0, 6, 0, 0, 0, 6, 0, 6, 6, 2, 2, 1, 2, 1, &
 0, 0, 8, 0, 8, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 8, 0, 8, 0, 0, 0, 6, 0, 6, 6, 2, 2, 1, 2, 1, &
 0, 0, 6, 0, 6, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 6, 0, 6, 0, 0, 0, 6, 0, 6, 6, 2, 2, 1, 2, 1, &
 6, 6, 6, 6, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 1, 2, 1, &
 0, 0, 6, 0, 6, 6, 0, 6, 0, 0, 0, 5, 0, 6, 5, 0, 0, 6, 0, 6, 0, 0, 0, 6, 0, 6, 6, 2, 2, 1, 2, 1, &
 8, 8, 6, 8, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 6, 6, 6, 6, 6, 8, 6, 2, 2, 1, 2, 1, &
 6, 6, 6, 6, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 2, 2, 1, 2, 1, &
 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 3, 4, 3, &
 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 3, 4, 3, &
 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, &
 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 3, 4, 3, &
 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3 &
 /), (/ 32, 32/) )

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Routine to return the FZtype and FZorder parameters for single or two-phase
! fundamental zone (FZ) computations; this includes all the FZ types from the 
! following paper:
!
! "Representation of Orientation and Disorientation data for Cubic, Hexagonal, 
! Tetragonal, and Orthorhombic Crystals", A. Heinz and P. Neumann, Acta Cryst. A47, 
! 780-789 (1991)
!
! this routine also allows for icosahedral symmetry, although this is not part 
! of the paper above.
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! SUBROUTINE: getFZtypeandorder
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief does Rodrigues point lie inside the relevant FZ
!
!> @param pgnum1 point group number for phase 1
!> @param FZtype FZ type
!> @param FZorder FZ order
!> @param pgnum2 point group number for phase 2 (optional)
!
!> @date 04/02/17 MDG 1.0 new routine, needed for two-phase disorientations
!--------------------------------------------------------------------------
recursive subroutine getFZtypeandorder(pgnum1,FZtype,FZorder)
!DEC$ ATTRIBUTES DLLEXPORT :: getFZtypeandorder

use typedefs

IMPLICIT NONE

integer(kind=irg),INTENT(IN)              :: pgnum1
integer(kind=irg),INTENT(OUT)             :: FZtype
integer(kind=irg),INTENT(OUT)             :: FZorder

integer(kind=irg)                         :: thisFZType
logical                                   :: twophase

! 0 -> x  (no symmetry or unbounded FZ for the cyclic symmetries)
! 1 -> a  mixed cubic-hexagonal FZ  -> FZtype = 6
! 2 -> b  mixed FZ -> FZtype = 7
! 3 -> c  octahedral FZ
! 4 -> d  tetrahedral FZ
! 5 -> e  24-sided prismatic FZ -> FZtype = 8
! 6 -> f  622 hexagonal dihedral FZ
! 7 -> g  422 octagonal dihedral FZ
! 8 -> h  32 trigonal dihedral FZ
! 9 -> i  222 dihedral FZ

! we reserve FZtype = 5 for icosahedral symmetry, which is under development on a separate code branch

FZtype = FZtarray(pgnum1)
FZorder = FZoarray(pgnum1)


end subroutine getFZtypeandorder

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! We define a number of logical routines, that decide whether or not 
! a point in Rodrigues representation lies inside the fundamental zone (FZ)
! for a given crystal symmetry. This follows the Morawiec@Field paper:
!
! A. Morawiec & D. P. Field (1996) Rodrigues parameterization for orientation 
! and misorientation distributions, Philosophical Magazine A, 73:4, 1113-1130, 
! DOI: 10.1080/01418619608243708
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! FUNCTION: IsinsideFZ
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief does Rodrigues point lie inside the relevant FZ
!
!> @param rod Rodrigues coordinates  (double precision)
!> @param FZtype FZ type
!> @param FZorder FZ order
!
!> @date 01/01/15 MDG 1.0 new routine, needed for dictionary indexing approach
!> @date 06/04/15 MDG 1.1 corrected infty to inftyd (double precision infinity)
!> @date 04/02/17 MDG 1.2 expanded FZ types to include misorientation FZs and icosahedral
!--------------------------------------------------------------------------
recursive function IsinsideFZ(rod,FZtype,FZorder) result(insideFZ)
!DEC$ ATTRIBUTES DLLEXPORT :: IsinsideFZ

use math

IMPLICIT NONE

real(kind=dbl), INTENT(IN)              :: rod(4)
integer(kind=irg),INTENT(IN)            :: FZtype
integer(kind=irg),INTENT(IN)            :: FZorder
logical                                 :: insideFZ

insideFZ = .FALSE.

! dealing with 180 rotations is needed only for 
! FZtypes 0 and 1; the other FZs are always finite.
  select case (FZtype)
    case (0)
      insideFZ = .TRUE.   ! all points are inside the FZ
    case (1)
      insideFZ = insideCyclicFZ(rod,FZtype,FZorder)        ! infinity is checked inside this function
    case (2)
      if (rod(4).ne.inftyd()) insideFZ = insideDihedralFZ(rod,FZorder)
    case (3)
      if (rod(4).ne.inftyd()) insideFZ = insideCubicFZ(rod,'tet')
    case (4)
      if (rod(4).ne.inftyd()) insideFZ = insideCubicFZ(rod,'oct')
    case (5) ! icosahedral symmetry
      if (rod(4).ne.inftyd()) insideFZ = insideIcosahedralFZ(rod)
    case (6) ! cubic-hexagonal misorientation FZ
      if (rod(4).ne.inftyd()) insideFZ = insideCubeHexFZ(rod)
    case (7)
!     if (rod(4).ne.inftyd) insideFZ = insideCubicFZ(rod,'oct')
    case (8)
!     if (rod(4).ne.inftyd) insideFZ = insideCubicFZ(rod,'oct')
  end select

end function IsinsideFZ


!--------------------------------------------------------------------------
!
! FUNCTION: insideIcosahedralFZ
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief does Rodrigues point lie inside cosahedral FZ?
!
!> @param rod Rodrigues coordinates  (double precision)
!
!> @date 03/24/17 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function insideIcosahedralFZ(rod) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: insideIcosahedralFZ

real(kind=dbl), INTENT(IN)        :: rod(4)
logical                           :: res

real(kind=dbl)                    :: dval, rv(3)
integer(kind=irg)                 :: i, j

res = .FALSE.

dval=0.32491969623290632616D0  ! sqrt(1-2/sqrt(5)))
rv(1:3) = rod(1:3)*rod(4)
j = 0
do i=1,12
  if (DOT_PRODUCT(IcoVertices(1:3,i),rv)+dval.ge.0.D0) j = j+1
end do
if (j.eq.12) res = .TRUE.

end function insideIcosahedralFZ

!--------------------------------------------------------------------------
!
! FUNCTION: insideCyclicFZ
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief does Rodrigues point lie inside cyclic FZ (for 2, 3, 4, and 6-fold)?
!
!> @param rod Rodrigues coordinates  (double precision)
!> @param FZtype symmetry type
!> @param FZorder depending on main symmetry axis
!
!> @date 05/12/14 MDG 1.0 original
!> @date 10/02/14 MDG 2.0 rewrite
!> @date 06/04/15 MDG 2.1 corrected infty to inftyd (double precision infinity)
!--------------------------------------------------------------------------
recursive function insideCyclicFZ(rod,FZtype,FZorder) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: insideCyclicFZ

use math

IMPLICIT NONE

real(kind=dbl), INTENT(IN)                :: rod(4)
integer(kind=irg), INTENT(IN)             :: FZtype
integer(kind=irg), INTENT(IN)             :: FZorder

logical                                   :: res

res = .FALSE.

if (rod(4).ne.inftyd()) then
  if ((FZtype.eq.1.).and.(FZorder.eq.2)) then
! check the y-component vs. tan(pi/2n)
    res = dabs(rod(2)*rod(4)) .le. BP(FZorder)
  else
! check the z-component vs. tan(pi/2n)
    res = dabs(rod(3)*rod(4)).le. BP(FZorder)
  end if
else
  if ((FZtype.eq.1.).and.(FZorder.eq.2)) then
    if(rod(2) .eq. 0.D0) res = .TRUE.
  else
    if (rod(3).eq.0.D0) res = .TRUE.
  end if
endif

end function insideCyclicFZ

!--------------------------------------------------------------------------
!
! FUNCTION: insideDihedralFZ
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief does Rodrigues point lie inside dihedral FZ (for 2, 3, 4, and 6-fold)?
!
!> @param rod Rodrigues coordinates (double precision)
!> @param order depending on main symmetry axis
!
!> @todo for now, we ignore here the fact that, among others, the 3m point group can be oriented in two ways;
!> @todo this should be fixable in the future with an additional optional argument
!
!> @date 05/12/14  MDG 1.0 original
!> @date 10/02/14  MDG 2.0 rewrite
!--------------------------------------------------------------------------
recursive function insideDihedralFZ(rod,order) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: insideDihedralFZ

use constants

IMPLICIT NONE

real(kind=dbl), INTENT(IN)                :: rod(4)
integer(kind=irg), INTENT(IN)             :: order

logical                                   :: res, c1, c2
real(kind=dbl)                            :: r(3)
real(kind=dbl),parameter                  :: r1 = 1.00D0
real(kind=dbl),allocatable                :: polygonvertex(:,:)
integer(kind=irg)                         :: inout
real(kind=dbl),parameter                  :: tol = 1.0D-8

if (rod(4).gt.1.5) then 
  res = .FALSE.
else
  r(1:3) = rod(1:3) * rod(4)

  ! first, check the z-component vs. tan(pi/2n)  (same as insideCyclicFZ)
  c1 = dabs(r(3)) .le. BP(order) + tol
  res = .FALSE.

  ! check the square boundary planes if c1=.TRUE.
  if (c1) then
    select case (order)
      case (2)
        c2 = (dabs(r(1)) .le. r1+tol).and.(dabs(r(2)) .le. r1 + tol)
      case (3)
        c2 =          dabs( LPssrt * r(1) + 0.5D0 * r(2)) .le. r1 + tol
        c2 = c2.and.( dabs( LPssrt * r(1) - 0.5D0 * r(2)) .le. r1 + tol)
        c2 = c2.and.( dabs(r(2)).le.r1 )
      case (4)
        c2 = (dabs(r(1)) .le. r1 + tol) .and. (dabs(r(2)) .le. r1 + tol)
        c2 = c2.and.((LPsr22 * dabs(r(1) + r(2)) .le. r1 + tol).and.(LPsr22 * dabs(r(1) - r(2)) .le. r1 + tol))
      case (6)
        c2 =          dabs( 0.5D0 * r(1) + LPssrt * r(2)) .le. r1 + tol
    
        c2 = c2.and.( dabs( LPssrt * r(1) + 0.5D0 * r(2)) .le. r1 + tol )

        c2 = c2.and.( dabs( LPssrt * r(1) - 0.5D0 * r(2)) .le. r1 + tol)

        c2 = c2.and.( dabs( 0.5D0 * r(1) - LPssrt * r(2)) .le. r1 + tol)
  
        c2 = c2.and.( dabs(r(2)) .le. r1 + tol )

        c2 = c2.and.( dabs(r(1)) .le. r1 + tol)

      ! add the 2-D quasi crystal type for 822, 1022, and 1222 rotational groups
      case (8)
        c2 = .FALSE.
        allocate(polygonvertex(order*2, 2))
        polygonvertex = 0.D0
        call getVertex(order, polygonvertex)
        inout = PNPOLY(r(1),r(2),polygonvertex(1:2*order,1),polygonvertex(1:2*order,2),order*2)
        if(inout .ge. 0) c2 = .TRUE.

      case (10)
        c2 = .FALSE.
        allocate(polygonvertex(order*2, 2))
        polygonvertex = 0.D0
        call getVertex(order, polygonvertex)
        inout = PNPOLY(r(1),r(2),polygonvertex(1:2*order,1),polygonvertex(1:2*order,2),order*2)
        if(inout .ge. 0) c2 = .TRUE.

      case(12)
        c2 = .FALSE.
        allocate(polygonvertex(order*2, 2))
        polygonvertex = 0.D0
        call getVertex(order, polygonvertex)
        inout = PNPOLY(r(1),r(2),polygonvertex(1:2*order,1),polygonvertex(1:2*order,2),order*2)
        if(inout .ge. 0) c2 = .TRUE.

    end select
    res = c2

  end if
end if

end function insideDihedralFZ

!--------------------------------------------------------------------------
!
! FUNCTION: insideCubicFZ
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief does Rodrigues point lie inside cubic FZ (octahedral or tetrahedral)?
!
!> @param rod Rodrigues coordinates  (double precision)
!> @param ot 'oct' or 'tet', depending on symmetry
!
!> @date 05/12/14 MDG 1.0 original
!> @date 10/02/14 MDG 2.0 rewrite
!> @date 01/03/15 MDG 2.1 correction of boundary error; simplification of octahedral planes
!> @date 06/04/15 MDG 2.2 simplified handling of components of r
!--------------------------------------------------------------------------
recursive function insideCubicFZ(rod,ot) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: insideCubicFZ

IMPLICIT NONE

real(kind=dbl), INTENT(IN)                :: rod(4)
character(3), INTENT(IN)                  :: ot

logical                                   :: res, c1, c2
real(kind=dbl)                            :: r(3)
real(kind=dbl),parameter                  :: r1  = 1.0D0
real(kind=dbl),parameter                  :: eps = 1.0D-4

r(1:3) = rod(1:3) * rod(4)

res = .FALSE.

! primary cube planes (only needed for octahedral case)
if (ot.eq.'oct') then
  c1 = (maxval(dabs(r)) - BP(4) .le. eps)
else 
  c1 = .TRUE.
end if

! octahedral truncation planes, both for tetrahedral and octahedral point groups
c2 = ((dabs(r(1))+dabs(r(2))+dabs(r(3))) - r1 .le. eps)

! if both c1 and c2, then the point is inside
if (c1.and.c2) res = .TRUE.

end function insideCubicFZ

!--------------------------------------------------------------------------
!
! FUNCTION: insideCubeHexFZ
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief does Rodrigues point lie inside combined cubic-hexagonal FZ?
!
!> @note For details on this test, see section 8 in "Representation of Orientation and
!> Disorientation data for Cubic, Hexagonal, Tetragonal, and Orthorhombic Crystals", A. Heinz
!> and P. Neumann, Acta Cryst. A47, 780-789 (1991)
!
!> @param rod Rodrigues coordinates  (double precision)
!
!> @date 04/01/17 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function insideCubeHexFZ(rod) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: insideCubeHexFZ

IMPLICIT NONE

real(kind=dbl), INTENT(IN)                :: rod(4)

logical                                   :: res
real(kind=dbl)                            :: r(3)
real(kind=dbl),parameter                  :: r1 = 0.414213562373095D0, r2 = 0.131652497587396D0, &
                                             alpha = 0.267949192431123D0, beta = 0.464101615137755D0
real(kind=dbl),parameter                  :: eps = 1.0D-3


r(1:3) = rod(1:3) * rod(4)

res = .FALSE.

if ( (r(2).ge.0.D0).and.(r(3).ge.0.D0) ) then 
  if ( ((alpha * (r(1)+r(3)) + r(2)) - beta .le. eps).and.( (alpha * (r(2)-r(3)) + r(1)) - beta .le. eps) ) then 
    if ( (r(1) - r1 .le. eps) .and. (r(2) - r1 .le. eps) .and. (r(3) - r2 .le. eps) ) res = .TRUE.
  end if
end if 

end function insideCubeHexFZ

!--------------------------------------------------------------------------
!
! FUNCTION: SampleRFZ
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Generate a uniform sampling of a Rodriguess FZ
!
!> @note This routine fills in a linked list FZlist of Rodrigues points that
!> are inside a specific fundamental zone determined by the sample point group;
!> this list can then be further dealt with in the calling program.
!>
!> Here's how you would use this routine in a main program:
!>
!> use so3
!>
!> integer(kind=irg)       :: FZcnt, nsteps, pgnum
!> type(FZpointd),pointer  :: FZlist, FZtmp
!>
!> nullify(FZlist)
!> FZcnt = 0
!> nsteps = 10
!> pgnum = 32
!> call sampleRFZ(nsteps, pgnum, FZcnt, FZlist)
!>
!> Then you can access all the entries in the list and, for instance, convert them to Euler angles...
!>
!> FZtmp => FZlist                        ! point to the top of the list
!> do i = 1, FZcnt                        ! loop over all entries
!>   eu = ro2eu(FZtmp%rod)                ! convert to Euler angles
!>   do something with eu                 ! for instance, write eu to a file
!>   FZtmp => FZtmp%next                  ! point to the next entry
!> end do
!>
!> If you just want to look at the first 10 entries on the list and show all other orientation representations:
!>
!> type(orientationtyped):: ot
!>
!> FZtmp => FZlist
!> do i = 1,10
!>   ot = init_orientation(FZtmp%rod,'ro')
!>   call print_orientation(ot)
!>   FZtmp => FZtmp%next
!> end do
!
!> @param nsteps number of steps along semi-edge in cubochoric grid
!> @param pgnum point group number to determine the appropriate Rodrigues fundamental zone
!> @param gridtype (input) 0 for origin-centered grid; 1 for grid with origin at box center
!> @param FZcnt (output) number of points inside fundamental zone
!> @param FZlist (output) linked list of points inside fundamental zone
!
!> @date 05/12/14 MDG 1.0 original
!> @date 10/02/14 MDG 2.0 rewrite, removed all globals, added function arguments
!> @date 09/15/15 MDG 2.1 removed explicit origin allocation; changed while to do loops.
!> @date 01/17/15 MDG 2.2 added gridtype option
!> @date 05/22/16 MDG 2.3 correction for monoclinic symmetry with twofold axis along b, not c !!!
!--------------------------------------------------------------------------
recursive subroutine SampleRFZ(nsteps, pgnum, gridtype, FZcnt, callable)
!DEC$ ATTRIBUTES DLLEXPORT :: SampleRFZ

use typedefs
use rotations

IMPLICIT NONE

integer(kind=irg), INTENT(IN)           :: nsteps
integer(kind=irg), INTENT(IN)           :: pgnum
integer(kind=irg), INTENT(IN)           :: gridtype
type(FZpointd),pointer                  :: FZlist               ! pointers
real(kind=dbl), allocatable             :: FZout(:,:)
integer(kind=irg),INTENT(OUT)           :: FZcnt                ! counts number of entries in linked list
external callable

real(kind=dbl)                          :: x, y, z, rod(4), delta, shift, sedge, ztmp
type(FZpointd), pointer                 :: FZtmp, FZtmp2
integer(kind=irg)                       :: FZtype, FZorder, i, j, k
logical                                 :: b

! cube semi-edge length s = 0.5D0 * LPsap
! step size for sampling of grid; total number of samples = (2*nsteps+1)**3
sedge = 0.5D0 * LPsap
delta = sedge / dble(nsteps)
if (gridtype.eq.0) then
  shift = 0.0D0
else
  shift = 0.5D0
end if

! set the counter to zero
FZcnt = 0

! make sure the linked lists are empty
if (associated(FZlist)) then
  FZtmp => FZlist%next
  FZtmp2 => FZlist
  do
    deallocate(FZtmp2)
    if (.not. associated(FZtmp) ) EXIT
    FZtmp2 => FZtmp
    FZtmp => FZtmp%next
  end do
  nullify(FZlist)
else
  nullify(FZlist)
end if

! determine which function we should call for this point group symmetry
FZtype = FZtarray(pgnum)
FZorder = FZoarray(pgnum)

! note that when FZtype is cyclic (1) and FZorder is 2, then we must rotate the
! rotation axis to lie along the b (y) direction, not z !!!!

! loop over the cube of volume pi^2; note that we do not want to include
! the opposite edges/facets of the cube, to avoid double counting rotations
! with a rotation angle of 180 degrees.  This only affects the cyclic groups.

 do i=-nsteps+1,nsteps
  x = (dble(i)+shift)*delta
  do j=-nsteps+1,nsteps
   y = (dble(j)+shift)*delta
   do k=-nsteps+1,nsteps
    z = (dble(k)+shift)*delta
! make sure that this point lies inside the cubochoric cell
    if (maxval( (/ abs(x), abs(y), abs(z) /) ).le.sedge) then

! convert to Rodrigues representation
      rod = cu2ro_d( (/ x, y, z /) )

! If insideFZ=.TRUE., then add this point to the linked list FZlist and keep
! track of how many points there are on this list
       b = IsinsideFZ(rod,FZtype,FZorder)
       if (b) then
        if (.not.associated(FZlist)) then
          allocate(FZlist)
          FZtmp => FZlist
! in DEBUG mode, there is a strange error here ...
! if gridtype = 0 then we must add the identity rotation here
! currently, the identity is not automatically recognized since the
! returns from several rotation routines also equal the identity when
! the point is actually invalid...
!          if (gridtype.eq.0) then
!            FZcnt = 1
!            FZtmp%rod = (/ 0.D0, 0.D0, 1.D0, 0.D0 /)
!            allocate(FZtmp%next)
!            FZtmp => FZtmp%next
!            nullify(FZtmp%next)
!          end if
        else
          allocate(FZtmp%next)
          FZtmp => FZtmp%next
        end if
        nullify(FZtmp%next)
! if monoclinic, then reorder the components !!!
!        if ((FZtype.eq.1).and.(FZorder.eq.2)) then
!          ztmp = rod(3)
!          rod(3) = rod(1)
!          rod(1) = rod(2)
!          rod(2) = ztmp
!        end if
        FZtmp%rod = rod
        FZtmp%gridpt(1:3) = (/i, j, k/)
        FZcnt = FZcnt + 1
       end if
    end if
  end do
 end do
end do

allocate(FZout(FZcnt,4))

FZtmp => FZlist

do i = 1,FZcnt
    FZout(i,1:4) = FZtmp%rod
    FZtmp => FZtmp%next
end do

call callable(FZout, FZcnt)

end subroutine SampleRFZ

!--------------------------------------------------------------------------
!
! SUBROUTINE: ReduceOrientationtoRFZ
!
!> @author Saransh Singh/Marc De Graef, Carnegie Mellon University
!
!> @brief Reduce an orientation (quaternions) to the Rodrigues Fundamental Zone
!
!> @param quin input quaternion in northern hemisphere
!> @param dict dict structure
!> @param FZtype Fundamental Zone type
!> @param FZorder Fundamental Zone order
!> @param quFZ output quaternion in fundamental zone
!> @param MFZ (optonal) apply MacKenzie cell
!
!> @date 07/29/16 MDG 1.0 original
!> @date 03/27/17 MDG 1.1 added checking of MacKenzie cell
!> @date 09/17/18 SS  2.0 modified for python interface for ODFpf package
!--------------------------------------------------------------------------
recursive subroutine ReduceOrientationtoRFZ(quin, Nqsym, Pm, FZtype, FZorder, quFZ)
!DEC$ ATTRIBUTES DLLEXPORT :: ReduceOrientationtoRFZ

use local
use rotations
use quaternions

IMPLICIT NONE

real(kind=dbl),INTENT(IN)               :: quin(4)
integer(kind=irg),INTENT(IN)            :: Nqsym
real(kind=dbl),INTENT(IN)               :: Pm(4,60)
integer(kind=irg),INTENT(IN)            :: FZtype
integer(kind=irg),INTENT(IN)            :: FZorder
real(kind=dbl),INTENT(OUT)              :: quFZ(4)

real(kind=dbl)                          :: Mu(4), qu(4), rod(4)
integer(kind=irg)                       :: i, j, Pmdims
real(kind=dbl)                          :: tol

tol = 5.0D-2

quFZ   = 0.D0
Pmdims = Nqsym
Mu = quin
Mu = Mu/NORM2(Mu)

if (Mu(1).lt.0.D0) Mu = -Mu
FZloop: do j=1,Pmdims
  qu = quat_mult_d(Pm(1:4,j),Mu)
  qu = qu/NORM2(qu)

  if (qu(1).lt.0.D0) qu = -qu
  rod = qu2ro_d(qu)

  !if(abs(rod(4)) .gt. 1.0D0+tol) rod(4) = inftyd()

  if (IsinsideFZ(rod,FZtype,FZorder))  EXIT FZloop

  ! we really should never get to the following line ...
  !if (j.eq.Pmdims) write (*,*) 'problem ... ',180.0*eu(1:3)/cPi,eu2ro(eu)
end do FZloop
quFZ = ro2qu_d(rod)

end subroutine ReduceOrientationtoRFZ

!--------------------------------------------------------------------------
!
! SUBROUTINE: getVertex
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief get vertices of RFZ for quasicrystals (dihedral symmetries)
!
!> @param order name of the Euler angle file (with usual path handling)
!> @param vertex the number of components in the returned linked list
!
!> @date 06/18/18 SS 1.0 original
!--------------------------------------------------------------------------
recursive subroutine getVertex(order, vertex)
!DEC$ ATTRIBUTES DLLEXPORT :: getVertex

use constants

IMPLICIT NONE

integer(kind=irg),INTENT(IN)            :: order
real(kind=dbl),INTENT(OUT)              :: vertex(2*order,2)

integer(kind=irg)                       :: ii
real(kind=dbl)                          :: th

do ii = 1,2*order
  th  = (dble(ii - 1)/dble(order) + 1.D0/2.D0/dble(order)) * cPi
  vertex(ii,1:2) = (/dcos(th), dsin(th)/)
end do

end subroutine getVertex
                                                                
!                                                                      
!     ..................................................................
!                                                                       
!        SUBROUTINE PNPOLY                                              
!                                                                       
!        PURPOSE                                                        
!           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            
!                                                                       
!        USAGE                                                          
!           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )                     
!                                                                       
!        DESCRIPTION OF THE PARAMETERS                                  
!           PX      - X-COORDINATE OF POINT IN QUESTION.                
!           PY      - Y-COORDINATE OF POINT IN QUESTION.                
!           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
!                     VERTICES OF POLYGON.                              
!           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
!                     VERTICES OF POLYGON.                              
!           N       - NUMBER OF VERTICES IN THE POLYGON.                
!           INOUT   - THE SIGNAL RETURNED:                              
!                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
!                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     
!                      1 IF THE POINT IS INSIDE OF THE POLYGON.         
!                                                                       
!        REMARKS                                                        
!           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      
!           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           
!           OPTIONALLY BE INCREASED BY 1.                               
!           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      
!           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    
!           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   
!           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              
!           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         
!           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM      
!           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   
!                                                                       
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
!           NONE                                                        
!                                                                       
!        METHOD                                                         
!           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  
!           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        
!           POINT IS INSIDE OF THE POLYGON.                             
!                                                                       
!     ..................................................................
!                                                                       
      RECURSIVE FUNCTION PNPOLY(PX,PY,XX,YY,N) RESULT(INOUT)

      IMPLICIT NONE

      REAL(KIND=DBL) PX, PY
      INTEGER(KIND=IRG) N

      REAL(KIND=DBL) X(200),Y(200),XX(N),YY(N)                                    
      LOGICAL MX,MY,NX,NY                                         
      INTEGER O,INOUT,I,J,MAXDIM                                                         
!      OUTPUT UNIT FOR PRINTED MESSAGES                                 
      DATA O/6/                                                         
      MAXDIM=200                                                        
      IF(N.LE.MAXDIM)GO TO 6                                            
      WRITE(O,7)                                                        
7     FORMAT('0WARNING:',I5,' TOO GREAT FOR THIS VERSION OF PNPOLY. RESULTS INVALID')                                                 
      RETURN                                                            
6     DO 1 I=1,N                                                        
      X(I)=XX(I)-PX                                                     
1     Y(I)=YY(I)-PY                                                     
      INOUT=-1                                                          
      DO 2 I=1,N                                                        
      J=1+MOD(I,N)                                                      
      MX=X(I).GE.0.0                                                    
      NX=X(J).GE.0.0                                                    
      MY=Y(I).GE.0.0                                                    
      NY=Y(J).GE.0.0                                                    
      IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2       
      IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3  
      INOUT=-INOUT                                                      
      GO TO 2                                                           
3     IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5                       
4     INOUT=0                                                           
      RETURN                                                            
5     INOUT=-INOUT                                                      
2     CONTINUE                                                          
      RETURN                                                            
      END FUNCTION PNPOLY

end module so3
