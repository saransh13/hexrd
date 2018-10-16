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
! EMsoft:Lambert.f90
!--------------------------------------------------------------------------
!
! MODULE: Lambert
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief everything that has to do with the modified Lambert projections
!
!> @details This module contains a number of projection functions for the modified
!> Lambert projection between square lattice and 2D hemisphere, hexagonal lattice
!> and 2D hemisphere, as well as the more complex mapping between a 3D cubic grid
!> and the unit quaternion hemisphere with positive scalar comoonent.  In addition, there
!> are some other projections, such as the stereographic one.  Each function is named
!> by the projection, the dimensionality of the starting grid, and the forward or inverse
!> character.  For each function, there is also a single precision and a double precision
!> version, but we use the interface formalism to have only a single call.  The Forward
!> mapping is taken to be the one from the simple grid to the curved grid.  Since the module
!> deals with various grids, we also add a few functions/subroutines that apply symmetry
!> operations on those grids.
!
!> @date 07/10/13   MDG 1.0 original
!> @date 07/12/13   MDG 1.1 added forward cube to ball to quaternion mappings
!> @date 08/01/13   MDG 1.2 added standard Lambert projection
!> @date 08/12/13   MDG 1.3 added inverse Lambert projections for Ball to Cube
!> @date 09/20/13   MDG 1.4 added ApplyLaueSymmetry
!> @date 08/29/15   MDG 1.5 small changes to hexagonal mapping routines; coordinate swap inside routines
!> @date 01/20/18   MDG 1.6 added Lambert interpolations routines
!--------------------------------------------------------------------------
module Lambert

use local

IMPLICIT NONE

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

contains

!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DSquareForwardSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief forward projection from 2D square to 3D hemisphere, single precision
!
!> @param xy 2D coordinates to be transformed (single precision)
!> @param ierr error status: 0=OK, 1=input point lies outside square grid bounds
!
!> @date 07/10/13   MDG 1.0 original
!> @date 08/31/15   MDG 1.1 coordinates are prescaled
!--------------------------------------------------------------------------
recursive function Lambert2DSquareForwardSingle(xy,ierr) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DSquareForwardSingle

IMPLICIT NONE

real(kind=sgl),INTENT(IN)               :: xy(2)
integer(kind=irg),INTENT(INOUT)         :: ierr
real(kind=sgl)                          :: res(3), q, qq, xy2(2)

xy2 = xy * sngl(LPssPio2)

ierr = 0
! check to make sure that the input point lies inside the square of edge length 2 sqrt(pi/2)
if (maxval(abs(xy2)).gt.LPssPio2) then
  res = (/ 0.0, 0.0, 0.0 /)
  ierr = 1
else
! Forward projection from square grid to Northern hemisphere.
! Equations (8) and (9) from D. Rosca, "New uniform grids on the sphere,"
! Astronomy & Astrophysics, 520, A63 (2010)

! deal with the origin first:
 if (maxval(abs(xy2)).eq.0.0) then
   res = (/ 0.0, 0.0, 1.0 /)
 else
  if (abs(xy2(1)).le.abs(xy2(2))) then
   q = 2.0*xy2(2)*LPsiPi*sqrt(LPsPi-xy2(2)*xy2(2))
   qq = xy2(1)*LPsPi*0.25/xy2(2)
   res = (/ q*sin(qq), q*cos(qq), 1.0-2.0*xy2(2)*xy2(2)*sngl(LPsiPi) /)
  else
   q = 2.0*xy2(1)*LPsiPi*sqrt(LPsPi-xy2(1)*xy2(1))
   qq = xy2(2)*LPsPi*0.25/xy2(1)
   res = (/ q*cos(qq), q*sin(qq), 1.0-2.0*xy2(1)*xy2(1)*sngl(LPsiPi) /)
  end if
  res = res/sqrt(sum(res*res))
 end if
end if

end function Lambert2DSquareForwardSingle


!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DSquareForwardDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief forward projection from 2D square to 3D hemisphere, double precision
!
!> @param xy 2D coordinates to be transformed (double precision)
!> @param ierr error status: 0=OK, 1=input point lies outside square grid bounds
!
!> @date 07/10/13   MDG 1.0 original
!> @date 08/31/15   MDG 1.1 coordinates are prescaled
!--------------------------------------------------------------------------
recursive function Lambert2DSquareForwardDouble(xy,ierr) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DSquareForwardDouble

IMPLICIT NONE

real(kind=dbl),INTENT(IN)               :: xy(2)
integer(kind=irg),INTENT(INOUT)         :: ierr
real(kind=dbl)                          :: res(3), q, qq, xy2(2)

xy2 = xy * LPssPio2

ierr = 0
! check to make sure that the input point lies inside the square of edge length 2 sqrt(pi)
if (maxval(dabs(xy2)).gt.LPssPio2) then
  res = (/ 0.D0, 0.D0, 0.D0 /)
  ierr = 1   ! input point does not lie inside square with edge length 2 sqrt(pi/2)
else
! Forward projection from square grid to Northern hemisphere.
! Equations (8) and (9) from D. Rosca, "New uniform grids on the sphere,"
! Astronomy & Astrophysics, 520, A63 (2010)

! deal with the origin first:
 if (maxval(abs(xy2)).eq.0.0) then
   res = (/ 0.D0, 0.D0, 1.D0 /)
 else
  if (dabs(xy2(1)).le.dabs(xy2(2))) then
   q = 2.D0*xy2(2)*LPsiPi*dsqrt(LPsPi-xy2(2)*xy2(2))
   qq = xy2(1)*LPsPi*0.25D0/xy2(2)
   res = (/ q*dsin(qq), q*dcos(qq), 1.D0-2.D0*xy2(2)*xy2(2)*LPsiPi /)
  else
   q = 2.D0*xy2(1)*LPsiPi*dsqrt(LPsPi-xy2(1)*xy2(1))
   qq = xy2(2)*LPsPi*0.25D0/xy2(1)
   res = (/ q*dcos(qq), q*dsin(qq), 1.D0-2.D0*xy2(1)*xy2(1)*LPsiPi /)
  end if
  res = res/dsqrt(sum(res*res))
 end if
end if

end function Lambert2DSquareForwardDouble


!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DSquareInverseSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief inverse projection from 3D hemisphere to 2D square, single precision
!
!> @note IMPORTANT: the calling routine must keep track of the sign of xyz(3);
!> this routine only uses the absolute value !
!
!> @param xyz 3D coordinates to be transformed (single precision)
!> @param ierr error status: 0=OK, 1=input point has norm different from 1
!
!> @date 07/10/13   MDG 1.0 original
!> @date 08/31/15   MDG 1.1 return scaled coordinates
!--------------------------------------------------------------------------
recursive function Lambert2DSquareInverseSingle(xyz,ierr) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DSquareInverseSingle

IMPLICIT NONE

real(kind=sgl),INTENT(IN)               :: xyz(3)
integer(kind=irg),INTENT(INOUT)         :: ierr
real(kind=sgl)                          :: res(2), q
real(kind=sgl),parameter                :: eps = 1.0E-6

ierr = 0
! check to make sure that the input point lies on the unit sphere
if (abs(1.0-sum(xyz**2)).gt.eps) then
  res = (/ 0.0, 0.0 /)
  ierr = 1
else
! intercept the points (0,0,+-1)
  if (abs(xyz(3)).eq.1.0) then
    res = (/ 0.0, 0.0 /)
  else
    if (abs(xyz(2)).le.abs(xyz(1))) then
      q = abs(xyz(1))/xyz(1) * sqrt(2.0*(1.0-abs(xyz(3))))
      res = (/ q * sngl(LPssPi2), q * atan(xyz(2)/xyz(1))/sngl(LPssPi2) /)
    else
      q = abs(xyz(2))/xyz(2) * sqrt(2.0*(1.0-abs(xyz(3))))
      res = (/  q * atan(xyz(1)/xyz(2))/sngl(LPssPi2), q * sngl(LPssPi2) /)
    end if
  end if
end if

res = res / sngl(LPssPio2)

end function Lambert2DSquareInverseSingle


!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DSquareInverseDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief inverse projection from 3D hemisphere to 2D square, double precision
!
!> @note IMPORTANT: the calling routine must keep track of the sign of xyz(3);
!> this routine only uses the absolute value !
!
!> @param xyz 3D coordinates to be transformed (double precision)
!> @param ierr error status: 0=OK, 1=input point has norm different from 1
!
!> @date 07/10/13   MDG 1.0 original
!> @date 08/31/15   MDG 1.1 return scaled coordinates
!--------------------------------------------------------------------------
recursive function Lambert2DSquareInverseDouble(xyz,ierr) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DSquareInverseDouble

IMPLICIT NONE

real(kind=dbl),INTENT(IN)               :: xyz(3)
integer(kind=irg),INTENT(INOUT) :: ierr
real(kind=dbl)                          :: res(2), q
real(kind=dbl),parameter                :: eps = 1.0D-12

ierr = 0
! check to make sure that the input point lies on the unit sphere
if (dabs(1.D0-sum(xyz**2)).gt.eps) then
  res = (/ 0.D0, 0.D0 /)
  ierr = 1
else
! intercept the points (0,0,+-1)
  if (dabs(xyz(3)).eq.1.D0) then
    res = (/ 0.D0, 0.D0 /)
  else
    if (dabs(xyz(2)).le.dabs(xyz(1))) then
      q = dabs(xyz(1))/xyz(1) * dsqrt(2.D0*(1.D0-dabs(xyz(3))))
      res = (/ q * LPssPi2, q * datan(xyz(2)/xyz(1))/LPssPi2 /)
    else
      q = dabs(xyz(2))/xyz(2) * dsqrt(2.D0*(1.D0-dabs(xyz(3))))
      res = (/  q * datan(xyz(1)/xyz(2))/LPssPi2, q * LPssPi2 /)
    end if
  end if
end if

res = res / LPssPio2

end function Lambert2DSquareInverseDouble

!--------------------------------------------------------------------------
! the functions below deal with the hexagonal to 2D hemisphere projection
!
! all derivations and equations can be found in
!
! D. Rosca and M. De Graef, "Area-preserving projections from hexagonal and triangular
! domains to the sphere and applications to electron back-scatter diffraction pattern simulations,"
! Modelling Simul. Mater. Sci. Eng. 21 (2013) 055021.
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DHexForwardSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief forward projection from 2D hexagon to 3D hemisphere, single precision
!
!> @param xy 2D coordinates to be transformed (single precision)
!> @param ierr error status: 0=OK, 1=input point outside hexagon
!
!> @date 07/10/13   MDG 1.0 original
!> @date 08/29/15   MDG 1.1 debug
!> @date 08/30/15   MDG 1.2 moved grid-to-cartesian coordinate transformation inside routine
!--------------------------------------------------------------------------
recursive function Lambert2DHexForwardSingle(xy,ierr) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DHexForwardSingle

IMPLICIT NONE

real(kind=sgl),INTENT(IN)       :: xy(2)
integer(kind=irg),INTENT(INOUT):: ierr

real(kind=sgl)                  :: res(3), q, XX, YY, xp, yp, XY2(2), xyc(2)
integer(kind=irg)               :: ks

 ierr = 0

 xyc = (/ xy(1) - 0.5 * xy(2), xy(2) * sngl(LPssrt) /) * sngl(LPspreg)

 if (maxval(abs(xyc)).eq.0.0) then
  res = (/ 0.0, 0.0, 1.0 /)
 else
! flip coordinates
  XY2 = (/ xyc(2), xyc(1) /)
! determine in which sextant this point lies
  ks = GetSextantSingle(XY2)

  select case (ks)
    case (0,3)
        q = XY2(2)*LPsprec/XY2(1)
        XX = LPspreb*XY2(1)*cos(q)
        YY = LPspreb*XY2(1)*sin(q)
    case (1,4)
        xp = XY2(1)+LPsrtt*XY2(2)
        yp = XY2(1)*LPspred/xp
        XX = LPsprea*xp*sin(yp)
        YY = LPsprea*xp*cos(yp)
    case (2,5)
        xp = XY2(1)-LPsrtt*XY2(2)
        yp = XY2(1)*LPspred/xp
        XX = LPsprea*xp*sin(yp)
        YY = -LPsprea*xp*cos(yp)
  end select
  q = XX**2+YY**2
! does the point lie outside the hexagon ?
  if (q.gt.4.0) then
    res = (/ 0.0, 0.0, 0.0 /)
    ierr = 1
  else
    res = (/ 0.50*XX*sqrt(4.0-q), 0.50*YY*sqrt(4.0-q), 1.0-0.5*q /)
  end if

! and flip the x and y coordinates
  xp = res(1)
  res(1) = res(2)
  res(2) = xp
 end if

end function Lambert2DHexForwardSingle

!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DHexForwardDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief forward projection from 2D hexagon to 3D hemisphere, double precision
!
!> @param xy 2D coordinates to be transformed (double precision)
!> @param ierr error status: 0=OK, 1=input point outside hexagon
!
!> @date 07/10/13   MDG 1.0 original
!> @date 08/29/15   MDG 1.1 debug
!> @date 08/30/15   MDG 1.2 moved grid-to-cartesian coordinate transformation inside routine
!--------------------------------------------------------------------------
recursive function Lambert2DHexForwardDouble(xy,ierr) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DHexForwardDouble

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: xy(2)
integer(kind=irg),INTENT(INOUT) :: ierr

real(kind=dbl)                  :: res(3), q, XX, YY, xp, yp, XY2(2), xyc(2)
integer(kind=irg)               :: ks

  ierr = 0

 xyc = (/ xy(1) - 0.5D0 * xy(2), xy(2) * LPssrt /) * LPspreg

 if (maxval(dabs(xyc)).eq.0.D0) then
  res = (/ 0.D0, 0.D0, 1.D0 /)
 else
! flip coordinates
  XY2 = (/ xyc(2), xyc(1) /)
! determine in which sextant this point lies
  ks = GetSextantDouble(XY2)

  select case (ks)
    case (0,3)
        q = XY2(2)*LPsprec/XY2(1)
        XX = LPspreb*XY2(1)*dcos(q)
        YY = LPspreb*XY2(1)*dsin(q)
    case (1,4)
        xp = XY2(1)+LPsrtt*XY2(2)
        yp = XY2(1)*LPspred/xp
        XX = LPsprea*xp*dsin(yp)
        YY = LPsprea*xp*dcos(yp)
    case (2,5)
        xp = XY2(1)-LPsrtt*XY2(2)
        yp = XY2(1)*LPspred/xp
        XX = LPsprea*xp*dsin(yp)
        YY = -LPsprea*xp*dcos(yp)
  end select
  q = XX**2+YY**2

! does the point lie outside the hexagon ?
  if (q.gt.4.D0) then
    res = (/ 0.D0, 0.D0, 0.D0 /)
    ierr = 1
  else
    res = (/ 0.5D0*XX*dsqrt(4.D0-q), 0.5D0*YY*dsqrt(4.D0-q), 1.D0-0.5D0*q /)
  end if

! and flip the x and y coordinates
  xp = res(1)
  res(1) = res(2)
  res(2) = xp
 end if

end function Lambert2DHexForwardDouble


!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DHexInverseSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief inverse projection from 3D hemisphere to 2D hexagon, single precision
!
!> @note Important: The calling program must keep track of the sign of xyz(3), since this
!> routine will take the absolute value |xyz(3)| for the z-component of the input vector!
!
!> @param xyz 3D coordinates to be transformed (single precision)
!> @param ierr error status: 0=OK, 1=input point not normalized
!
!> @date 07/10/13   MDG 1.0 original
!> @date 05/04/15   MDG 1.1 correction to abs(X)/X for sign(X)
!> @date 05/09/15   MDG 1.2 added slight offset to XX,YY to avoid interpolation issues
!> @date 08/30/15   MDG 1.3 moved grid-to-cartesian coordinate transformation inside routine
!--------------------------------------------------------------------------
recursive function Lambert2DHexInverseSingle(xyz,ierr) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DHexInverseSingle

IMPLICIT NONE

real(kind=sgl),INTENT(IN)       :: xyz(3)
integer(kind=irg),INTENT(INOUT) :: ierr

real(kind=sgl)                  :: res(2), q, qq, XX, YY, xxx, yyy, sgnX, XYZ2(3), xy(2)
integer(kind=irg)               :: ks
real(kind=sgl),parameter        :: eps = 1.0E-7, eps2 = 1.0E-4

ierr = 0
! check to make sure that the input point lies on the unit sphere
if (abs(1.0-sum(xyz**2)).gt.eps) then
  res = (/ 0.0, 0.0 /)
  ierr = 1
else
 if (abs(xyz(3)).eq.1.0) then
  res = (/ 0.0, 0.0 /)
 else
! flip x and y components, and take the | | of the third component.
  XYZ2 = (/ xyz(2), xyz(1), abs(xyz(3)) /)

! first do the Lambert projection
  q = sqrt(2.0/(1.0+XYZ2(3)))
  XX = q * XYZ2(1)+eps2
  YY = q * XYZ2(2)+eps2

! determine in which sextant this point lies
  ks = GetSextantSingle( (/ XX, YY /) )
  sgnX = sqrt(XX**2+YY**2)
  if (XX.lt.0.0) sgnX=-sgnX

! then perform the inverse to the hexagonal grid
  select case (ks)
    case (0,3)
        q = LPspree * sgnX
        xxx = q * LPssPio2
        if (XX.eq.0.0) then
          yyy = LPspref * LPsPi * 0.5
        else
          yyy = q * LPspref * atan(YY/XX)
        end if
    case (1,4)
        q = LPsprea * sgnX
        qq= atan((YY-LPsrtt*XX)/(XX+LPsrtt*YY))
        xxx = q * LPsrtt *( LPsPi/6.0 - qq )
        yyy = q * ( 0.5*LPsPi + qq )
    case (2,5)
        q = LPsprea * sgnX
        qq= atan((YY+LPsrtt*XX)/(XX-LPsrtt*YY))
        xxx = q * LPsrtt *( LPsPi/6.0 + qq )
        yyy = q * ( -0.5*LPsPi + qq )
  end select
! and flip the coordinates
  res = (/ yyy, xxx /)
 end if
end if

! and finally, transform the coordinates back to the hexagonal grid
xy = res
res = (/ xy(1) + xy(2) * sngl(LPsisrt), xy(2) * 2.0 * sngl(LPsisrt) /) / LPspreg

end function Lambert2DHexInverseSingle


!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DHexInverseDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief inverse projection from 3D hemisphere to 2D hexagon, double precision
!
!> @note Important: The calling program must keep track of the sign of xyz(3), since this
!> routine will take the absolute value |xyz(3)| for the z-component of the input vector!
!
!> @param xyz 3D coordinates to be transformed (double precision)
!> @param ierr error status: 0=OK, 1=input point not normalized
!
!> @date 07/10/13   MDG 1.0 original
!> @date 05/04/15   MDG 1.1 correction to abs(X)/X for sign(X)
!> @date 05/09/15   MDG 1.2 added slight offset to XX,YY to avoid interpolation issues
!> @date 08/30/15   MDG 1.3 moved grid-to-cartesian coordinate transformation inside routine
!--------------------------------------------------------------------------
recursive function Lambert2DHexInverseDouble(xyz,ierr) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: Lambert2DHexInverseDouble

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: xyz(3)
integer(kind=irg),INTENT(INOUT) :: ierr

real(kind=dbl)                  :: res(2), q, qq, XX, YY, xxx, yyy, sgnX, XYZ2(3), xy(2)
integer(kind=irg)               :: ks
real(kind=dbl),parameter        :: eps = 1.0E-12, eps2 = 1.0E-4

ierr = 0
! check to make sure that the input point lies on the unit sphere
if (dabs(1.0-sum(xyz**2)).gt.eps) then
  res = (/ 0.D0, 0.D0 /)
  ierr = 1
else
 if (dabs(xyz(3)).eq.1.D0) then
  res = (/ 0.D0, 0.D0 /)
 else
! flip x and y components, and take the | | of the third component.
  XYZ2 = (/ xyz(2), xyz(1), dabs(xyz(3)) /)

! first do the Lambert projection
  q = dsqrt(2.D0/(1.D0+XYZ2(3)))
  XX = q * XYZ2(1)+eps2
  YY = q * XYZ2(2)+eps2

! determine in which sextant this point lies
  ks = GetSextantDouble( (/ XX, YY /) )
  sgnX = dsqrt(XX**2+YY**2)
  if (XX.lt.0.D0) sgnX=-sgnX

! then perform the inverse to the hexagonal grid
  select case (ks)
    case (0,3)
        q = LPspree * sgnX
        xxx = q * LPssPio2
        if (XX.eq.0.0) then
          yyy = LPspref * LPsPi * 0.5D0
        else
          yyy = q * LPspref * datan(YY/XX)
        end if
    case (1,4)
        q = LPsprea * sgnX
        qq= datan((YY-LPsrtt*XX)/(XX+LPsrtt*YY))
        xxx = q * LPsrtt *( LPsPi/6.D0 - qq )
        yyy = q * ( 0.5D0*LPsPi + qq )
    case (2,5)
        q = LPsprea * sgnX
        qq= datan((YY+LPsrtt*XX)/(XX-LPsrtt*YY))
        xxx = q * LPsrtt *( LPsPi/6.D0 + qq )
        yyy = q * ( -0.5D0*LPsPi + qq )
  end select
! and flip the coordinates back
    res = (/ yyy, xxx /)
end if
end if

! and finally, transform the coordinates back to the hexagonal grid
xy = res
res = (/ xy(1) + xy(2) * LPsisrt, xy(2) * 2.D0 * LPsisrt /) / LPspreg

end function Lambert2DHexInverseDouble

!--------------------------------------------------------------------------
!
! FUNCTION: GetSextantSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief determine to which sextant a point in hexagonal coordinates belongs
!
!> @param xy 2D coordinates to be considered (single precision)
!
!> @date 11/21/12    MDG 1.0 original
!> @date 08/29/15    MDG 1.1 debug
!--------------------------------------------------------------------------
recursive function GetSextantSingle(xy) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: GetSextantSingle

IMPLICIT NONE

real(kind=sgl),INTENT(IN)       :: xy(2)
integer(kind=irg)               :: res

real(kind=sgl)                  :: xx

xx = abs(xy(1)*LPsisrt)        ! |x| / sqrt(3)

if (xy(1).ge.0.0) then
  if (abs(xy(2)).le.xx) then
    res = 0
  else
    if (xy(2).gt.xx) then
      res = 1
    else
      res = 5
    end if
  end if
else
  if (abs(xy(2)).le.xx) then
    res = 3
  else
    if (xy(2).gt.xx) then
      res = 2
    else
      res = 4
    end if
  end if
end if

end function GetSextantSingle

!--------------------------------------------------------------------------
!
! FUNCTION: GetSextantDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief determine to which sextant a point in hexagonal coordinates belongs
!
!> @param xy 2D coordinates to be considered (double precision)
!
!> @date 11/21/12    MDG 1.0 original
!> @date 08/29/15    MDG 1.1 debug
!--------------------------------------------------------------------------
recursive function GetSextantDouble(xy) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: GetSextantDouble

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: xy(2)
integer(kind=irg)               :: res

real(kind=dbl)                  :: xx

xx = dabs(xy(1)*LPsisrt)       ! |x| / sqrt(3)

if (xy(1).ge.0.D0) then
  if (dabs(xy(2)).le.xx) then
    res = 0
  else
    if (xy(2).gt.xx) then
      res = 1
    else
      res = 5
    end if
  end if
else
  if (dabs(xy(2)).le.xx) then
    res = 3
  else
    if (xy(2).gt.xx) then
      res = 2
    else
      res = 4
    end if
  end if
end if

end function GetSextantDouble

end module Lambert