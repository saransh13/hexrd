! ###################################################################
! Copyright (c) 2013-2017, Marc De Graef/Carnegie Mellon University
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

module typedefs

use local 
use constants

IMPLICIT NONE

! the "orientation" type contains entries for all rotation and orientation representations
type orientationtype
  real(kind=sgl)        :: eulang(3)            ! Bunge Euler angles in radians
  real(kind=sgl)        :: om(3,3)              ! 3x3 matrix
  real(kind=sgl)        :: axang(4)             ! axis-angle pair (angle in rad, component 4; axis in direction cosines)
  real(kind=sgl)        :: rodrigues(4)         ! Rodrigues vector (stored as direction cosines and length, to allow for Infinity)
  real(kind=sgl)        :: quat(4)              ! quaternion representation (q(1) is scalar part, q(2:4) vector part)
  real(kind=sgl)        :: homochoric(3)        ! homochoric representation according to Frank's paper
  real(kind=sgl)        :: cubochoric(3)        ! cubic grid representation (derived from homochoric)
  real(kind=sgl)        :: stereographic(3)     ! 3D stereographic  [added 10/05/17]
end type orientationtype


! double precision version
type orientationtyped
  real(kind=dbl)        :: eulang(3)            ! Bunge Euler angles in radians
  real(kind=dbl)        :: om(3,3)              ! 3x3 matrix
  real(kind=dbl)        :: axang(4)             ! axis-angle pair (angle in rad, component 4; axis in direction cosines)
  real(kind=dbl)        :: rodrigues(4)         ! Rodrigues vector (stored as direction cosines and length, to allow for Infinity)
  real(kind=dbl)        :: quat(4)              ! quaternion representation (q(1) is scalar part, q(2:4) vector part)
  real(kind=dbl)        :: homochoric(3)        ! homochoric representation according to Frank's paper
  real(kind=dbl)        :: cubochoric(3)        ! cubic grid representation (derived from homochoric)
  real(kind=dbl)        :: stereographic(3)     ! 3D stereographic  [added 10/05/17]
end type orientationtyped

! type definition for linked list of Rodrigues-Frank vectors (used in so3.module)
type FZpointd
        real(kind=dbl)          :: rod(4)        ! Rodrigues-Frank vector [nx, ny, nz, tan(omega/2) ]
        real(kind=dbl)          :: trod(4)       ! second Rodrigues-Frank vector; can be used for coordinate transformations
        integer(kind=irg)       :: gridpt(3)     ! coordinates of grid point ! added on 06/19/18 by SS
        type(FZpointd),pointer  :: next          ! link to next point
end type FZpointd

! Dimension of the state 
integer, parameter :: ns = 4 
!DEC$ ATTRIBUTES DLLEXPORT :: ns

! Default seed vector 
integer, parameter, dimension(ns) :: default_seed = (/ 521288629, 362436069, 16163801, 1131199299 /) 
!DEC$ ATTRIBUTES DLLEXPORT :: default_seed

! A data type for storing the state of the RNG 
type :: rng_t 
  integer, dimension(ns) :: state = default_seed 
end type rng_t

end module typedefs
