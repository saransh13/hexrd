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
! EMsoft:symmetry.f90
!--------------------------------------------------------------------------
!
! MODULE: symmetry
!
!> @author Saransh Singh, Lawrence Livermore National Lab
!
!> @brief everything symmetry related for sampling
!
!> @todo verify that this is correct for point groups with multiple settings, eg, 3m, 32, ...
!
!> @date 09/17/18 SS 1.0 original
!--------------------------------------------------------------------------
module symmetry

use local

! following are used to define the quaternion symmetry operators

real(kind=dbl),parameter        :: sq22=0.7071067811865475244D0 ! sqrt(2)/2
real(kind=dbl),parameter        :: sq32=0.8660254037844386467D0 ! sqrt(3)/2
real(kind=dbl),parameter        :: half=0.5D0                   ! 1/2

real(kind=dbl),parameter        :: SYM_Qsymop(4,35) = reshape( (/ &
                                1.D0, 0.D0, 0.D0, 0.D0, &       ! 1: identity operator
                                0.D0, 1.D0, 0.D0, 0.D0, &       ! 2: 180@[100]
                                0.D0, 0.D0, 1.D0, 0.D0, &       ! 3: 180@[010]
                                0.D0, 0.D0, 0.D0, 1.D0, &       ! 4: 180@[001]
                                sq22, sq22, 0.D0, 0.D0, &       ! 5: 90@[100]
                                sq22, 0.D0, sq22, 0.D0, &       ! 6: 90@[010]
                                sq22, 0.D0, 0.D0, sq22, &       ! 7: 90@[001]
                                sq22,-sq22, 0.D0, 0.D0, &       ! 8: 270@[100]
                                sq22, 0.D0,-sq22, 0.D0, &       ! 9: 270@[010]
                                sq22, 0.D0, 0.D0,-sq22, &       !10: 270@[001]
                                0.D0, sq22, sq22, 0.D0, &       !11: 180@[110]
                                0.D0,-sq22, sq22, 0.D0, &       !12: 180@[-110]
                                0.D0, 0.D0, sq22, sq22, &       !13: 180@[011]
                                0.D0, 0.D0,-sq22, sq22, &       !14: 180@[0-11]
                                0.D0, sq22, 0.D0, sq22, &       !15: 180@[101]
                                0.D0,-sq22, 0.D0, sq22, &       !16: 180@[-101]
                                half, half, half, half, &       !17: 120@[111]
                                half,-half,-half,-half, &       !18: 120@[-1-1-1]
                                half, half,-half, half, &       !19: 120@[1-11]
                                half,-half, half,-half, &       !20: 120@[-11-1]
                                half,-half, half, half, &       !21: 120@[-111]
                                half, half,-half,-half, &       !22: 120@[1-1-1]
                                half,-half,-half, half, &       !23: 120@[-1-11]
                                half, half, half,-half, &       !24: 120@[11-1]
                                sq32, 0.D0, 0.D0, half, &       !25:  60@[001]   (hexagonal/trigonal operators start here)
                                half, 0.D0, 0.D0, sq32, &       !26: 120@[001]
                                0.D0, 0.D0, 0.D0, 1.D0, &       !27: 180@[001]  (duplicate from above, but useful to keep it here)
                               -half, 0.D0, 0.D0, sq32, &       !28: 240@[001]
                               -sq32, 0.D0, 0.D0, half, &       !29: 300@[001]
                                0.D0, 1.D0, 0.D0, 0.D0, &       !30: 180@[100]
                                0.D0, sq32, half, 0.D0, &       !31: 180@[xxx]
                                0.D0, half, sq32, 0.D0, &       !32: 180@[xxx]
                                0.D0, 0.D0, 1.D0, 0.D0, &       !33: 180@[010]
                                0.D0,-half, sq32, 0.D0, &       !34: 180@[xxx]
                                0.D0,-sq32, half, 0.D0  &       !35: 180@[xxx]
                                /), (/4,35/) )

integer(kind=irg),parameter       :: PGrot(36) = (/1,1,3,3,3,6,6,6,9,9,9,12,12,12,12,16,16, &
                                                  18,18,18,21,21,21,24,24,24,24,28,28,30,30,30,33,34,35,36/)
contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: Symmetry_INIT
!
!> @author Saransh Singh/Marc De Graef, Lawrence Livermore National Lab/
!> Carnegie Mellon University
!
!> @brief initialize the symmetry operations in quaternion form for periodic as well as
!> aperiodic crystals (point groups 532, 822, 1022 and 1222)
!
!> @details For all details, see following paper:
!
!> @date 12/31/14 MDG 1.0 original
!> @date 01/06/15 MDG 1.1 added optional argument full
!> @date 02/06/15 MDG 1.2 removed full again after extensive testing; no need to use 2M operators
!> @date 01/04/18 MDG 1.3 added icosahedral symmetry operator to handle quasi-crystal computations.
!> @date 06/XX/18 SS  1.4 added octagonal, decagonal and dodecagonal symmetries
!> @date 09/17/18 SS  2.0 modified for python interfacing
!--------------------------------------------------------------------------
recursive subroutine Symmetry_INIT(Nqsym, Pm, pgnum)
!DEC$ ATTRIBUTES DLLEXPORT :: Symmetry_INIT

use local
use typedefs
use error

IMPLICIT NONE

real(kind=dbl),INTENT(OUT)              :: Pm(4,60)
integer(kind=irg),INTENT(OUT)           :: Nqsym
integer(kind=irg),INTENT(IN)            :: pgnum

integer(kind=irg)                       :: i, prot
real(kind=dbl)                          :: y1, y2


! here we need to analyze the rotational symmetry group, and copy the appropriate
! quaternion symmetry operators into the dict%Pm array

! first get the number of the rotational point group that corresponds to the crystal point group
prot = PGrot(pgnum)

! possible values for prot are: (/1,3,6,9,12,16,18,21,24,28,30/)
! corresponding to the point groups 1, 2, 222, 4, 422, 3, 32, 6, 622, 23, and 432 respectively

! identity operator is part of all point groups
Pm = 0.D0                  ! initialize all entries to zero
Pm(1:4,1) = SYM_Qsymop(1:4,1)

! select statement for each individual rotational point group (see typedefs.f90 for SYM_Qsymop definitions)
select case (prot)
        case(1)         ! 1 (no additional symmetry elements)
                Nqsym = 1
                Pm(1:4,2) = -Pm(1:4,1)

        case(3)         ! 2  (we'll assume that the two-fold axis lies along the e_y-axis)
                Nqsym = 2
                Pm(1:4,2) = SYM_Qsymop(1:4,3)

        case(6)         ! 222
                Nqsym = 4
                do i=2,4
                  Pm(1:4,i) = SYM_Qsymop(1:4,i)
                end do

        case(9)         ! 4
                Nqsym = 4
                Pm(1:4,2) = SYM_Qsymop(1:4,4)
                Pm(1:4,3) = SYM_Qsymop(1:4,7)
                Pm(1:4,4) = SYM_Qsymop(1:4,10)

        case(12)        ! 422
                Nqsym = 8
                Pm(1:4,2) = SYM_Qsymop(1:4,4)
                Pm(1:4,3) = SYM_Qsymop(1:4,7)
                Pm(1:4,4) = SYM_Qsymop(1:4,10)
                Pm(1:4,5) = SYM_Qsymop(1:4,2)
                Pm(1:4,6) = SYM_Qsymop(1:4,3)
                Pm(1:4,7) = SYM_Qsymop(1:4,11)
                Pm(1:4,8) = SYM_Qsymop(1:4,12)

        case(16)        ! 3
                Nqsym = 3
                Pm(1:4,2) = SYM_Qsymop(1:4,26)
                Pm(1:4,3) = SYM_Qsymop(1:4,28)


        case(18)        ! 32 (needs special handling)
                Nqsym = 6
                Pm(1:4,2) = SYM_Qsymop(1:4,26)
                Pm(1:4,3) = SYM_Qsymop(1:4,28)
                Pm(1:4,4) = SYM_Qsymop(1:4,30)
                Pm(1:4,5) = SYM_Qsymop(1:4,32)
                Pm(1:4,6) = SYM_Qsymop(1:4,34)

        case(21)        ! 6
                Nqsym = 6
                do i=25,29
                  Pm(1:4,i-23) = SYM_Qsymop(1:4,i)
                end do

        case(24)        ! 622
                Nqsym = 12
                do i=25,35
                  Pm(1:4,i-23) = SYM_Qsymop(1:4,i)
                end do

        case(28)        ! 23
                Nqsym = 12
                do i=2,4
                  Pm(1:4,i) = SYM_Qsymop(1:4,i)
                end do
                do i=17,24
                  Pm(1:4,4+(i-16)) = SYM_Qsymop(1:4,i)
                end do

        case(30)        ! 432
                Nqsym = 24
                do i=2,24
                  Pm(1:4,i) = SYM_Qsymop(1:4,i)
                end do
        case(33)        ! 532
                Nqsym = 60
                do i=2,60
                  Pm(1:4,i) = SYM_Qsymop(1:4,35+i)
                end do
        case(34)  ! 822
                Nqsym = 16
                do i = 1,15
                  Pm(1:4,i+1) = SYM_Qsymop(1:4,95+i)
                end do
        case(35)  ! 1022
                Nqsym = 20
                do i = 1,19
                  Pm(1:4,i+1) = SYM_Qsymop(1:4,110+i)
                end do
        case(36)  ! 1222
                Nqsym = 24
                do i = 1,23
                  Pm(1:4,i+1) = SYM_Qsymop(1:4,129+i)
                end do

        case default    ! this should never happen ...
                call FatalError('InitDictionaryIndexing','unknown rotational point group number')
end select

end subroutine Symmetry_INIT

end module symmetry