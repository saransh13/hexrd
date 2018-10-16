! ###################################################################
! Copyright (c) 2014-2018, Marc De Graef/Carnegie Mellon University
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
! EMsoft:local.f90
!--------------------------------------------------------------------------
!
! MODULE: local
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief definitions of single and double precision, general constants and variables
!
!> @details
!> defines the kind-parameters for short and long integers, and single/double
!> precision reals, as well as a few other global variables (path variables etc).
!
!> @date 1/8/99   MDG 1.0 original
!> @date 5/6/01   MDG 2.0 f90
!> @date 11/27/01 MDG 2.1 added sgl and dbl kinds
!> @date 12/08/01 MDG 2.2 added CTEMsoft subroutine
!> @date 03/19/13 MDG 3.0 rewrite of entire package
!> @date 05/16/13 MDG 3.1 added stdout
!> @date 01/10/14 MDG 4.0 new version
!> @date 06/05/14 MDG 4.1 added comments about global variables in Release 3.0
!> @date 03/29/14 MDG 4.2 added path variable for .xtal files (which are now HDF5 files)
!> @date 04/05/15 MDG 4.3 added additional path variables and an EMsoft_path_init routine
!> @date 11/20/15 MDG 4.4 added UserXXX variables to EMSoftConfig.json file
!> @date 02/16/16 MDG 4.5 added path delimiter conversion routines for Windows implementation
!> @date 02/25/16 MDG 4.6 added automatic generation of .config/EMsoft/EMsoftConfig.json file
!> @date 07/02/16 MDG 5.0 replaced all remaining globals by function calls to please the intel compiler
!> @date 10/25/16 MDG 5.1 added functionality for initializing EMsoftConfig.json on Windows
!> @date 08/16/17 MDG 5.2 replaced all string constants by string constants declared in stringconstants.f90
!> @date 10/28/17 MDG 6.0 changed handling of config parameters to a single structure read in at program start
!--------------------------------------------------------------------------

!---------------------------
! THIS FILE IS AUTOMATICALLY GENERATED DURING CMAKE TIME. THE ORIGINAL FILE
! LOCATED AT @CTEMSoftLib_SOURCE_DIR@/local.f90.in
! YOU NEED TO MAKE CHANGES TO THAT FILE. ANY CHANGES MADE TO THIS FILE WILL
! BE OVER WRITTEN THE NEXT TIME CMAKE IS EXECUTED.
!---------------------------

module local

!> @note This module must be "use"d by every program, subroutine, and function.
!> These are the only global variables used by the EMsoft package.

! The entire EMsoft package should be processor independent.  This can
! be accomplished by the use of the "kind" parameters.

! Define the "kind" parameters for single and double precision reals,
!> single precision real kind parameter
  integer,parameter                     :: sgl = SELECTED_REAL_KIND(p=6,r=37)
!> double precision real kind parameter
  integer,parameter                     :: dbl = SELECTED_REAL_KIND(p=13,r=200)
!DEC$ ATTRIBUTES DLLEXPORT :: sgl
!DEC$ ATTRIBUTES DLLEXPORT :: dbl

! Define the "kind" parameters for short and regular integers,
!> short integer kind parameter
  integer,parameter                     :: ish = SELECTED_INT_KIND(3)
!> long integer kind parameter
  integer,parameter                     :: irg = SELECTED_INT_KIND(9)
!> long long kind parameter
  integer,parameter                     :: ill = SELECTED_INT_KIND(12)
!DEC$ ATTRIBUTES DLLEXPORT :: ish
!DEC$ ATTRIBUTES DLLEXPORT :: irg
!DEC$ ATTRIBUTES DLLEXPORT :: ill

!> character type used for json routines
  integer,parameter :: jsonCK = selected_char_kind('DEFAULT')
!DEC$ ATTRIBUTES DLLEXPORT :: jsonCK

!> standard string length for filenames
  integer(kind=irg),parameter           :: fnlen = 512
!DEC$ ATTRIBUTES DLLEXPORT :: fnlen

!> reserved IO unit identifiers for postscript (20) and data (21-23)
  integer(kind=irg), parameter          :: psunit = 20, dataunit = 21, dataunit2 = 22, dataunit3 = 23
!DEC$ ATTRIBUTES DLLEXPORT :: psunit
!DEC$ ATTRIBUTES DLLEXPORT :: dataunit
!DEC$ ATTRIBUTES DLLEXPORT :: dataunit2
!DEC$ ATTRIBUTES DLLEXPORT :: dataunit3

end module local
