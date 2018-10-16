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
! EMsoft:error.f90
!--------------------------------------------------------------------------
!
! MODULE: error
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief error handling routines
!
!> @details  Just a simple routine that reports errors and aborts the program
!
!> @date 01/05/99 MDG 1.0 original
!> @date 05/19/01 MDG 2.0 f90 version
!> @date 11/27/01 MDG 2.1 added kind support
!> @date 03/19/13 MDG 3.0 minor changes
!> @date 06/05/14 MDG 4.0 minor changes
!--------------------------------------------------------------------------

module error

use local
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit

public

interface WriteValue
        module procedure WriteValueIntShort
        module procedure WriteValueIntLong
        module procedure WriteValueIntLongLong
        module procedure WriteValueRealSingle
        module procedure WriteValueRealDouble
        module procedure WriteValueRealComplex
        module procedure WriteValueString
end interface WriteValue

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: FatalError
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Write error message and abort program
!
!> @param var1 routine name string
!> @param var2 explanation string
!> @param stdout optional output unit identifier
!
!> @date   10/13/98 MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/19/13 MDG 3.0 updated for new io routines
!> @date   06/05/14 MDG 4.0 added stdout argument
!> @date   03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
recursive subroutine FatalError(var1,var2)
!DEC$ ATTRIBUTES DLLEXPORT :: FatalError


IMPLICIT NONE

character(*), INTENT(IN)                :: var1  !< first part of error message (routine name)
character(*), INTENT(IN)                :: var2  !< second part of error message (brief explanation)

 call Message(' ----> Fatal error in routine '//var1//': '//var2, frm='(//A//)') 
 stop '  Progam ended abnormally'

end subroutine

!--------------------------------------------------------------------------
!
! SUBROUTINE: Message
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief dump a message to standard output
!
!> @details Simple routine to print a string on the standard output, with optional formatting
!> instructions, for instance if one wants an empty line before (frm='(/A)') or after (frm='(A/)') 
!> the string.  Note that one can include the name of the optional variable in the subroutine
!> call, as in:
!> call Message('this is a string', frm='(//A//)' , stdout = 22)
!> this makes it clear that frm and stdout are optional variables.
! 
!> @param mess message string
!> @param frm optional string formatting command
!> @param stdout optional output unit identifier
!
!> @date 01/05/99 MDG 1.0 original
!> @date 05/19/01 MDG 2.0 f90 version
!> @date 03/19/13 MDG 3.0 made argument optional and introduced default format
!> @date 06/05/14 MDG 4.0 added stdout and mess as mandatory arguments
!> @date 03/29/18 MDG 4.1 removed stdout argument
!--------------------------------------------------------------------------
recursive subroutine message(mess,frm)
!DEC$ ATTRIBUTES DLLEXPORT :: Message

character(*),INTENT(IN)                 :: mess         !< message string
character(*),OPTIONAL,INTENT(IN)        :: frm          !< optional formatting string

! default format or not ?
if (PRESENT(frm)) then
 write (6,fmt=frm) trim(mess)
else    ! default output format: a simple string
 write (6,fmt="(A)") trim(mess)
end if 

end subroutine Message

! ###################################################################
! writing routines
! ###################################################################

! ###################################################################
! 
!  subroutine WriteValueString 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a string
!
!> @param Qstring question string
!> @param out_string output string
!> @param frm optional formatting argument
!
!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
recursive subroutine WriteValueString(Qstring, out_string, frm)
!DEC$ ATTRIBUTES DLLEXPORT :: WriteValueString

character(*),INTENT(IN)                         :: Qstring 
character(*),INTENT(IN)                         :: out_string
character(*),INTENT(IN),OPTIONAL                :: frm

! send Qstring to the output only if it is non-zero length
if (len(Qstring).ne.0) call Message(Qstring, frm = "(A$)")


if (PRESENT(frm)) then 
  call Message(out_string, frm = frm)
else
 call Message(out_string, frm = "(A)")
end if

end subroutine WriteValueString

! ###################################################################
! 
!  subroutine WriteValueIntShort 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write one or more short integers
!
!> @param Qstring question string
!> @param out_int output string
!> @param num optional number of integers
!> @param frm optional formatting argument
!
!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
recursive subroutine WriteValueIntShort(Qstring, out_int, num, frm)
!DEC$ ATTRIBUTES DLLEXPORT :: WriteValueIntShort

character(*), INTENT(IN)                        :: Qstring
integer(kind=ish),INTENT(IN)                    :: out_int(*)
character(*),INTENT(IN),OPTIONAL                :: frm
integer(kind=irg),INTENT(IN),OPTIONAL           :: num

if (len(Qstring).ne.0) call Message(Qstring, frm = "(A$)")

! one or more than one values expected ?
if (PRESENT(num)) then
 if (PRESENT(frm)) then
  write (stdout, fmt=frm) (out_int(i),i=1,num)
 else
  write (stdout,*) (out_int(i),i=1,num)
 end if
else
 if (PRESENT(frm)) then
  write (stdout, fmt=frm) out_int(1)
 else
  write (stdout,*) out_int(1)
 end if
end if

end subroutine WriteValueIntShort

! ###################################################################
! 
!  subroutine WriteValueIntLong 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write one or more 4-byte integers
!
!> @param Qstring question string
!> @param out_int output string
!> @param num optional number of integers
!> @param frm optional formatting argument
!
!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
recursive subroutine WriteValueIntLong(Qstring, out_int, num, frm)
!DEC$ ATTRIBUTES DLLEXPORT :: WriteValueIntLong

character(*), INTENT(IN)                        :: Qstring
integer(kind=irg),INTENT(IN)                    :: out_int(*)
character(*),INTENT(IN),OPTIONAL                :: frm
integer(kind=irg),INTENT(IN),OPTIONAL           :: num

if (len(Qstring).ne.0) call Message(Qstring, frm = "(A$)")

! one or more than one values expected ?
if (PRESENT(num)) then
 if (PRESENT(frm)) then
  write (stdout, fmt=frm) (out_int(i),i=1,num)
 else
  write (stdout,*) (out_int(i),i=1,num)
 end if
else
 if (PRESENT(frm)) then
  write (stdout, fmt=frm) out_int(1)
 else
  write (stdout,*) out_int(1)
 end if
end if

end subroutine WriteValueIntLong

! ###################################################################
! 
!  subroutine WriteValueIntLongLong
!
!> @author Saransh, Carnegie Mellon University
!
!> @brief write one or more 8-byte integers
!
!> @param Qstring question string
!> @param out_int output string
!> @param num optional number of integers
!> @param frm optional formatting argument
!
!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
recursive subroutine WriteValueIntLongLong(Qstring, out_int, num, frm)
!DEC$ ATTRIBUTES DLLEXPORT :: WriteValueIntLongLong

character(*), INTENT(IN)                        :: Qstring
integer(kind=ill),INTENT(IN)                    :: out_int(*)
character(*),INTENT(IN),OPTIONAL                :: frm
integer(kind=irg),INTENT(IN),OPTIONAL           :: num

if (len(Qstring).ne.0) call Message(Qstring, frm = "(A$)")

! one or more than one values expected ?
if (PRESENT(num)) then
 if (PRESENT(frm)) then
  write (stdout, fmt=frm) (out_int(i),i=1,num)
 else
  write (stdout,*) (out_int(i),i=1,num)
 end if
else
 if (PRESENT(frm)) then
  write (stdout, fmt=frm) out_int(1)
 else
  write (stdout,*) out_int(1)
 end if
end if

end subroutine WriteValueIntLongLong


! ###################################################################
! 
!  subroutine WriteValueRealSingle 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write one or more single precision reals
!
!> @param Qstring question string
!> @param out_real output string
!> @param num optional number of integers
!> @param frm optional formatting argument
!> @param stdout optional output unit identifier
!
!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
recursive subroutine WriteValueRealSingle(Qstring, out_real, num, frm)
!DEC$ ATTRIBUTES DLLEXPORT :: WriteValueRealSingle

character(*), INTENT(IN)                        :: Qstring
real(kind=sgl),INTENT(IN)                       :: out_real(*)
character(*),INTENT(IN),OPTIONAL                :: frm
integer(kind=irg),INTENT(IN),OPTIONAL           :: num

if (len(Qstring).ne.0) call Message(Qstring, frm = "(A$)")

! one or more than one values expected ?
if (PRESENT(num)) then
 if (PRESENT(frm)) then
  write (stdout, fmt=frm) (out_real(i),i=1,num)
 else
  write (stdout,*) (out_real(i),i=1,num)
 end if
else
 if (PRESENT(frm)) then
  write (stdout, fmt=frm) out_real(1)
 else
  write (stdout,*) out_real(1)
 end if
end if

end subroutine WriteValueRealSingle



! ###################################################################
! 
!  subroutine WriteValueRealDouble 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write one or more double precision reals
!
!> @param Qstring question string
!> @param out_real output string
!> @param num optional number of integers
!> @param frm optional formatting argument
!
!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
recursive subroutine WriteValueRealDouble(Qstring, out_real, num, frm)
!DEC$ ATTRIBUTES DLLEXPORT :: WriteValueRealDouble

character(*), INTENT(IN)                        :: Qstring
real(kind=dbl),INTENT(IN)                       :: out_real(*)
character(*),INTENT(IN),OPTIONAL                :: frm
integer(kind=irg),INTENT(IN),OPTIONAL           :: num

if (len(Qstring).ne.0) call Message(Qstring, frm = "(A$)")

! one or more than one values expected ?
if (PRESENT(num)) then
 if (PRESENT(frm)) then
  write (stdout, fmt=frm) (out_real(i),i=1,num)
 else
  write (stdout,*) (out_real(i),i=1,num)
 end if
else
 if (PRESENT(frm)) then
  write (stdout, fmt=frm) out_real(1)
 else
  write (stdout,*) out_real(1)
 end if
end if

end subroutine WriteValueRealDouble


! ###################################################################
! 
!  subroutine WriteValueRealComplex 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write one or more single precision complex numbers
!
!> @param Qstring question string
!> @param out_real output string
!> @param num optional number of integers
!> @param frm optional formatting argument
!
!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
recursive subroutine WriteValueRealComplex(Qstring, out_cmplx, num, frm)
!DEC$ ATTRIBUTES DLLEXPORT :: WriteValueRealComplex

character(*), INTENT(IN)                        :: Qstring
complex(kind=sgl),INTENT(IN)                    :: out_cmplx(*)
character(*),INTENT(IN),OPTIONAL                :: frm
integer(kind=irg),INTENT(IN),OPTIONAL           :: num

if (len(Qstring).ne.0) call Message(Qstring, frm = "(A$)")

! one or more than one values expected ?
if (PRESENT(num)) then
 if (PRESENT(frm)) then
  write (stdout, fmt=frm) (out_cmplx(i),i=1,num)
 else
  write (stdout,*) (out_cmplx(i),i=1,num)
 end if
else
 if (PRESENT(frm)) then
  write (stdout, fmt=frm) out_cmplx(1)
 else
  write (stdout,*) out_cmplx(1)
 end if
end if

end subroutine WriteValueRealComplex




recursive subroutine PrintMatrixd(s,a)
!DEC$ ATTRIBUTES DLLEXPORT :: PrintMatrixd

use local

IMPLICIT NONE

real(kind=dbl)   :: a(3,3)
integer(kind=irg):: i,j
character(4)     :: s

write (stdout,"(A/)") s
do i=1,3
  write (stdout,"(3(F12.5,2x))") (a(i,j),j=1,3)
end do
write (stdout,"(/)")

end subroutine

recursive subroutine PrintMatrixcd(s,a)
!DEC$ ATTRIBUTES DLLEXPORT :: PrintMatrixcd

use local

IMPLICIT NONE

complex(kind=dbl)   :: a(3,3)
integer(kind=irg):: i,j
character(4)     :: s

write (stdout,"(A/)") s
do i=1,3
  write (stdout,*) (a(i,j),j=1,3)
end do
write (stdout,"(/)")

end subroutine


end module

