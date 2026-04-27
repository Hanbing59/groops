!***********************************************/
!**
!* @file dtm2020Wrapper.F90
!*
!* @brief Fortran Wrapper.
!*
!* @author Andreas Kvas
!* @date   2022-04-21
!*
!/***********************************************/

SUBROUTINE dtm2020initWrapper(inputFile) bind(C, name="dtm2020initWrapper")

  use, intrinsic                              :: ISO_C_BINDING

  IMPLICIT NONE
  character(len=1,kind=C_char), intent(in)    :: inputFile(*)
  character(len=:), allocatable               :: str1
  integer                                     :: i, nchars

  ! CONVERT: c_char -> char necessary to call the msisinit subroutine
  i = 1
  do
     if (inputFile(i) == c_null_char) exit
     i = i + 1
  end do
  nchars = i - 1  ! Exclude null character from Fortran string
  allocate(character(len=nchars) :: str1)
  str1 = transfer(inputFile(1:nchars), str1)

  ! initialise the model
  open(unit=42, file=str1)
  call lecdtm(42)
  close(42)

  deallocate(str1)

END SUBROUTINE dtm2020initWrapper

!**********************************************

SUBROUTINE dtm2020calcWrapper(day,f,fbar,akp,akpMean,alti,hl,alat,xlon,tz,tinf,ro,d,wmm) bind(C, name="dtm2020calcWrapper")

  use, intrinsic       :: ISO_C_BINDING
  real, dimension(2) :: f_input, fbar_input
  real, dimension(4) :: akp_input

  f_input(1) = f
  f_input(2) = 0

  fbar_input(1) = fbar
  fbar_input(2) = 0

  akp_input(1) = akp
  akp_input(2) = 0
  akp_input(3) = akpMean
  akp_input(4) = 0

  call dtm3(day,f_input,fbar_input,akp_input,alti,hl,alat,xlon,tz,tinf,ro,d,wmm)

END SUBROUTINE dtm2020calcWrapper

!**********************************************
