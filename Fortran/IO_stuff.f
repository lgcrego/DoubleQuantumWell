module IO_routines

use constants_and_parameters

implicit none

!   module variables 

    public :: create_timely_dir

    contains

!========================================
 subroutine create_timely_dir( instance ) 
!========================================
character(len=*) , intent(in) :: instance

! variaveis locais ...
integer      :: datatoday
character(10):: time
character(8) :: date
character(5) :: zone 

select case (instance)

   case('make') 
         call date_and_time( date=date , time=time , zone=zone , values=values )
         call date_and_time( DATE=date,ZONE=zone )
         call date_and_time( TIME=time )
         call date_and_time( VALUES=values )
         
         read( date , '(I8)' ) datatoday
         
         call system('mkdir ' // adjustl(trim( dirname( datatoday ) ) ))

   case('move')
         call system('mv *.dat ' // adjustl(trim( dirname( datatoday ) ) ))

end select

end subroutine create_timely_dir
!
!
!
!========================
 function dirname(number)
!========================
    integer,intent(in) :: number
    character(len=8)   :: dirname

    ! Cast the (rounded) number to string using 6 digits and
    ! leading zeros
    ! This is the same w/o leading zeros  
    write (dirname, '(I8)')  number

    ! This is for one digit (no rounding)
    !write (dirname, '(F4.1)')  number

end function
!
!
!
end module IO_routines
