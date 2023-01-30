module IO_routines

use constants_and_parameters
use griding                  , only : grade , sumtrap 

implicit none

!   module variables 

    public :: create_timely_dir , write_pop , write_energies , write_eigenenergies

    contains

!
!
!
!==================================
 subroutine write_pop( tempo, aux ) 
!==================================
implicit none
real*8  , intent(in) :: tempo , aux(:)

!local variables
real*8 :: pop1, pop2

pop1 = sumtrap( 1 , sum(p(1:2)) + p(3)/2 , x(:,1) ,  aux(:))            ! Psi square for first region
open( 365 , file='pop1.dat' , position='append' )
    write(365 , 13 ) tempo , pop1
close( 365 )

pop2 = sumtrap( sum(p(1:2)) + p(3)/2 , grid_size , x(:,1) , aux(:))     ! Psi square for second region
open( 366 , file='pop2.dat' , position='append' )
    write(366 , 13 ) tempo , pop2
close( 366 )

open( 367 , file='poptotal.dat' , position='append' )
    write(367 , 13 ) tempo , pop1+pop2
close( 367 )

13 format(7d15.4)

end subroutine write_pop
!
!
!
!====================================================================
 subroutine write_energies( tempo, energy_wf  , energy_t0 ) 
!====================================================================
implicit none
real*8  , intent(in) :: tempo , energy_wf  , energy_t0

open( 327 , file='energy_wf.dat' , position='append' )
    write(327 , 13 ) tempo , energy_wf , (energy_wf - energy_t0)
close( 327 )

13 format(3d15.4)

end subroutine write_energies
!
!
!
!==================================================
 subroutine write_eigenenergies( tempo , energias ) 
!==================================================
implicit none
real*8  , intent(in) :: tempo , energias(:)
open( 36 , file='eigen-energy.dat' , position='append' )
    write( 36 , 13 ) tempo , rymev*energias( 1 ) ,  rymev*energias( 2 ) , rymev*energias( 3 ), rymev*energias( 4 ) , rymev*energias( 5 ) , rymev*energias( 6 )
close( 36 )

13 format(7d15.4)

end subroutine write_eigenenergies
!
!
!

!========================================
 subroutine create_timely_dir( instance ) 
!========================================
character(len=*) , intent(in) :: instance

! variaveis locais ...
integer       :: date_time(8)
character(2)  :: day, hour, minute
character(4)  :: year 
character(17) , save :: nome 

! local parameters ...
character(len=3) :: month(12)=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
 
select case (instance)

   case('make') 
         call date_and_time(values=date_time) 
         
         write(day, '(i2)') date_time(3)
         write(year,'(i4)') date_time(1)
         write(hour,'(i2)') date_time(5)
         write(minute,'(i2)') date_time(6)
         
         nome = trim(day//'-'//month(date_time(2))//'-'//year//'-'//hour//':'//minute)

         call system('mkdir ' //  nome )

   case('move')
         call system('mv *.dat ' // nome )

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
