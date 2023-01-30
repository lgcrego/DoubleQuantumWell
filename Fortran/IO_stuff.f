module IO_routines

use constants_and_parameters
use griding                  , only : grade , sumtrap 

implicit none

!   module variables 

    public :: create_timely_dir , diagram_p_x , GenerateGnuPlotScript , write_pop , write_energies , write_eigenenergies

    contains

!=================================
 Subroutine diagram_p_x( psi_t )
!=================================
implicit none
complex*16 , intent(in)  ::  psi_t(:)

! local variable
complex*16 , allocatable :: psi_tlinha(:)
real*8     , allocatable :: vetor_re_p(:)
real*8     , allocatable :: vetor_img_p(:)
real*8     , allocatable :: vetor_x(:)
real*8     :: medx , medp_re , a(grid_size) , medp_img 
integer    :: i , n , k 

forall( i = 1 : grid_size-1) a(i) = x(i+1,1) - x(i,1)
a(grid_size) = a(grid_size-1)

if ( .not. allocated(psi_tlinha))  allocate( psi_tlinha( grid_size ) )
if ( .not. allocated(vetor_re_p))   allocate( vetor_re_p( grid_size ) )
if ( .not. allocated(vetor_img_p))  allocate( vetor_img_p( grid_size ) )
if ( .not. allocated(vetor_x))      allocate( vetor_x( grid_size ) )

k = 1
do n = 1 , N_of_domains
    do i = 1 , p(n)
        psi_tlinha( k ) = -zi*(psi_t( k + 1 ) - psi_t( k ))/a(k)  
        k = k + 1
    enddo
enddo

vetor_re_p(:)  = real( conjg( psi_t(:) )*psi_tlinha(:) )
medp_re   = sumtrap( 1 , grid_size , x(:,1) , vetor_re_p )

vetor_img_p(:) = imag( conjg( psi_t(:) )*psi_tlinha(:) )
medp_img  = sumtrap( 1 , grid_size , x(:,1) , vetor_img_p )


forall( i = 1 : grid_size-1 ) vetor_x(i) = conjg(psi_t(i))*x(i,1)*aba*psi_t(i)
medx = sumtrap( 1 , grid_size , x(:,1) , vetor_x )

open( 60 , file='diagramox-real.dat' , position='append' )
write( 60 , 13 ) medx , medp_re
close( 60 )

open( 61 , file='diagrampx-img.dat' , position='append' )
write( 61 , 13 ) medx , medp_img
close( 61 )

deallocate(psi_tlinha,vetor_re_p,vetor_img_p,vetor_x)

13 format(7d15.4)

end subroutine diagram_p_x
!
!
!
!=====================================================
subroutine GenerateGnuPlotScript( psi_t , energy_wf )
!=====================================================
implicit none
complex*16 , intent(in) ::  psi_t(:)
real*8     , intent(in) :: energy_wf
integer :: i

!========================================================================
! Não esquecer de descomentar a linha 113 para fazer o video
!========================================================================
! ESCREVE AS FUNÇÕES QUE SERÃO DESENHADAS NO GNUPLOT.
!========================================================================                                                                                                                                       
open(20 , file="plot.gnu" )
 write(20,fmt='(a)')        'clear'
 write(20,fmt='(a)')        'reset'
 write(20,fmt='(a)')        'unset key'
 write(20,fmt='(a)')        "set title'Dinâmica do Pacote'"
 write(20,fmt='(a)')        "set terminal pngcairo size 1000,600 enhanced font 'Verdana,20'"
 write(20,fmt='(a,f6.1,a,f6.1,a)') "set xrange[-100:350]"
 write(20,fmt='(a)')        "set yrange[0:300]"
 write(20,fmt='(a)')        "set y2range[0:300]"
 write(20,fmt='(a)')        "set xlabel'Comprimento'"
 write(20,fmt='(a)')        "set y2label'Energy'"
 write(20,fmt='(a)')        "set ylabel '20*{/Symbol y} + E_{package}' textcolor rgb 'blue'"
 write(20,fmt='(a)')        "set xtics 100"
!\/ ja estava comentado
! write(20,fmt='(a)')        'set ytics 300 textcolor rgb "blue"'
!/\ ja estava comentado
! write(20,fmt='(a)')        "set y2tics 250"
 write(20,fmt='(a)')        "set pointsize 0.005"
 write(20,fmt='(a)')        'set border linewidth 1.0'
 write(20,fmt='(a)')        "set samples 300000"
 write(20,fmt='(a,i0,a)')        "system('mkdir -p animacao')"


!\/ ja estava comentado
!    write(*,*)  '   Frame = ' , tlinha
!/\ ja estava comentado
!    write(20,fmt='(a,i0,a,i0,a)') "outfile = sprintf('animacao/geral%05.0f.png',",step,")"
    write(20,fmt='(a)')      "set output outfile"
    write(20,fmt='(a,i3,a,i3,a)') "plot '-' u 1:2 w l lw 4 lc 3 , '-' u 1:4 w l lw 4 lc 0"

!=================================================
do i = 1 , grid_size
    write( 20 , 14 ) x(i,1)*aba , 20*conjg(psi_t(i))*psi_t(i) + energy_wf 
enddo

 write(20,fmt='(a)') "unset output outfile"
 write(20,fmt='(a)') 'e'
 write(20,fmt='(a)') ' '

14 format(7e15.4)

end subroutine GenerateGnuPlotScript
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
