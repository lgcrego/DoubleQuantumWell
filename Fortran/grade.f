module griding

use constants_and_parameters

implicit none

!   module variables 

    public :: grade , sumtrap

    contains

!==================
 subroutine Grade
!==================
! variaveis locais ...
integer               :: i , n , k
integer , allocatable :: domain(:)
real*8  , allocatable :: x1(:) , x2(:) 

if ( .not. allocated(x1)) allocate( x1(N_of_domains) )
if ( .not. allocated(x2)) allocate( x2(N_of_domains) )
!  |                                                        |
!  |                                                        |
!  |                                                        |
!  |_______________         _______          _______________|
!                 |         |     |         |               
!                 |         |     |         |              
!         1       |    2    |  3  |    4    |       5      
!                 |    l1   |  l0 |    l2   |              
!                 |________ |     |________ |              
! domain 1 ...
x1(1) = -(lim_left+l0+l1)
x2(1) = zero
p(1)  = nint((x2(1)-x1(1))*aba/dl)

! doamin 2 ...
x1(2) = zero
x2(2) = l1
p(2)  = nint((x2(2)-x1(2))*aba/dl)

! domain 3 ...
x1(3) = l1
x2(3) = (l1+l0)
p(3)  = nint((x2(3)-x1(3))*aba/dl)

! domain 4 ...
x1(4) = (l1+l0)
x2(4) = (l1 + l0 + l2)
p(4)  = nint((x2(4)-x1(4))*aba/dl)

! domain 5 ...
x1(5) = (l1+l0+l2)
x2(5) = lim_right+l0+l1
p(5)  = nint((x2(5)-x1(5))*aba/dl)

if (.not. allocated( x )) allocate( x(sum(p),2) )
if (.not. allocated( domain )) allocate( domain(sum(p)) )

! make dynamic grid ...
k = 1
do n = 1 , N_of_domains
    do i = 1 , p(n)
        x(k,1) = (x1(n)) + float(i)*(abs( x2(n) - x1(n) )) / float(p(n))
        x(k,2) = n

        domain(k) = n
        if( x(k,1) == (l0+l1+l2)) parameter_n = k
        k = k + 1
    enddo
enddo

grid_size = size(x(:,1)) 


end subroutine Grade

!
!===========================
 function sumtrap(i1,i2,x,y)
!===========================
integer , intent(in) :: i1 , i2
real*8  , intent(in) :: x(:)
real*8  , intent(in) :: y(:)

real*8  :: sumtrap

!------------------------------------------------------------------------------
! CALCULA A INTEGRAL DA FUNCAO Y(I) PELO METODO DO TRAPEZIO COM PASSO VARIAVEL
!------------------------------------------------------------------------------

sumtrap  = sum( (x(i1+1:i2)-x(i1:i2-1)) * (y(i1+1:i2)+y(i1:i2-1)) ) / two


end function sumtrap

end module griding
