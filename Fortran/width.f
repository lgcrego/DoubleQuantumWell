module dinamic

use constants_and_parameters
use QW_dyn , only : QW_dynamics

implicit none 

public :: QW_push_pull

contains

!=========================================================
 Subroutine QW_push_pull( l22_maximo , l22_minimo , freq )
!=========================================================
real*8  , intent(in) :: l22_maximo , l22_minimo , freq

! local variables
integer :: step , i , m 
real*8  :: ti , Ac , Bc , l22 , tf , tempo , sinal

Ac = (l22_maximo+l22_minimo) / two
Bc = (l22_maximo-l22_minimo) / two

tempo = 0.d0
step = 0
ti    = 0.d0
n_dt  = ciclo*( nint(200.d0/dl) + 1 )

do m = 1 , ciclo

   !=================================
   !       Expandindo o poço
   !=================================
   do i = 0 , nint( 100.d0 / dl )

        step = step + 1
        tf  = acos(one - dl*float(i)/50.d0) / freq 
        l22 = Ac - Bc*cos(tf * freq)
        l2  = l22/ABA
   
        dt = abs(tf - ti)
        sinal = 1.d0
        tempo = tempo + dt*freq

        print*,  "width of second well is " , l22 , "time is ", tempo

        call QW_dynamics( step , tempo , sinal ) 

        ti = tf
   enddo

   !==========================
   !   Comprimindo o Poço
   !==========================
   do i = nint( 100.d0 / dl )  , nint( 200.d0 / dl ) 

        step = step + 1
   
        tf = acos(dl*float(i)/50.d0 - three) / freq
        l22 = Ac - Bc*cos(tf * freq)
        l2  = l22/ABA
   
        dt = abs(tf - ti)
        sinal = -1.d0
        tempo = tempo + dt*freq
        
        print*,  "width of second well is " , l22 , "tempo is ", tempo

        call QW_dynamics( step , tempo , sinal ) 

        ti = tf
   enddo

enddo

end subroutine

end module
