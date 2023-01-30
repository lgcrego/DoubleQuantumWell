program QuantumPiston

use constants_and_parameters
use IO_routines    , only : create_timely_dir
use dinamic        , only : QW_push_pull
use estado_inicial , only : get_Psi_t0

implicit none

real*8       :: l22_maximo , l22_minimo , freq 

!==============================================
!       Escolha dos parametros do poco 

 call setup( l22_maximo , l22_minimo , freq )

!==============================================

call get_Psi_t0

call QW_push_pull( l22_maximo , l22_minimo , freq )

call create_timely_dir( instance = 'move' )

end program QuantumPiston
!
!
!
!==================================================
 subroutine setup( l22_maximo , l22_minimo , freq )
!==================================================

use constants_and_parameters
use griding        , only : grade
use rootfinding    , only : Roots, func_dwell
use IO_routines    , only : create_timely_dir

implicit none
real*8 , intent(out) :: l22_maximo , l22_minimo , freq 

!local variables 
real*8               :: l00 , l11 , l22 , v00 , v11 , de
real*8 , allocatable :: energias(:)
integer              :: regime 

!  Declarar Parametros
v00 = 250.d0
v11 = 250.d0
l00 = 45.d0                    ! Barreira
l11 = 100.d0                   ! Poço Esquerda

!    Admensionalisando
v0 = v00/RYMEV
v1 = v11/RYMEV
l1 = l11/ABA
l0 = l00/ABA

!===================================================================
!       Escolha dos parametros do poco para dinamica.
!===================================================================

!---------------------------------
!  Escolha 1 para ressonancia
regime = 1    ! 0 para caso adiabatico
              ! 1 para caso ressonante
              ! 2 para caso nonadiabatico

ciclo       = 2
nivel       = 1
N_all_roots = 40

!---------------------------------
!  Poço Direita (compressível)
l22_maximo = 150.d0 + 100.d0*float((nivel - 1))
l22_minimo = 50.d0
dl         = 0.5d0/float(nivel)

l22 = l11
l2  = l22/aba

call grade

if(.not. allocated( energias ) ) allocate( energias , source = Roots( func_dwell ) )

de = abs( energias(2*nivel-1) - energias(2*nivel) )

freq = de

deallocate( energias )
!---------------------------------

if( regime == 0 ) freq = freq/50.d0
if( regime == 1 ) freq = freq
if( regime == 2 ) freq = freq*50.d0

call create_timely_dir( instance = 'make' )

end subroutine setup
