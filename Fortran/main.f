program dwell

use constants_and_parameters
use IO_routines              , only : create_timely_dir
use rootfinding              , only : Roots,&
                                      func_dwell,&
                                      func_swell,&
                                      func_cont_dwell
use griding                  , only : grade
use g_state                  , only : Coeficientes
use wfunction                , only : coefic
use wfunction_cont           , only : coefic_cont
use basis                    , only : dinamic_eletron
use dinamic                  , only : dinamic_well

implicit none

real*8 , allocatable :: energias(:)
real*8       :: l00 , l11 , l22 , v00 , v11 , freq , l22_maximo , l22_minimo , de
integer      :: resso 
character(4) :: instance

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
resso = 1     ! 0 para caso adiabatico
              ! 1 para caso ressonante
              ! 2 para caso nonadiabatico

ciclo       = 1
nivel       = 1
N_all_roots = 40

!---------------------------------
!       Poço Direita
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

if( resso == 0 ) freq = freq/50.d0
if( resso == 1 ) freq = freq
if( resso == 2 ) freq = freq*50.d0

if( .false. ) call create_timely_dir( instance = 'make' )

call dinamic_well( n_all_roots, l22_maximo , l22_minimo , freq , ciclo )

if( .false. ) call create_timely_dir( instance = 'move' )

end program
