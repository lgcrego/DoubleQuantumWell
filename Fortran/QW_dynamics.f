module QW_dyn

use constants_and_parameters
use rootfinding    , only : Roots , func_dwell , func_cont_dwell 
use griding        , only : grade
use wfunction      , only : coefic
use wfunction_cont , only : coefic_cont
use propagator     , only : quantum_propagation
use IO_routines    , only : write_eigenenergies

implicit none 

public :: QW_dynamics

contains

!===============================================
 Subroutine QW_dynamics( step , tempo , sinal )
!===============================================
integer , intent(in) :: step 
real*8  , intent(in) :: tempo 
real*8  , intent(in) :: sinal

! local variables
real*8  , allocatable :: energias_cont(:) , energias_dw(:) , energias(:)

!=====================================
!       Calculo tamanho da base
!=====================================
call grade

if(.not. allocated( energias ))         allocate( energias( N_all_roots ) , source = 0.d0 )
if(.not. allocated( phi_adiabatic ))    allocate( phi_adiabatic( grid_size , N_all_roots ) , source = 0.d0 )

!=======================================
!       Poço duplo abaixo de v0
!=======================================
if(.not. allocated( energias_dw ) )     allocate( energias_dw , source = Roots( Func_dwell ))
call coefic( energias_dw , N_of_roots )

energias( 1 : N_of_roots ) = energias_dw(:)  !<== bound states
deallocate( energias_dw )

!=======================================
!       Poço duplo acima de v0
!=======================================
if(.not. allocated( energias_cont ) )   allocate( energias_cont , source = Roots( Func_cont_dwell ))
call coefic_cont( energias_cont , N_of_roots_cont )

energias( N_of_roots+1 : N_all_roots ) = energias_cont(:)  !<== all (bound+continuum) states
deallocate( energias_cont )

!====================================
!       Cálculo da Dinâmica
!====================================
call quantum_propagation( energias , step , tempo , sinal )

call write_eigenenergies( tempo , energias )  

deallocate( energias , x , phi_adiabatic )

end subroutine QW_dynamics

end module
