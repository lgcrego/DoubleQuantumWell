module dinamic

use constants_and_parameters
use rootfinding    , only : Roots , func_dwell , func_cont_dwell , func_swell
use griding        , only : grade
use g_state        , only : Coeficientes
use wfunction      , only : coefic
use wfunction_cont , only : coefic_cont
use basis          , only : dinamic_eletron
use IO_routines    , only : write_eigenenergies

implicit none 

public :: dinamic_well

contains

!===============================================================================
 Subroutine dinamic_well( n_all_roots , l22_maximo , l22_minimo , freq , ciclo )
!===============================================================================
integer , intent(in) :: n_all_roots , ciclo 
real*8  , intent(in) :: l22_maximo , l22_minimo , freq

! local variables
integer :: passo , i , j , m , k
real*8  :: ti , Ac , Bc , l22 , tf , tempo , sinal
real*8  , allocatable :: energias_sw(:) , energias_cont(:) , energias_dw(:) , energias(:)

!=====================================================================
! Estado Inicial localizado no nivel=nivel do poço simples da esquerda
!=====================================================================
l2 = 0.d0

call grade

if(.not. allocated( energias_sw ))  allocate( energias_sw , source = Roots( Func_swell ) )

! calculates the wavefunction of the initial state = wf_gstate
call Coeficientes( energias_sw , nivel )

! single well not necessary anymore
deallocate( energias_sw )

!=================================
!       Expandindo o poço
!=================================
Ac = (l22_maximo+l22_minimo)/2.d0
Bc = (l22_maximo-l22_minimo)/2.d0

tempo = 0.d0
passo = 0
ti    = 0.d0
n_dt  = ciclo*( nint(200.d0/dl) + 1 )

do k = 1 , 4
do m = 1 , ciclo

   do i = 0 , nint( 100.d0 / dl )

        passo = passo + 1
        tf  = one*( acos( 1.d0 - (dl*float(i))/50.d0 ))/(freq) 
        l22 = Ac - Bc*cos( tf * freq )
        l2  = l22/ABA
   
        dt = abs( tf - ti )
        sinal = 1.d0
        tempo = tempo + dt*freq

        print*,  "width of second well is " , l22 , "time is ", tempo

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
        call dinamic_eletron( energias , passo , tempo , sinal , k )
        
        call write_eigenenergies( tempo , energias )  

        deallocate( energias , x , phi_adiabatic )

        ti = tf
   enddo

   !==========================
   !   Comprimindo o Poço
   !==========================

   do j = nint( 100.d0 / dl )  , nint( 200.d0 / dl ) 

        passo = passo + 1
   
        tf = one*( acos( (dl*float(j))/50.d0 - 3.d0 ) )/(freq)
        l22 = Ac - Bc*cos( tf * freq )
        l2  = l22/ABA
   
        dt = abs( tf - ti )
        sinal = -1.d0
        tempo = tempo + dt*freq
        
        print*,  "width of second well is " , l22 , "tempo is ", tempo

        !======================================
        !       Calculo do tamanho da base.
        !======================================
        call grade
   
        !===================================================
        !       Energia abaixo de v0 e acima de v0
        !===================================================
        if(.not. allocated( energias ))         allocate( energias( N_all_roots ) )
        if(.not. allocated( phi_adiabatic ))    allocate( phi_adiabatic( grid_size , N_all_roots ) , source = 0.d0 )
        !===================================================
        !       Poço Duplo abaixo de v0
        !===================================================
        if(.not. allocated( energias_dw ) )     allocate( energias_dw     , source = Roots( func_dwell )   )
            call coefic( energias_dw , N_of_roots )
            energias(1:N_of_roots) = energias_dw
        deallocate( energias_dw )
        !===================================================
        ! Poço Duplo acima de v0
        !===================================================
        if(.not. allocated( energias_cont ) )   allocate( energias_cont  , source = Roots( func_cont_dwell )       )
            call coefic_cont( energias_cont , N_of_roots_cont )
            energias(N_of_roots+1 : N_all_roots ) = energias_cont
        deallocate( energias_cont )
        !===================================================
        ! Cálculo da Dinâmica
        !===================================================
        call dinamic_eletron( energias , passo , tempo , sinal , k)
   
        call write_eigenenergies( tempo , energias )  

        deallocate( energias , x , phi_adiabatic )

        ti = tf

   enddo

enddo
enddo

13 format(7e15.4)

end subroutine

end module
