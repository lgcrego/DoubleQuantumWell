module basis
use f95_precision
use blas95
use lapack95

use constants_and_parameters
use RootFinding    , only : Roots   
use griding        , only : grade , sumtrap
use g_state        , only : wf_gstate , coeficientes
use wfunction      , only : wavefunc , coefic
use wfunction_cont , only : wavefunc_cont , coefic_cont
use IO_routines    , only : GenerateGnuPlotScript , diagram_p_x , write_pop , write_energies
use Joule          , only : ensemble

implicit none

public :: dinamic_eletron 

contains

!==============================================================
 Subroutine dinamic_eletron( enp , step , tempo , sinal , level )
!==============================================================
real*8     , intent(in)  :: enp(:)
real*8     , intent(in)  :: tempo , sinal
integer    , intent(in)  :: step , level

! local variable
complex*16 , allocatable :: wave_t(:) , temporal_operator(:)
real*8     , allocatable :: aux(:) , mat_rho(:)
real*8                   :: norm , energy_wf
integer                  :: n , i 

if ( .not. allocated(coefi_phi_t))        allocate( coefi_phi_t( N_all_roots , 2 ) , source = (0.d0,0.d0) )
if ( .not. allocated(coefi_x_t))          allocate( coefi_x_t( grid_size , 2 ) )
if ( .not. allocated(wave_t))             allocate( wave_t( grid_size ) )
if ( .not. allocated(aux))                allocate( aux( grid_size ) )
if ( .not. allocated(temporal_operator) ) allocate(temporal_operator( N_all_roots ))
if ( .not. allocated(mat_rho) )           allocate(mat_rho( N_all_roots ) , source=0.d0)

if( step == 1 .or. step == n_dt+1 ) then

        do n = 1 , N_all_roots
           aux(:)               = phi_adiabatic(:,n)*psi(:,nivel)
           coefi_phi_t(n,1)     = sumtrap( 1 , grid_size , x(:,1) , aux )  ! <== C_n(t)
        enddo
        
        energy_zero_package = rymev*sum( enp(:)*coefi_phi_t(:,1)*conjg(coefi_phi_t(:,1)))
        
        !================================
        ! Funcao na base phi (adiabatica)
        ! wave_t(i) = sum_n  [C_n(t)*phi_n(t)]
        !================================
        forall( i = 1 : grid_size )  wave_t(i) = sum( coefi_phi_t(:,1)*phi_adiabatic(i,:) )
        
        !========================
        ! normalization of wave_t
        !========================
        aux(:)           = conjg(wave_t(:))*wave_t(:)
        norm             = sumtrap( 1 , grid_size , x(:,1) , aux )
        coefi_phi_t(:,1) = coefi_phi_t(:,1)/sqrt(norm)
        
        forall( i = 1 : grid_size )  wave_t(i) = sum( coefi_phi_t(:,1)*phi_adiabatic(i,:) )
        
        do n = 1 , N_all_Roots
           mat_rho(n) = coefi_phi_t(n,1)*conjg(coefi_phi_t(n,1))  ! <== trace = 1
        enddo

        coefi_phi_t(:,2) = coefi_phi_t(:,1)
        
        if( MecStat ) call ensemble( mat_rho , tempo , enp , step , coefi_phi_t ,sinal , level)
        
        if( diagram ) call diagram_p_x( wave_t )

        if( verbose ) call write_pop(tempo , aux)
        
        !===========================
        ! Mudança de Base.  phi -> x
        !===========================
        forall( i = 1 : grid_size ) coefi_x_t(i,1) = sum(coefi_phi_t(:,1)*phi_adiabatic(i,:))
   
        deallocate( psi )

end if


if( step > 1 ) then
        !===========================
        ! Mudança de Base.  x -> phi
        !===========================
        forall( n = 1 : N_all_roots ) coefi_phi_t(n,2) = sum(coefi_x_t(:,1)*phi_adiabatic(:,n))

        !===================================
        ! Evolução temporal (adiabatic part)
        !===================================
        do n = 1 , N_all_roots
           temporal_operator(n) = exp(-zi*enp(n)*dt)
           coefi_phi_t(n,2)     = temporal_operator(n)*coefi_phi_t(n,2)
        enddo

        !==========================
        ! Funcao na base phi
        ! wave_t(i) = sum_n  [C_n(t+dt)*phi_n(t)]
        !==========================
        forall( i = 1 : grid_size ) wave_t(i) = sum( coefi_phi_t(:,2)*phi_adiabatic(i,:) )
  
        !========================
        ! normalization of wave_t
        !========================
        aux(:)         = conjg(wave_t(:))*wave_t(:)
        norm             = sumtrap( 1 , grid_size , x(:,1) , aux )
        coefi_phi_t(:,2) = coefi_phi_t(:,2)/sqrt(norm)
        
        forall( i = 1 : grid_size )  wave_t(i) = sum( coefi_phi_t(:,2)*phi_adiabatic(i,:) )
        
        do n = 1 , N_all_Roots
           mat_rho(n) = coefi_phi_t(n,2)*conjg(coefi_phi_t(n,2))
        enddo

        if( MecStat ) call ensemble( mat_rho , tempo , enp , step , coefi_phi_t , sinal , level )

        deallocate( mat_rho )

        !======================
        ! Calculo da População.
        !======================
        aux(:) = wave_t(:)*conjg(wave_t(:))

        if( verbose ) call write_pop(tempo , aux)

        !==========================
        ! Mudança de Base  phi -> x
        !==========================
        forall(i=1:grid_size)  coefi_x_t(i,2) = sum(coefi_phi_t(:,2)*phi_adiabatic(i,:))
    
        !if( sinal == -1.d0 .and. l2 == 150.d0/aba ) then
        !        coefi_x_t(:,1)          = coefi_x_t(:,2)
        !print*, "deu certo?"
        !        coefi_phi_t(:,1)        = conjg(coefi_phi_t(:,2))
        !else
        coefi_x_t(:,1)   = coefi_x_t(:,2)
        coefi_phi_t(:,1) = coefi_phi_t(:,2)
        !end if

        if( diagram ) call diagram_p_x( wave_t )

end if

energy_wf = rymev*sum( enp(:)*coefi_phi_t(:,1)*conjg(coefi_phi_t(:,1)))

call write_energies( tempo, energy_wf  , energy_zero_package )

! mencoder mf://*.png -mf w=1000:h=600:fps=20:type=png -ovc copy -oac copy -o din_wave.avi
if( video ) call GenerateGnuPlotScript( wave_t , energy_wf ) 

end subroutine

end module
