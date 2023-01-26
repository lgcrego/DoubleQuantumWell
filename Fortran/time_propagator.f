module propagator

use f95_precision
use blas95
use lapack95
use constants_and_parameters
use RootFinding    , only : Roots   
use griding        , only : grade , sumtrap
use wfunction      , only : wavefunc , coefic
use wfunction_cont , only : wavefunc_cont , coefic_cont
use IO_routines    , only : GenerateGnuPlotScript , diagram_p_x , write_pop , write_energies
use Joule          , only : ensemble

implicit none

public :: quantum_propagation 

private

! module variables
complex*16 , allocatable :: psi_t(:) , temporal_operator(:)
real*8     , allocatable :: aux(:) , mat_rho(:)


contains

!============================================================
 Subroutine quantum_propagation( enp , step , tempo , sinal )
!============================================================
real*8  , intent(in)  :: enp(:)
real*8  , intent(in)  :: tempo , sinal
integer , intent(in)  :: step 

! local variable
real*8  :: norm , energy_wf
integer :: n , i 

call allocate_stuff

!================================
! setting up initial wave packet 
!================================
if( step == 1 .or. step == n_dt+1 ) then

        do n = 1 , N_all_roots
           aux(:)           = phi_adiabatic(:,n)*psi(:,nivel)
           coefi_phi_t(n,1) = sumtrap( 1 , grid_size , x(:,1) , aux )  ! <== C_n(t0)
        enddo
        
        energy_t0 = rymev*sum( enp(:)*coefi_phi_t(:,1)*conjg(coefi_phi_t(:,1)))
        
        !================================
        ! Funcao na base phi (adiabatica)
        ! psi_t0(i) = sum_n  [C_n(t0)*phi_n(t0)]
        !================================
        forall( i = 1 : grid_size )  psi_t(i) = sum( coefi_phi_t(:,1)*phi_adiabatic(i,:) )
        
        !========================
        ! normalization of psi_t
        !========================
        aux(:)           = conjg(psi_t(:))*psi_t(:)
        norm             = sumtrap( 1 , grid_size , x(:,1) , aux )
        coefi_phi_t(:,1) = coefi_phi_t(:,1)/sqrt(norm)
        
        forall( i = 1 : grid_size )  psi_t(i) = sum( coefi_phi_t(:,1)*phi_adiabatic(i,:) )
        
        coefi_phi_t(:,2) = coefi_phi_t(:,1)
        
        do n = 1 , N_all_Roots
           mat_rho(n) = coefi_phi_t(n,1)*conjg(coefi_phi_t(n,1))  ! <== trace = 1
        enddo

        if( MecStat ) call ensemble( mat_rho , tempo , enp , step , coefi_phi_t ,sinal )
        
        if( diagram ) call diagram_p_x( psi_t )

        if( verbose ) call write_pop(tempo , aux)
        
        !===========================
        ! Mudança de Base.  phi -> x
        !===========================
        forall( i = 1 : grid_size ) coefi_x_t(i,1) = sum(coefi_phi_t(:,1)*phi_adiabatic(i,:))
   
        deallocate( psi )

end if


!========================
! propagating wave packet 
!========================
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
        ! psi_t(i) = sum_n  [C_n(t+dt)*phi_n(t)]
        !==========================
        forall( i = 1 : grid_size ) psi_t(i) = sum( coefi_phi_t(:,2)*phi_adiabatic(i,:) )
  
        !========================
        ! normalization of psi_t
        !========================
        aux(:)         = conjg(psi_t(:))*psi_t(:)
        norm             = sumtrap( 1 , grid_size , x(:,1) , aux )
        coefi_phi_t(:,2) = coefi_phi_t(:,2)/sqrt(norm)
        
        forall( i = 1 : grid_size )  psi_t(i) = sum( coefi_phi_t(:,2)*phi_adiabatic(i,:) )
        
        do n = 1 , N_all_Roots
           mat_rho(n) = coefi_phi_t(n,2)*conjg(coefi_phi_t(n,2))
        enddo

        if( MecStat ) call ensemble( mat_rho , tempo , enp , step , coefi_phi_t , sinal )

        deallocate( mat_rho )

        !======================
        ! Calculo da População.
        !======================
        aux(:) = psi_t(:)*conjg(psi_t(:))

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

        if( diagram ) call diagram_p_x( psi_t )

end if

energy_wf = rymev*sum( enp(:)*coefi_phi_t(:,1)*conjg(coefi_phi_t(:,1)))

call write_energies( tempo, energy_wf  , energy_t0 )

! mencoder mf://*.png -mf w=1000:h=600:fps=20:type=png -ovc copy -oac copy -o din_wave.avi
if( video ) call GenerateGnuPlotScript( psi_t , energy_wf ) 

end subroutine
!
!
!
!===========================
 subroutine allocate_stuff
!===========================
implicit none

if ( .not. allocated(coefi_phi_t))        allocate( coefi_phi_t( N_all_roots , 2 ) , source = (0.d0,0.d0) )
if ( .not. allocated(coefi_x_t))          allocate( coefi_x_t( grid_size , 2 ) )
if ( .not. allocated(psi_t))              allocate( psi_t( grid_size ) )
if ( .not. allocated(aux))                allocate( aux( grid_size ) )
if ( .not. allocated(temporal_operator) ) allocate( temporal_operator( N_all_roots ))
if ( .not. allocated(mat_rho) )           allocate( mat_rho( N_all_roots ) , source=0.d0)

end subroutine allocate_stuff
!
!
!
end module
