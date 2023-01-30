module estado_inicial

use constants_and_parameters
use griding                  , only : grade , sumtrap
use RootFinding              , only : Roots , func_swell

implicit none

!   module variables 

    public :: get_Psi_t0

    contains

!====================
subroutine get_Psi_t0
!====================
implicit none

!local variables
real*8 , allocatable :: energias_sw(:)

l2 = 0.d0

call grade

if(.not. allocated( energias_sw ))  allocate( energias_sw , source = Roots( Func_swell ) )

!=====================================================================
! calculates the wavefunction of the initial state 
! Estado Inicial localizado no nivel=nivel do poço simples da esquerda

call Coeficientes( energias_sw , nivel )
!=====================================================================

! single well not necessary anymore
deallocate( energias_sw )

end subroutine get_Psi_t0
!
!
!
!============================================
 subroutine Coeficientes( e_n_swell , root ) 
! Estado Inicial localizado no nivel=root 
! e_n_swell(root) = energia do nivel=root 
!============================================
real*8  , intent(in) :: e_n_swell(:)
integer , intent(in) :: root

!   variaveis locais
real*8 :: A,B,C,D,An,Bn,Cn,Dn,psinormalizado,c1,c2
real*8 :: k2,k1,k2l1,k1l1,Am(500),Bm(500),Cm(500),Dm(500)
real*8 , allocatable :: psi_2(:) 

 k1 = sqrt( v0 - e_n_swell(root) )
 k2 = sqrt( e_n_swell(root) )

 k1l1 = k1*l1
 k2l1 = k2*l1

 c1 = ( k2*cos(k2l1)/k1 + sin(k2l1) )
 c2 = ( k2*( -k2*sin(k2l1)/k1 + cos(k2l1) ))

 D = 1
 C = D*exp(-k1l1)/c1
 B = k2*C/k1
 A = B       

 if( .not. allocated(psi_2) ) allocate( psi_2( grid_size ) )
 
 where( x(:,2) == 1 ) psi_2 = (A*exp( k1*x(:,1) ) )**2
 where( x(:,2) == 2 ) psi_2 = (B*cos( k2*x(:,1)) + C*sin(k2*x(:,1)) )**2
 where( x(:,2) == 3 ) psi_2 = (D*exp(-k1*x(:,1)) )**2
 where( x(:,2) == 4 ) psi_2 = (D*exp(-k1*x(:,1)) )**2
 where( x(:,2) == 5 ) psi_2 = (D*exp(-k1*x(:,1)) )**2

 psinormalizado = sumtrap( 1 , grid_size , x(:,1) , psi_2 )

 psinormalizado = one/sqrt( psinormalizado )

 Dn = psinormalizado*D
 Cn = psinormalizado*C
 Bn = psinormalizado*B
 An = psinormalizado*A

 Dm(root) = Dn
 Cm(root) = Cn
 Bm(root) = Bn
 Am(root) = An

! wave-function do estado inicial ...
 call  wf_t0( e_n_swell , root , An , Bn , Cn , Dn )

end subroutine Coeficientes
!
!
!
!====================================================
 subroutine wf_t0( e_n_swell , root , A , B , C , D )
!====================================================
real*8  , intent(in) :: A , B , C , D
real*8  , intent(in) :: e_n_swell(:)
integer , intent(in) :: root

!   variveis locais
real*8  :: rk1 , rk2
real*8  :: norm_psi
integer :: i
real*8  , allocatable   :: Psi_2(:) 

 if( .not. allocated(Psi) )   allocate( Psi( size(x(:,1)) , N_of_roots ) , source = 0.d0 )
 if( .not. allocated(Psi_2) ) allocate( Psi_2(grid_size) )

 rk1 = sqrt( v0 - e_n_swell(root) )
 rk2 = sqrt( e_n_swell(root) )

 where( x(:,2) == 1 ) Psi(:,root) = A*exp( rk1*x(:,1) )
 where( x(:,2) == 2 ) Psi(:,root) = B*cos( rk2*x(:,1) ) + C*sin( rk2*x(:,1) )
 where( x(:,2) == 3 ) Psi(:,root) = D*exp(-rk1*x(:,1))
 where( x(:,2) == 4 ) Psi(:,root) = D*exp(-rk1*x(:,1))
 where( x(:,2) == 5 ) Psi(:,root) = D*exp(-rk1*x(:,1))

 Psi_2 = Psi(:,root) * Psi(:,root)
 
 norm_psi = sumtrap( 1 , grid_size , x(:,1) , Psi_2 )

 Psi(:,root) = Psi(:,root) / sqrt(norm_psi)

if( verbose ) then
  do i = 1 , grid_size
      write(100 + root , 13) x(i,1) , Psi(i,root)
  end do
  endif

13  format(7e15.4)

end subroutine wf_t0
!
!
!
end module estado_inicial
