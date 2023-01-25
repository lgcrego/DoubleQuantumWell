module wfunction

use constants_and_parameters
use RootFinding                     , only : Roots
use griding                         , only : sumtrap , grade

implicit none

public :: wavefunc , coefic

    contains

!===========================================
 subroutine COEFIC( enp_dwell , Nroots ) 
!===========================================
real*8  , intent(in)    :: enp_dwell(:)
integer , intent(in)    :: Nroots

!   variaveis locais
integer               :: root
real*8                :: A,B,C,D,E,F,G,H,An,Bn,Cn,Dn,En,Fn,Gn,Hn,psinormalizado
real*8                :: Am(grid_size),Bm(grid_size),Cm(grid_size),Dm(grid_size),Em(grid_size),Fm(grid_size),Gm(grid_size),Hm(grid_size)
real*8                :: k2,k1,k3,w10,w11,w12,w13,w14,w15,w16,k2l1,k1l1,k3l1,k3m,k2m,k2n,k3n,k1n,k1m
real*8                :: w1,w2,w3,w4,w5,w6,w7,w8,w9,m,n
real*8 , allocatable  :: psi_2(:)


do root = 1 , Nroots
    k1 = sqrt( v0 - enp_dwell(root))
    k2 = sqrt( enp_dwell(root) )
    k3 = sqrt( v1 - enp_dwell(root))

    m  = l1 + l0
    n  = l1 + l0 + l2

    k1l1 = k1*l1
    k1m  = k1*m
    k1n  = k1*n
    k2l1 = k2*l1
    k2m  = k2*m
    k2n  = k2*n
    k3l1 = k3*l1
    k3m  = k3*m
    k3n  = k3*n


    w1  =  k2/cos(k2l1)
    w2  = -k2*tan(k2l1)*exp(k1l1) - k1*exp(k1l1)
    w3  = -k2*tan(k2l1)*exp(-k3l1) + k3*exp(-k3l1)
    w4  = -k1*exp(-k3m) - k3*exp(-k3m)
    w5  =  k1*cos(k2m) + k2*sin(k2m)
    w6  =  k1*sin(k2m) - k2*cos(k2m)
    w7  =  k2/cos(k2n)
    w8  = -k2*exp(-k3n)*tan(k2n) + k3*exp(-k3n)
    w9  =  k2*cos(k2l1)/k1 + sin(k2l1)
    w10 =  w1*exp(k1l1)/w9 + w2
    w11 =  w1*exp(-k3l1)/w9 + w3
    w12 = -w11*exp(k1m)/w10 + exp(-k3m)
    w13 =  w4*cos(k2m)/w12 + w5
    w14 =  w4*sin(k2m)/w12 + w6
    w15 = -cos(k2n)*w14/w13 + sin(k2n)
    w16 =  w7*exp(-k3n)/w15 + w8


    H = one
    G = H*exp(-k3n) / w15
    F = -w14*G / w13
    E = (F*cos(k2m) + G*sin(k2m)) / w12
    D = -E*w11 / w10
    C = (D*exp(k1l1) + E*exp(-k3l1)) / w9
    B = k2*(C) / k1
    A = B

    if( .not. allocated(psi_2) ) allocate( psi_2(grid_size) )
    
    where( x(:,2) == 1 ) psi_2 = (A*exp( k1*x(:,1) ))**2               
    where( x(:,2) == 2 ) psi_2 = ( B*cos( k2*x(:,1) ) + C*sin( k2*x(:,1) ))**2
    where( x(:,2) == 3 ) psi_2 = ( D*exp( k1*x(:,1) ) + E*exp(-k1*x(:,1) ))**2
    where( x(:,2) == 4 ) psi_2 = ( F*cos( k2*x(:,1) ) + G*sin( k2*x(:,1) ) )**2
    where( x(:,2) == 5 ) psi_2 = ( H*exp( -k1*x(:,1)) )**2

    psinormalizado = sumtrap( 1 , grid_size , x(:,1) , psi_2 )

    psinormalizado = one /sqrt(psinormalizado)

    Hn = psinormalizado*H
    Fn = psinormalizado*F
    Gn = psinormalizado*G
    En = psinormalizado*E
    Dn = psinormalizado*D
    Cn = psinormalizado*C
    Bn = psinormalizado*B 
    An = psinormalizado*A  
 
    Hm(root) = Hn
    Fm(root) = Fn
    Gm(root) = Gn
    Em(root) = En
    Dm(root) = Dn
    Cm(root) = Cn
    Bm(root) = Bn
    Am(root) = An

 ! Norma
 call  wavefunc(enp_dwell,root,An,Bn,Cn,Dn,En,Fn,Gn,Hn)

enddo

end subroutine COEFIC

!=====================================================================
 subroutine wavefunc(enp_dwell , root,A,B,C,D,E,F,G,H )
!=====================================================================
real*8  , intent(in) :: enp_dwell(:)
integer , intent(in) :: root 
real*8  , intent(in) :: A , B , C , D , E , F , G , H

!   variveis locais
real*8  :: rk1 , rk2 , rk3 
real*8  :: norm_phi
real*8 , allocatable :: phi(:,:)
real*8 , allocatable :: phi_2(:)

if( .NOT. allocated(Phi))    allocate( Phi( size(x(:,1)) , N_of_roots ) , source = 0.d0 )
if( .NOT. allocated(Psi))    allocate( Psi( size(x(:,1)) , N_of_roots ) , source = 0.d0 )
if( .NOT. allocated(phi_2) ) allocate( phi_2(grid_size) )

rk1 = sqrt( v0 - enp_dwell(root))
rk2 = sqrt( enp_dwell(root) )
rk3 = sqrt( v1 - enp_dwell(root))

where( x(:,2) == 1 ) phi(: , root) = A*exp( rk1*x(:,1) )
where( x(:,2) == 2 ) phi(: , root) = B*cos( rk2*x(:,1) ) + C*sin( rk2*x(:,1) )
where( x(:,2) == 3 ) phi(: , root) = D*exp( rk1*x(:,1) ) + E*exp(-rk3*x(:,1) )
where( x(:,2) == 4 ) phi(: , root) = F*cos( rk2*x(:,1) ) + G*sin( rk2*x(:,1) )
where( x(:,2) == 5 ) phi(: , root) = H*exp( -rk3*x(:,1))

phi_2           = phi(:,root) * phi(:,root)
norm_phi        = sumtrap( 1 , grid_size , x(:,1) , phi_2 )
phi(:,root)     = phi(:,root) / sqrt(norm_phi)

!if( root < 5 ) then
!do i = 1 , grid_size

!    write(200 + root , 13) x(i,1) , Phi(i , root)

!enddo
!endif

phi_adiabatic( : , root ) = phi( : , root )

deallocate( phi )

13  format(7d15.4)

end subroutine wavefunc
!
!
!
end module wfunction
