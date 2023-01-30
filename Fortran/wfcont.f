module wfunction_cont

use constants_and_parameters
use RootFinding                     , only : Roots
use griding                         , only : sumtrap , grade

implicit none

public :: wavefunc_cont , coefic_cont

    contains

!==========================================================
 subroutine COEFIC_cont( e_n_cont_dwell , Nroots_cont )
!==========================================================
real*8  , intent(in) :: e_n_cont_dwell(:)
integer , intent(in) :: Nroots_cont 

!   variaveis locais
integer               :: root
real*8                :: Am(grid_size),Bm(grid_size),Cm(grid_size),Dm(grid_size),Em(grid_size),Fm(grid_size),Gm(grid_size),Hm(grid_size),Im(grid_size),Jm(grid_size)
real*8                :: k2,k1,k2l1,k1l1,k2m,k2n,k1n,k1m, m,n , psicnormalizado , k1right,k1left
real*8                :: w1 ,w2 , w3 , w4 , w5 , w6 , w7 , w8 , w9 , w10 ,w11 ,w12 ,w13
real*8                :: A , B, C, D, E, F, G, H, I, J, An, Bn, Cn, Dn, En, Fn, Gn, Ino, Jn,Hn
real*8 , allocatable  :: psic_2(:)

do root = 1 , Nroots_cont
      k1 = sqrt( e_n_cont_dwell(root) - v0)
      k2 = sqrt( e_n_cont_dwell(root) )
  
      m = l1 + l0
      n = l1 + l0 + l2
      k1l1 = k1*l1
      k1m  = k1*m
      k1n  = k1*n
      k2l1 = k2*l1
      k2m  = k2*m
      k2n  = k2*n
      k1left     = k1*(lim_left+l0+l1)
      k1right    = k1*(lim_right+l0+l1)

      w1 = cos(k2l1)*tan(k1left)*k2/k1 + sin(k2l1)
      w2 = -k2*sin(k2l1)*tan(k1left)*k2/k1 + k2*cos(k2l1)
      w3 = w2*cos(k1l1)/w1 + k1*sin(k1l1)
      w4 = w2*sin(k1l1)/w1 - k1*cos(k1l1)
      w5 = -cos(k1m)*w4/w3 + sin(k1m)
      w6 = k1*sin(k1m)*w4/w3 + k1*cos(k1m)
      w7 = w6*cos(k2m)/w5 + k2*sin(k2m)
      w8 = w6*sin(k2m)/w5 - k2*cos(k2m)
      w9 = -cos(k2n)*w8/w7 + sin(k2n)
      w10 = k2*sin(k2n)*w8/w7 + k2*cos(k2n)
      w11 = w10*cos(k1n)/w9 + k1*sin(k1n)
      w12 = w10*sin(k1n)/w9 - k1*cos(k1n)
      w13 = -w12*cos(k1right)/w11 + sin(k1right)

      J = one
      I = -J*w12/w11
      H = (I*cos(k1n) + J*sin(k1n))/w9
      G = -H*w8/w7
      F = (G*cos(k2m) + H*sin(k2m))/w5
      E = -F*w4/w3
      D = (E*cos(k1l1) + F*sin(k1l1))/w1
      C = D*k2*tan(k1left)/k1
      B = C/tan(k1left)
      A = B*tan(k1left)

      if( .not. allocated(psic_2) ) allocate( psic_2(grid_size) )
             
      where( x(:,2) == 1 ) psic_2 = (A*cos( k1*x(:,1) ) + B*sin( k1*x(:,1) ))**2
      where( x(:,2) == 1 ) psic_2 = (C*cos( k2*x(:,1) ) + D*sin( k2*x(:,1) ))**2
      where( x(:,2) == 1 ) psic_2 = (E*cos( k1*x(:,1) ) + F*sin( k1*x(:,1) ))**2
      where( x(:,2) == 1 ) psic_2 = (G*cos( k2*x(:,1) ) + H*sin( k2*x(:,1) ))**2
      where( x(:,2) == 1 ) psic_2 = (I*cos( k1*x(:,1) ) + J*sin( k1*x(:,1) ))**2  
      
      psicnormalizado = sumtrap( 1 , grid_size , x(:,1) , psic_2 )

      psicnormalizado = one /sqrt(psicnormalizado)

      Ino = psicnormalizado*I
      Jn = psicnormalizado*J
      Hn = psicnormalizado*H
      Fn = psicnormalizado*F
      Gn = psicnormalizado*G
      En = psicnormalizado*E
      Dn = psicnormalizado*D
      Cn = psicnormalizado*C
      Bn = psicnormalizado*B 
      An = psicnormalizado*A  
 
      Jm(root) = Jn
      Im(root) = Ino
      Hm(root) = Hn
      Fm(root) = Fn
      Gm(root) = Gn
      Em(root) = En
      Dm(root) = Dn
      Cm(root) = Cn
      Bm(root) = Bn
      Am(root) = An

     ! Norma

     call  wavefunc_cont(e_n_cont_dwell , root,An,Bn,Cn,Dn,En,Fn,Gn,Hn,Ino, Jn)
enddo
    
end subroutine COEFIC_cont

!===================================================================
 subroutine wavefunc_cont(e_n_cont_dwell , root,A,B,C,D,E,F,G,H,I,J)
!===================================================================
real*8      , intent(in)  :: e_n_cont_dwell(:)
integer     , intent(in)  :: root
real*8      , intent(in)  :: A , B, C, D, E, F, G, H, I, J

!   variveis locais
real*8  :: rk1 , rk2 
real*8  :: norm_phi_cont
real*8  , allocatable :: phi_2_cont(:) , phi_cont(:,:)

if( .not. allocated(Phi_cont))    allocate( Phi_cont( size(x(:,1)) , N_of_roots_cont ) , source = 0.d0)
if( .not. allocated(phi_2_cont) ) allocate( phi_2_cont(grid_size) )

rk1 = sqrt( e_n_cont_dwell(root) - v0)
rk2 = sqrt( e_n_cont_dwell(root) )

where( x(:,2) == 1 ) phi_cont(: , root) = A*cos( rk1*x(:,1) ) + B*sin( rk1*x(:,1) )
where( x(:,2) == 2 ) phi_cont(: , root) = C*cos( rk2*x(:,1) ) + D*sin( rk2*x(:,1) )
where( x(:,2) == 3 ) phi_cont(: , root) = E*cos( rk1*x(:,1) ) + F*sin( rk1*x(:,1) ) 
where( x(:,2) == 4 ) phi_cont(: , root) = G*cos( rk2*x(:,1) ) + H*sin( rk2*x(:,1) )
where( x(:,2) == 5 ) phi_cont(: , root) = I*cos( rk1*x(:,1) ) + J*sin( rk1*x(:,1) ) 

phi_2_cont       = phi_cont(:,root)*phi_cont(:,root)
norm_phi_cont    = sumtrap( 1 , grid_size , x(:,1) , phi_2_cont )
phi_cont(:,root) = phi_cont(:,root) / sqrt(norm_phi_cont)

!do k = 1 , grid_size
!    write(240 + root , 13) x(k,1)*aba , Phi_cont(k , root)
!enddo

phi_adiabatic( : , N_of_roots + root ) = phi_cont( : , root )

deallocate( phi_cont )

13  format(7d15.4)

end subroutine wavefunc_cont
!
!
!
end module wfunction_cont
