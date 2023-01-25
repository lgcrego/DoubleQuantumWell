module RootFinding

use constants_and_parameters
use omp_lib

implicit none

public :: Roots , func_dwell , func_swell , zbrent , func_cont_dwell

private

! module variables ...
 character(len=10) :: my_name

contains

!======================
 function Roots( func ) 
!======================

! define interface for user provided function func(e) ...
interface
    real*8 function func(e)
       real*8 , optional , intent(in) :: e
    end function func
end interface    

! local variables
real*8  :: f1 , f2 , x1 ,x2 , step , dummy
integer ::  i , j , l , i1 , i2

! local parameters
integer :: nsteps = 500000 , extra_steps = 200000

real*8 , allocatable :: temp(:) , Roots(:) 

!--------------------------------------------------------------------
! identify func() ...
dummy = func()

select case (my_name)
       
       case( "func_swell" , "func_dwell" )
           i1 = 1
           i2 = nsteps -1

       case("cont_dwell")
           i1 = nsteps
           i2 = nsteps + extra_steps -1

end select
!--------------------------------------------------------------------
! calculo das raízes de func

if(.not. allocated(temp) ) allocate( temp(grid_size) , source = 0.d0 )

step = v0 / float(nsteps)

j = 1
l = 1
do i = i1 , i2

    x1 = float(i)*step       + tiny_epsilon
    x2 = (float(i) + 1)*step - tiny_epsilon

    f1 = func(x1)
    f2 = func(x2)

    if( f1*f2 <= zero ) then  ! found j_th root ...
    !----------------------------------------
       if( x1 < v0 ) then
          temp(j) = zbrent(x1,x2,func)
          j = j + 1
       else
          if( l <= N_all_roots - N_of_roots ) then
             temp(l) = zbrent(x1,x2,func)
             l = l + 1
          endif
       endif
    !----------------------------------------
    endif

end do

if( l == 1 ) then  ! <== no states in the continuum
        N_of_roots = j - 1   
        if(.not. allocated( Roots ) )allocate( Roots(N_of_Roots) , source = temp(1:N_of_Roots) )
else               ! <== there are states in the continuum
        N_of_roots_cont = l - 1
        if(.not. allocated( Roots ) )allocate( Roots(N_of_Roots_cont) , source = temp(1:N_of_Roots_cont) )
end if

deallocate( temp )

if (verbose) write(*,*) "number of roots = ", N_of_roots , N_of_roots_cont , N_all_roots
!----------------------------------------------------------

28 format(5d15.8)

end function Roots
!
!
!=================================
 function ZBRENT( x1 , x2 , func )
!=================================
real*8, intent(in) :: x1 , x2

! define interface for user provided function func(e) ...
interface
    real*8 function func(e)
        real*8 , optional , intent(in) :: e
    end function func
end interface  

!   variaveis locais
real*8   :: a,b,fa,fb,fc,c,d,e,tol1,xm,s,p,q,r,zbrent
integer  :: iter

! local parameters
real*8,  parameter :: tol = 1.d-18
real*8,  parameter :: epsi=1.d-10
integer, parameter :: imax=1000

a  = x1
b  = x2
fa = func(a)
fb = func(b)

if(fb*fa > zero) Stop 'root must be bracketed for ZBRENT'
     fc = fb
     do iter=1 , imax
          if(fb*fc > zero) then
               c  = a
               fc = fa
               d  = b-a
               e  = d
          endif

          if(abs(fc) < abs(fb)) then
               a  = b
               b  = c
               c  = a
               fa = fb
               fb = fc
               fc = fa
          endif

     tol1 = two*epsi*abs(b)+half*tol
     xm = half*(c-b)

          if(abs(xm) <= tol1 .or. fb == zero) then
               zbrent = b
          return
          endif

          if(abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
               s = fb/fa
               if(a == c) then
                    p = two*xm*s
                    q = one - s
               else
                    q = fa/fc
                    r = fb/fc
                    p = s*(two*xm*q*(q-r)-(b-a)*(r-one))
                    q = (q-one)*(r-one)*(s-one)
               endif

               if(p > zero) q = -q
                    p = abs(p)
                    if(two*p < min(3*xm*q - abs(tol1*q),abs(e*q))) then
                         e = d
                         d = p/q
                    else
                         d = xm
                         e = d
                    endif
               else
                    d = xm
                    e = d
               endif
               a  = b
               fa = fb

               if(abs(d) > tol1) then
                    b = b + d
               else
                    b = b + sign(tol1,xm)
               endif
     fb = func(b)
enddo

stop 'zbrent exceeding maximum iterations'

zbrent = b

return

end function ZBRENT
! 
!======================
 function func_swell(e)
!======================
real*8 , optional , intent(in) :: e

! variaveis locais
real*8 :: k1 , k2 , k1l1 , k2l1 , c1 , c2
real*8 :: func_swell

if( .not. present(e) ) then
   my_name = "func_swell"
   func_swell = one
   return
   end if

 k1 = sqrt(v0 - e)
 k2 = sqrt(e)

 k1l1 = k1*l1
 k2l1 = k2*l1

 c1 = ( k2*cos(k2l1)/k1 + sin(k2l1) )
 c2 = ( k2*( -k2*sin(k2l1)/k1 + cos(k2l1) ))

 func_swell = k1*( c2*exp(-k1l1) + c1*k1*exp(-k1l1) )

end function func_swell
!
!
!======================
 function func_dwell(e)
!======================
real*8 , optional , intent(in) :: e

!   variaveis locais
real*8 :: m,n,w1,w2,w3,w4,w5,w6,w7,w8,w9
real*8 :: w10,w11,w12,w13,w14,w15,w16
real*8 :: k1,k2,k3,k1m,k2l1,k1l1,k3l1,k3m,k2m,k2n,k3n,k1n
real*8 :: func_dwell

if( .not. present(e) ) then
   my_name = "func_dwell"
   func_dwell = one
   return
   end if
 
 k1=sqrt(v0-e)       
 k2=sqrt(e)
 k3=sqrt(v1-e)
 m = l1 + l0
 n = l1 + l0 + l2

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

 func_dwell = k1*w9*w10*w12*w13*w15*w16

end function func_dwell
!
!===========================
 function func_cont_dwell(e)
!===========================
real*8 , optional , intent(in) :: e

!   variaveis locais
real*8 :: m,n,k1,k2,k2l1,k1l1,k2m,k2n,k1n,k1m, func_cont_dwell , k1left, k1right
real*8 :: w1 ,w2 , w3 , w4 , w5 , w6 , w7 , w8 , w9 , w10 ,w11 ,w12 ,w13

if( .not. present(e) ) then
   my_name = "cont_dwell"
   func_cont_dwell = one
   return
   end if

 k1=sqrt(e-v0)
 k2=sqrt(e)
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

 func_cont_dwell = cos(k1left)*k1*w1*w3*w5*w7*w9*w11*w13

end function func_cont_dwell

end module RootFinding
