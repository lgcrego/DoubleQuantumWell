module Joule

use f95_precision
use blas95
use lapack95
use constants_and_parameters

implicit none

public :: ensemble 

contains

!=================================================================================
 Subroutine ensemble( mat_rho , tempo , energy , step , coefi_phi_t , sinal , level)
!=================================================================================
real*8     , intent(in) :: mat_rho(:)
real*8     , intent(in) :: tempo , sinal
real*8     , intent(in) :: energy(:)
integer    , intent(in) :: step , level
complex*16 , intent(in) :: coefi_phi_t(:,:)

!local variables
real*8  ,allocatable  :: entropy(:) , trace(:) , dentropy(:) , soma(:)!, trabalho(:), calor(:), e_package(:)
real*8  :: e_package , calor , trabalho , force
integer :: n , m

if(.not. allocated(ptaumn) )       allocate( ptaumn( 4,4 ) )
if(.not. allocated(entropy) )      allocate( entropy( 10 ) )
if(.not. allocated(dentropy) )     allocate( dentropy( 10 ) )
if(.not. allocated(trace) )        allocate( trace( n_all_roots ) )
if(.not. allocated(ener_zero) )    allocate( ener_zero( n_all_roots ) )
if(.not. allocated(rho_zero) )     allocate( rho_zero( n_all_roots ) )
if(.not. allocated(rho_inic) )     allocate( rho_inic( 10 ) )
if(.not. allocated(entropy_zero) ) allocate( entropy_zero( 10 ) )
if(.not. allocated(c_phi_zero) )   allocate( c_phi_zero( n_all_roots ) )
if(.not. allocated(soma) )         allocate( soma( n_all_roots ) )
if ( .not. allocated(force_nad))   allocate( force_nad( 2 ) )

if(.not. allocated(sum_entropy) )  allocate( sum_entropy( 2 ) , source = 0.d0 )
!if(.not. allocated(calor) )        allocate( calor( n_all_roots ) , source = 0.d0 )
!if(.not. allocated(trabalho) )     allocate( trabalho( n_all_roots ) , source = 0.d0 )
!if(.not. allocated(e_package) )    allocate( e_package( n_all_roots ) )

!============================================
!       VERIFICA O TRACO DE PHO
!============================================
open( 368 , file='coefi.dat' , position='append' )
     write(368 , 13 ) tempo , mat_rho(1) , mat_rho(2) , mat_rho(3) , mat_rho(4) , mat_rho(5)
close( 368 )

forall( n = 1 : N_all_roots ) trace(n) = mat_rho( n )

open( 76 , file='traco.dat' , position='append' )
    write( 76 , 13 ) tempo , sum(trace(:))
close( 76 ) 

deallocate( trace )

!============================================
!       CALCULA O TRACO DE PHO H
!============================================
if( step == 1 ) then
trabalhototal   = 0.d0
calortotal      = 0.d0
sum_entropy(1)  = 0.d0
force_nad(1)    = 0.d0
force_nad(2)    = 0.d0
        
rho_zero(:)  = mat_rho(:)
ener_zero(:) = energy(:)
do n = 1 , 10
   entropy_zero(n) = -mat_rho(n)*log(mat_rho(n))               
enddo

write( 150 , 13 ) sum(entropy_zero(:))

beta = 0.35d0

print*, energy(1) , energy(2) , energy(3) , energy(4) 
rho_inic(1) = exp(-beta*energy(1))/sum(exp(-beta*energy(1:4)))!0.993d0
rho_inic(2) = exp(-beta*energy(2))/sum(exp(-beta*energy(1:4)))!0.6916d-3
rho_inic(3) = exp(-beta*energy(3))/sum(exp(-beta*energy(1:4)))!0.1261d-5    
rho_inic(4) = exp(-beta*energy(4))/sum(exp(-beta*energy(1:4)))!0.2342d-5

print*, rho_inic(:)
stop

end if

if( step > 1 ) then

    !=============================
    ! Calculo da forca_na
    !=============================
    force = 0.d0
    do n = 1 , n_all_roots-1
       do m = n + 1 , n_all_roots
          force = force - 2.d0*250.d0*(real(conjg(coefi_phi_t(m,2))*coefi_phi_t(n,2)))*phi_adiabatic(parameter_n,m)*phi_adiabatic(parameter_n,n)!*(energy(n) - energy(m))/(ener_zero(n) - ener_zero(m))
       enddo
    enddo
    
    force_nad(1) = force*float(nivel)*dl*sinal/aba
    force_nad(2) = force_nad(2) + force_nad(1)
    
    open( 174 , file='force_nad_ins.dat' , position='append' )
    write(174 , 13 ) tempo , force
    close( 174 )
    
    open( 175 , file='work_fric.dat' , position='append' )
        write(175 , 13 ) tempo , force_nad(2)
    close( 175 )
    
    !=====================
    !       trabalhos
    !=====================
    trabalho  = rymev*sum( rho_Zero(:)*( energy(:) - ener_zero(:) ) )
    e_package = rymev*sum( energy( : )*mat_rho(:)  ) - rymev*sum( ener_zero(: )*rho_zero(:) )
    calor     = rymev*sum( ener_zero(:)*(mat_rho(:) - rho_zero(:)))
    
    open( 80 , file='trabalho.dat', position='append'  )
        write( 80 , 13 ) tempo , trabalho
    close( 80 ) 
    
    open( 81 , file='calor.dat' , position='append' )
        write( 81 , 13 ) tempo , calor
    close( 81 )
    
    open( 82 , file='termoextra.dat' , position='append' )
        write( 82 , 13 ) tempo , e_package - calor - trabalho
    close( 82 )
    
    open( 83 , file='de_package.dat' , position='append' )
        write( 83 , 13 ) tempo , e_package
    close( 83 ) 
    
    forall( n = 1 : 4 ) ptaumn(level,n) = mat_rho(n) 
    
    !==========================================================
    !       Jarzynski Equality
    !==========================================================
    if(step == n_dt ) then
    open( 120 , file='diagram.dat', position='append'  )
        write(120,13) energy(1) - energy(4) , ptaumn(1,4)*rho_inic(4)
        write(120,13) energy(2) - energy(4) , ptaumn(2,4)*rho_inic(4)
        write(120,13) energy(3) - energy(4) , ptaumn(3,4)*rho_inic(4)
    
        write(120,13) energy(4) - energy(3) , ptaumn(4,3)*rho_inic(3)
        write(120,13) energy(2) - energy(3) , ptaumn(2,3)*rho_inic(3)
        write(120,13) energy(1) - energy(3) , ptaumn(1,3)*rho_inic(3)
    
        write(120,13) energy(4) - energy(2) , ptaumn(4,2)*rho_inic(2)
        write(120,13) energy(3) - energy(2) , ptaumn(3,2)*rho_inic(2)
        write(120,13) energy(1) - energy(2) , ptaumn(1,2)*rho_inic(2)
    
        write(120,13) energy(4) - energy(1) , ptaumn(4,1)*rho_inic(1)
        write(120,13) energy(3) - energy(1) , ptaumn(3,1)*rho_inic(1)
        write(120,13) energy(2) - energy(1) , ptaumn(2,1)*rho_inic(1)
        write(120,13) 0.d0 , rho_inic(1)*ptaumn(1,1) + rho_inic(2)*ptaumn(2,2)+rho_inic(3)*ptaumn(3,3)+ rho_inic(4)*ptaumn(4,4)
    close( 120 )
    end if
    
    !==========================================================
    !       DIAGONALIZA RHO PARA OBTER OS AUTO ESTADOS
    !==========================================================
    !call heev( mat_rho , eigen_value ,'N' , 'U' , info )
    
    do n = 1 , 10
        entropy(n)      = -mat_rho(n)*log(mat_rho(n))               
        entropy_zero(n) = -rho_zero(n)*log(mat_rho(n))               
    enddo
    sum_entropy(1)  = sum_entropy(1) + sum( entropy(:) - entropy_zero(:) )
    
    open( 75 , file='deentropia.dat', position='append'  )
        write( 75 , 13 ) tempo , calor/sum(entropy(:) - entropy_zero(:) )
    close( 75 ) 
    
    open( 73 , file='kbt.dat', position='append'  )
        write( 73 , 13 ) tempo , sum_entropy(1)
    close( 73 ) 
    
    trabalhototal   = trabalho + trabalhototal
    calortotal      = calor + calortotal
    
    open( 90 , file='workall.dat', position='append'  )
        write( 90 , 13 ) tempo , trabalhototal
    close( 90 ) 
    
    open( 91 , file='heatall.dat' , position='append' )
        write( 91 , 13 ) tempo , calortotal
    close( 91 )
    
    open( 92 , file='allenerg.dat' , position='append' )
        write( 92 , 13 ) tempo , calortotal + trabalhototal + energy_zero_package - rymev*sum( energy( : )*mat_rho( : ) )
    close( 92 )

endif

rho_zero(:)  = mat_rho(:)
ener_zero(:) = energy(:)

13 format(7d15.4)

end subroutine

end module
