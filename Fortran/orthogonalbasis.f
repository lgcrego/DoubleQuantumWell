module basis
use f95_precision
use blas95
use lapack95

use constants_and_parameters
use RootFinding                , only : Roots   
use griding                    , only : grade , sumtrap
use g_state                    , only : wf_gstate , coeficientes
use wfunction                  , only : wavefunc , coefic
use wfunction_cont             , only : wavefunc_cont , coefic_cont

implicit none

public :: dinamic_eletron 

contains

!============================================================================
 Subroutine dinamic_eletron( enp ,nivel , step , tempo , freq , sig , level )
!============================================================================
real*8     , intent(in)  :: enp(:)
real*8     , intent(in)  :: tempo , freq , sig
integer    , intent(in)  :: step , nivel , level

! local variable
complex*16 , allocatable :: wave_t(:)  , temporal_operator(:)
real*8     , allocatable :: coefi(:) , mat_rho(:)
real*8                   :: norm , pop1 , pop2 , energy_wf
integer                  :: n , i , j , m

if ( .not. allocated(coefi_phi_t)) allocate( coefi_phi_t( N_all_roots , 2 ) , source = (0.d0,0.d0) )
if ( .not. allocated(coefi_x_t))   allocate( coefi_x_t( grid_size , 2 ) )
if ( .not. allocated(wave_t))      allocate( wave_t( grid_size ) )
if ( .not. allocated(coefi))       allocate( coefi( grid_size ) )

if ( .not. allocated(temporal_operator) ) allocate(temporal_operator( N_all_roots ))
if ( .not. allocated(mat_rho) )           allocate(mat_rho( N_all_roots ) , source=0.d0)

if( step == 1 .or. step == n_dt+1 ) then

        do n = 1 , N_all_roots
           coefi(:)             = phi_adiabatic(:,n)*phi_adiabatic(:,4)!psi(:,nivel)
           temporal_operator(n) = exp(-zi*enp(n)*dt)
           coefi_phi_t(n,1)     = sumtrap( 1 , grid_size , x(:,1) , coefi )
        enddo
        
        print*, coefi_phi_t(1,1)*conjg(coefi_phi_t(1,1)) , "c1"
        print*, coefi_phi_t(2,1)*conjg(coefi_phi_t(2,1)) , "c2"
        print*, coefi_phi_t(3,1)*conjg(coefi_phi_t(3,1)) , "c3"
        print*, coefi_phi_t(4,1)*conjg(coefi_phi_t(4,1)) , "c4"
        
        energy_zero_package = rymev*sum( enp(:)*coefi_phi_t(:,1)*conjg(coefi_phi_t(:,1)))
        
        !====================
        ! Funcao na base phi
        !====================
        forall( i = 1 : grid_size )  wave_t(i) = sum( coefi_phi_t(:,1)*phi_adiabatic(i,:) )
        
        !======================
        ! Calculo da População.
        !======================
        coefi(:)         = conjg(wave_t(:))*wave_t(:)
        norm             = sumtrap( 1 , grid_size , x(:,1) , coefi )
        coefi_phi_t(:,1) = coefi_phi_t(:,1)/sqrt(norm)
        
        forall( i = 1 : grid_size )  wave_t(i) = sum( coefi_phi_t(:,1)*phi_adiabatic(i,:) )
        
        do n = 1 , N_all_Roots
           mat_rho(n) = coefi_phi_t(n,1)*conjg(coefi_phi_t(n,1))
        enddo
        
        coefi_phi_t(:,2) = coefi_phi_t(:,1)
        
        call ensemble( mat_rho , tempo , enp , step , coefi_phi_t ,sig,level)
        
        !call diagram_p_x( wave_t )
        !do i = 1 , grid_size
        !        write(113,14) x(i,1) , wave_t(i)*conjg(wave_t(i)) , flaoat(step)
        !enddo
        
        open( 365 , file='pop1.dat' , position='append' )
            pop1 = sumtrap( 1 , sum(p(1:2)) + p(3)/2 , x(:,1) ,  coefi(:))               ! Psi square for first region
            write(365 , 13 ) tempo , pop1
        close(365)
        
        open( 366 , file='pop2.dat' , position='append' )
            pop2 = sumtrap( sum(p(1:2)) + p(3)/2 , grid_size , x(:,1) , coefi(:))     ! Psi square for second region
            write(366 , 13 ) tempo , pop2
        close(366)
        
        open( 367 , file='poptotal.dat' , position='append' )
            write(367 , 13 ) tempo , pop1+pop2
        close(367)
        
        
        !=============================
        ! Mudança de Base.  phi -> x
        !=============================
        forall( i = 1 : grid_size ) coefi_x_t(i,1) = sum(coefi_phi_t(:,1)*phi_adiabatic(i,:))
   
        deallocate( psi )

end if


if( step > 1 ) then
        !=============================
        ! Mudança de Base.  x -> phi
        !=============================
        forall( n = 1 : N_all_roots ) coefi_phi_t(n,2) = sum(coefi_x_t(:,1)*phi_adiabatic(:,n))

        !=============================
        ! Evolução temporal
        !=============================
        do n = 1 , n_all_roots
           temporal_operator(n) = exp(-zi*enp(n)*dt)
           coefi_phi_t(n,2)     = temporal_operator(n)*coefi_phi_t(n,2)
        enddo

        !==========================
        ! Funcao na base phi
        !==========================
        forall( i = 1 : grid_size ) wave_t(i) = sum( coefi_phi_t(:,2)*phi_adiabatic(i,:) )
  
        coefi(:)         = conjg(wave_t(:))*wave_t(:)
        norm             = sumtrap( 1 , grid_size , x(:,1) , coefi )
        coefi_phi_t(:,2) = coefi_phi_t(:,2)/sqrt(norm)
        
        forall( i = 1 : grid_size )  wave_t(i) = sum( coefi_phi_t(:,2)*phi_adiabatic(i,:) )
        
        !do i = 1 , grid_size
        !        write(113,14) x(i,1) , wave_t(i)*conjg(wave_t(i)) , float(step)
        !enddo

        do n = 1 , N_all_Roots
           mat_rho(n) = coefi_phi_t(n,2)*conjg(coefi_phi_t(n,2))
        enddo

        open( 368 , file='coefi.dat' , position='append' )
             write(368 , 13 ) tempo , mat_rho(1) , mat_rho(2) , mat_rho(3) , mat_rho(4) , mat_rho(5)
        close( 368 )

        call ensemble( mat_rho , tempo , enp , step , coefi_phi_t , sig, level )

        deallocate( mat_rho )

        !=============================
        ! Calculo da População.
        !=============================
        coefi(:) = wave_t(:)*conjg(wave_t(:))

        pop1 = sumtrap( 1 , sum(p(1:2)) + p(3)/2 , x(:,1) ,  coefi(:))            ! Psi square for first region
        open( 365 , file='pop1.dat' , position='append' )
            write(365 , 13 ) tempo , pop1
        close( 365 )

        pop2 = sumtrap( sum(p(1:2)) + p(3)/2 , grid_size , x(:,1) , coefi(:))     ! Psi square for second region
        open( 366 , file='pop2.dat' , position='append' )
            write(366 , 13 ) tempo , pop2
        close( 366 )

        open( 367 , file='poptotal.dat' , position='append' )
            write(367 , 13 ) tempo , pop1+pop2
        close( 367 )

        !=============================
        ! Mudança de Base.  phi -> x
        !=============================
        forall(i=1:grid_size)  coefi_x_t(i,2) = sum(coefi_phi_t(:,2)*phi_adiabatic(i,:))
    
        !if( sig == -1.d0 .and. l2 == 150.d0/aba ) then
        !        coefi_x_t(:,1)          = coefi_x_t(:,2)
        !print*, "deu certo?"
        !        coefi_phi_t(:,1)        = conjg(coefi_phi_t(:,2))
        !else
        coefi_x_t(:,1)          = coefi_x_t(:,2)
        coefi_phi_t(:,1)        = coefi_phi_t(:,2)
        !end if
        !    call diagram_p_x( wave_t )
end if

13 format(7d15.4)
14 format(7e15.4)

! mencoder mf://*.png -mf w=1000:h=600:fps=20:type=png -ovc copy -oac copy -o din_wave.avi

energy_wf = rymev*sum( enp(:)*coefi_phi_t(:,1)*conjg(coefi_phi_t(:,1)))

open( 327 , file='energy_wf.dat' , position='append' )
    write(327 , 13 ) tempo , energy_wf
close( 327 )

open( 325 , file='denergy.dat' , position='append' )
    write(325 , 13 ) tempo , (energy_wf - energy_zero_package)
close( 325 )

!========================================================================
! Não esquecer de descomentar a linha 113 para fazer o video
!========================================================================
! ESCREVE AS FUNÇÕES QUE SERÃO DESENHADAS NO GNUPLOT.
!========================================================================
!open(20 , file="plot.gnu" )
! write(20,fmt='(a)')        'clear'
! write(20,fmt='(a)')        'reset'
! write(20,fmt='(a)')        'unset key'
! write(20,fmt='(a)')        "set title'Dinâmica do Pacote'"
! write(20,fmt='(a)')        "set terminal pngcairo size 1000,600 enhanced font 'Verdana,20'"
! write(20,fmt='(a,f6.1,a,f6.1,a)') "set xrange[-100:350]"
! write(20,fmt='(a)')        "set yrange[0:300]"
! write(20,fmt='(a)')        "set y2range[0:300]"
! write(20,fmt='(a)')        "set xlabel'Comprimento'"
! write(20,fmt='(a)')        "set y2label'Energy'"
! write(20,fmt='(a)')        "set ylabel '20*{/Symbol y} + E_{package}' textcolor rgb 'blue'"
! write(20,fmt='(a)')        "set xtics 100"
!\/ ja estava comentado
! write(20,fmt='(a)')        'set ytics 300 textcolor rgb "blue"'
!/\ ja estava comentado
! write(20,fmt='(a)')        "set y2tics 250"
! write(20,fmt='(a)')        "set pointsize 0.005"
! write(20,fmt='(a)')        'set border linewidth 1.0'
! write(20,fmt='(a)')        "set samples 300000"
! write(20,fmt='(a,i0,a)')        "system('mkdir -p animacao')"


!\/ ja estava comentado
!    write(*,*)  '   Frame = ' , tlinha
!/\ ja estava comentado
!    write(20,fmt='(a,i0,a,i0,a)') "outfile = sprintf('animacao/geral%05.0f.png',",step,")"
!    write(20,fmt='(a)')      "set output outfile"
!    write(20,fmt='(a,i3,a,i3,a)') "plot '-' u 1:2 w l lw 4 lc 3 , '-' u 1:4 w l lw 4 lc 0"

!=================================================
!=================================================
!do i = 1 , grid_size
!    write( 20 , 14 ) x(i,1)*aba , 20*conjg(wave_t(i))*wave_t(i) + energy_wf , u(i)
!enddo
!
! write(20,fmt='(a)') "unset output outfile"
! write(20,fmt='(a)') 'e'
! write(20,fmt='(a)') ' '

!==========================================

end subroutine
!
!
!
!=================================
 Subroutine diagram_p_x( wave_t )
!=================================
complex*16 , intent(in)  ::  wave_t(:)

! local variable
complex*16 , allocatable :: wave_tlinha(:)
real*8     , allocatable :: vetor_re_p(:)
real*8     , allocatable :: vetor_img_p(:)
real*8     , allocatable :: vetor_x(:)
real*8     :: medx , medp_re , a(grid_size) , medp_img , medp_square
integer    :: i , n , k , j

forall( i = 1 : grid_size-1) a(i) = x(i+1,1) - x(i,1)
a(grid_size) = a(grid_size-1)

if ( .not. allocated(wave_tlinha))  allocate( wave_tlinha( grid_size ) )
if ( .not. allocated(vetor_re_p))   allocate( vetor_re_p( grid_size ) )
if ( .not. allocated(vetor_img_p))  allocate( vetor_img_p( grid_size ) )
if ( .not. allocated(vetor_x))      allocate( vetor_x( grid_size ) )

k = 1
do n = 1 , N_of_domains
    do i = 1 , p(n)
        wave_tlinha( k ) = -zi*(wave_t( k + 1 ) - wave_t( k ))/a(k)  
        k = k + 1
    enddo
enddo

vetor_re_p(:)  = real( conjg( wave_t(:) )*wave_tlinha(:) )
medp_re   = sumtrap( 1 , grid_size , x(:,1) , vetor_re_p )

vetor_img_p(:) = imag( conjg( wave_t(:) )*wave_tlinha(:) )
medp_img  = sumtrap( 1 , grid_size , x(:,1) , vetor_img_p )


forall( i = 1 : grid_size-1 ) vetor_x(i) = conjg(wave_t(i))*x(i,1)*aba*wave_t(i)
medx = sumtrap( 1 , grid_size , x(:,1) , vetor_x )

open( 60 , file='diagramox-real.dat' , position='append' )
write( 60 , 13 ) medx , medp_re
close( 60 )

open( 61 , file='diagrampx-img.dat' , position='append' )
write( 61 , 13 ) medx , medp_img
close( 61 )

deallocate(wave_tlinha,vetor_re_p,vetor_img_p,vetor_x)

13 format(7d15.4)

end subroutine
!
!
!
!=================================================================================
 Subroutine ensemble( mat_rho , tempo , energy , step , coefi_phi_t , sig , level)
!=================================================================================
real*8     , intent(in) :: mat_rho(:)
real*8     , intent(in) :: tempo , sig
real*8     , intent(in) :: energy(:)
integer    , intent(in) :: step , level
complex*16 , intent(in) :: coefi_phi_t(:,:)

!local variables
real*8  ,allocatable  :: entropy(:) , trace(:) , dentropy(:) , soma(:)!, trabalho(:), calor(:), e_package(:)
real*8  :: kbt  , e_package , calor , trabalho , force
integer :: info , n , m

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
    
    force_nad(1) = force*float(nivel)*dl*sig/aba
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
