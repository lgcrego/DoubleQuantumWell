module constants_and_parameters
!	Declarar Unidades

real*8 , parameter  :: half = 0.5d0, zero = 0.0d0, one = 1.0d0, two = 2.0d0, five = 5.0d0, four = 4.0d0, three = 3.0d0, eight = 8.0d0 , six = 6.0d0

real*8 , parameter :: EM0 = 9.1091d-28                     !massa do eletron kg
real*8 , parameter :: EMB = 0.067*EM0
real*8 , parameter :: EP0 = 12.9
real*8 , parameter :: PI  = 3.14159265358979
real*8 , parameter :: EC  = 4.80298d-10                    !Carga do eletron em esu
real*8 , parameter :: HB  = 1.0545d-27                     !Hbar em ergs*segundo
real*8 , parameter :: AB  = HB**two*EP0/(EMB*EC**two)      !1.01d-6  cm
real*8 , parameter :: RY  = EC**two/(two*EP0*AB)           !8.77d-15 erg

real*8 , parameter :: ns = 1d-9
real*8 , parameter :: FS =  1d-15                          !fentossegundo
real*8 , parameter :: as =  1d-18                          !atossegundo
real*8 , parameter :: zs =  1d-21                          !zeptossegundo
real*8 , parameter :: HB_ev = 6.582d-18

real*8 , parameter :: ONEA   = 1.0d-8                      !Um angstron equivale a 1.0e-8 cm
real*8 , parameter :: ONEMEV = 1.6021d-15                  !Um meV equivale a 1.602e-15 joules
real*8 , parameter :: ONEFS  = 1d-15
real*8 , parameter :: RYMEV  = RY/ONEMEV                   !transforma meV em erg
real*8 , parameter :: ABA    = AB/ONEA
real*8 , parameter :: hbar   = 1.97329d7
real*8 , parameter :: mass   = 0.511003d6
real*8 , parameter :: radius = (mass)/(two*(hbar**two))
real*8 , parameter :: Cnorm  = 139.104041547323 

real*8  :: v0 , v1 , l0 , l1 , l2 , dt , dl , work0 , energy_zero_package , calortotal , trabalhototal , freqress , beta
integer :: nivel , grid_size , N_of_roots, N_of_roots_cont , N_all_roots , passo , ciclo , parameter_n , n_dt

real*8 , parameter :: lim_left = 5000.d0/aba
real*8 , parameter :: lim_right = 5000.d0/aba

integer , dimension(5)  :: p
real*8  , allocatable   :: x(:,:)
real*8  , allocatable   :: Psi(:,:)
real*8  , allocatable   :: ptaumn(:,:)
real*8  , allocatable   :: Phi_adiabatic(:,:)

complex*16  , allocatable   :: coefi_x_t(:,:) , c_phi_zero(:) 
complex*16  , allocatable   :: coefi_phi_t(:,:)
real*8      , allocatable   :: energ_phi_t(:,:),ener_zero(:) , rho_zero(:) ,entropy_zero(:) , rho_inic(:), rho_inic2(:)
real*8      , allocatable   :: work_vet(:) ,arg_c_phi_zero(:)
real*8      , allocatable   :: calor_vet(:)
real*8      , allocatable   :: force_nad(:)
real*8      , allocatable   :: sum_entropy(:)

integer    , parameter :: N_of_domains = 5

complex*16 , parameter :: zi = (0.0d0,1.d0)

integer, dimension(8) :: values

end module constants_and_parameters

