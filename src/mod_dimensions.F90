module mod_dimensions
   integer, parameter :: nrsamp=2000     ! Ensemble size
   integer, parameter :: ndim=801        ! Dimension of state
   integer, parameter :: nrmes=80        !80 160    ! number of measurements
   integer, parameter :: nrobs=3         ! number of mes at each time
   integer, parameter :: nrtobs=nrmes*nrobs   ! total number of mes 
   real, parameter   ::  sigma=10.0          ! Prantl number
   real, parameter   ::  r=28.0              ! Normalized Rayleigh number
   real, parameter   ::  b=8.0/3.0           ! Nondimensional wave number

   type variances
      real dyn
      real mes
      real ini
   end type variances

end module mod_dimensions
