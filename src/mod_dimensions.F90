module mod_dimensions
   integer, parameter :: neq=3      ! Number of equations
   integer, parameter :: nrobs=3    ! Number of mes at each time
   integer nrsamp                   ! Ensemble size
   integer nrmes                    ! Number of measurement times
   integer ndim                     ! Number of grid points in time
   integer nrtobs

   type variances
      real dyn
      real mes
      real ini
   end type variances

end module mod_dimensions
