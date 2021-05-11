module m_readinfile
   use mod_dimensions
   integer mode              ! 0-Ensemble, 1-EnKF, 2-ES, 3-EnKS
   integer lag               ! lag used in EnKS (0 gives EnKF)
   character(len=8) fname    ! Name of output file
   real obsdt                ! time between measurements
   real dt                   ! timestep for model outputs
contains
subroutine readinfile(var,dt)
   use mod_dimensions  
   implicit none
   type(variances), intent(out) :: var
   real,            intent(out) :: dt
   real tfin

! The weigths should here be specified as follows:
   open(10,file='infile.in')
      read(10,*)nrsamp  ; print *,'nrsamp  :',nrsamp      ! Ensemble size
      read(10,*)nrmes   ; print *,'nrmes   :',nrmes       ! Number of measurement times
      read(10,*)ndim    ; print *,'ndim    :',ndim        ! Number of grid points in time
      read(10,*)tfin    ; print *,'tfin    :',tfin        ! Final time 
      read(10,*)mode    ; print *,'mode    :',mode        ! 0-Ensemble, 1-EnKF, 2-ES, 3-EnKS 
      read(10,*)lag     ; print *,'EnKS lag:',lag         ! lag for EnKS
      read(10,*)var%dyn ; print *,'var_dyn :',var%dyn     ! Model error variance
      read(10,*)var%ini ; print *,'var_ini :',var%ini     ! Initial error variance
      read(10,*)var%mes ; print *,'var_mes :',var%mes     ! Measurement error variance
   close(10)

   nrtobs=nrmes*nrobs          ! total number of mes 

   dt=(tfin-0.0)/float(ndim-1) ! timestep
   print *,'dt=',dt

   obsdt=float(ndim-1)*dt/float(nrmes)
   print *,'obsdt=',obsdt

   select case (mode)
   case(0) 
      fname='Ensemble'
   case(1) 
      fname='EnKF'
      lag=0
   case(2) 
      fname='ES'
   case(3) 
      fname='EnKS'
      if (lag == 0) stop 'You must have lag > 0 for EnKS'
   end select

   print '(a,f10.6)','Time step dt=',dt
   print '(a,f10.2)','Final time T=',tfin
   print '(a,g13.4)','var%dyn=',var%dyn
   print '(a,g13.4)','var%ini=',var%ini
   print '(a,g13.4)','var%mes=',var%mes
   
end subroutine
end module
