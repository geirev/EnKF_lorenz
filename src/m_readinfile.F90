module m_readinfile
contains
subroutine readinfile(tfin,var,dt,mode)
   use mod_dimensions  
   implicit none
   real,            intent(out) :: tfin
   type(variances), intent(out) :: var
   real,            intent(out) :: dt
   integer,         intent(out) :: mode

! The weigths should here be specified as follows:
   open(10,file='infile.in')
      read(10,*)tfin              ! Final time 
      read(10,*)var%dyn        ! Model error variance
      read(10,*)var%ini        ! Initial error variance
      read(10,*)var%mes        ! Measurement error variance
      read(10,*)mode           ! Operation: 0 - pure ensemble int
                               !            1 - ensemble Kalman filter
                               !            2 - ensemble Kalman smoother
                               !            3 - Lagged ensemble Kalman smoother
   close(10)

   dt=(tfin-0.0)/float(ndim-1)

   print '(a,f10.6)','Time step dt=',dt
   print '(a,f10.2)','Final time T=',tfin
   print '(a,g13.4)','var%dyn=',var%dyn
   print '(a,g13.4)','var%ini=',var%ini
   print '(a,g13.4)','var%mes=',var%mes
   
end subroutine
end module
