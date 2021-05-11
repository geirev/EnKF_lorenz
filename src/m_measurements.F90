module m_measurements
contains
subroutine measurements(ref,obs,ih,dt,vars,obsdt)
   use mod_dimensions
   use m_random
   implicit none
   type(variances), intent(in) :: vars
   integer, intent(out) :: ih(0:nrmes)       ! Measurement matrix
   real, intent(out) :: obs(nrobs,nrmes)   ! data
   real, intent(in)  :: ref(neq,ndim)      ! state variable
   real, intent(in)  :: dt
   real, intent(in)  :: obsdt
   integer m

   ih(0)=1 ! used as start time but not measurement there
! Setting measurement pointer to time indexes
   do m=1,nrmes
      ih(m)=nint(obsdt*real(m)/dt+1.0)
   enddo
   call random(obs,nrobs*nrmes)
   do m=1,nrmes
      obs(:,m)=ref(:,ih(m))+sqrt(vars%mes)*obs(:,m)
   enddo

end subroutine
end module
