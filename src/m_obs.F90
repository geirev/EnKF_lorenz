module m_obs
contains
subroutine obs(ref,d,ih,v,dt,vars,deltaobs)
   use mod_dimensions
   use m_random
   implicit none
   type(variances), intent(in) :: vars
   integer, intent(out) :: ih(nrmes)     ! Measurement matrix
   real, intent(out) :: v(nrmes,nrmes)   ! Measurement weigth matrix
   real, intent(out) :: d(3,nrmes)       ! data
   real, intent(in)  :: ref(3,ndim)      ! state variable
   real, intent(in)  :: dt
   real, intent(out) :: deltaobs

   integer, parameter :: naux=33*nrmes
   real, allocatable :: aux(:)

   integer i,ii

   allocate(aux(naux))

! Generating the measurement matrix

   deltaobs=float(ndim-1)*dt/float(nrmes)
   print *,'delta_obs=',deltaobs

   do i=1,nrmes
      ii=nint(deltaobs*float(i)/dt+1.0)
      ih(i)=ii
   enddo
   if (ih(nrmes) > ndim) then
      print *,'OBS: ih(nrmes) >ndim',ih(nrmes),ndim
   endif

   
   call random(d,3*nrmes)
   do i=1,nrmes
      d(:,i)=ref(:,ih(i))+sqrt(vars%mes)*d(:,i)
   enddo

! Generating the measurement error covariance matrix
   v = 0.0
   do i = 1,nrmes
      v(i,i)=vars%mes
   enddo

   deallocate(aux)

end subroutine obs
end module m_obs
