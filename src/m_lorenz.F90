module m_lorenz
contains
subroutine lorenz(na,nb,dt,sol,modvar)
! Geir Evensen 22/3-90
   use mod_dimensions
   use m_random
   use m_fex
   implicit none

   real, intent(in)     :: dt          ! Time step
   integer, intent(in)  :: na          ! Start time
   integer, intent(in)  :: nb          ! end time
   real, intent(in)     :: modvar      ! model error variance
   real, intent(inout)  :: sol(neq,ndim) ! Solution vector

   integer, parameter :: lrw=100
   integer, parameter :: liw=50

   real tt(neq)
   real y(neq)
   real atol(neq)
   real rwork(lrw)
   integer iwork(liw)

   integer itask,istate,iopt,mf,iout,itol,k,jex
   real t1,t2,rtol


   mf=10            ! method flag
   itol=2           ! itol=2 mean atol array
   rtol=1.e-5       ! Relative tolerance parameter
   atol(:)=1.e-5    ! Absolute tolerance parameter
   itask=1          ! normal computation with output in tout
   istate=1         ! integer flag
   iopt=0           ! no optional inputs

   y(:)=sol(:,na)
   if (nb > ndim) stop 'nb>ndim'
   do iout=na+1,nb
      t1=real(iout-2)*dt   
      t2=real(iout-1)*dt   
      iwork(6)=2000
      iwork(7)=1

      do k=1,100
         call slsode(fex,neq,y,t1,t2,itol,rtol,atol,itask,&
                    istate,iopt,rwork,lrw,iwork,liw,jex,mf)
         if (istate == -1) then
            istate=1
            write(*,*)'istate=-1   t=',t2
         else
            exit
         endif
      enddo

      if (istate <= 0) then
         write(*,*)'istate=',istate
         stop
      endif

      if (modvar /= 0.0 ) then
         call random(tt,3)
         sol(1,iout)=y(1)+1.0*sqrt(0.1491)*tt(1)
         sol(2,iout)=y(2)+1.0*sqrt(0.9048)*tt(2)
         sol(3,iout)=y(3)+1.0*sqrt(0.9180)*tt(3)
      else
         sol(1,iout)=y(1)
         sol(2,iout)=y(2)
         sol(3,iout)=y(3)
      endif

   enddo
end subroutine lorenz
end module m_lorenz
