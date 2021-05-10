program main 
   use mod_dimensions       ! Parameter and type  declarations
   use m_ensmean
   use m_ensvar
   use m_random
   use m_lorenz
   use m_obs
   use m_readinfile

   implicit none
   type(variances) vars
   character(len=8) fname

   integer :: lag=0
   real A(3,ndim,nrsamp)
   real ave(3,ndim)          ! predicted state variable
   real var(3,ndim)          ! predicted state variable
   real ref(3,ndim)          ! Reference solution
   real fg(3)                ! first guess initial condition
   real d(3,nrmes)           ! vector of measurements
  
   integer ih(nrmes)      ! Measurement pointers
   real v(nrmes,nrmes)    ! Measurement error covariance matrix
   real dt,tfin
   integer na,nb,nc
   real deltaobs
   integer i,j,m
   integer mode

   real Renkf(nrobs,nrobs),Denkf(nrobs,nrsamp),Eenkf(nrobs,nrsamp),Senkf(nrobs,nrsamp),aveenkf(3),innovenkf(3)
   real Res(nrtobs,nrtobs),Des(nrtobs,nrsamp),Ees(nrtobs,nrsamp),Ses(nrtobs,nrsamp),avees(nrtobs),innoves(nrtobs)

! Reading initial parameters
   call readinfile(tfin,vars,dt,mode)
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
      lag=75
   end select

! Simulating reference case
   ref(1,1)= 1.508870
   ref(2,1)=-1.531271
   ref(3,1)= 25.46091
   na=1
   nb=ndim
   call lorenz(na,nb,dt,ref,0.0)

! Defining initial cond for central forecast
   call random(fg,3)
   fg(:)=sqrt(vars%ini)*fg(:)+ref(:,1)


! Generating initial ensemble
   do j=1,nrsamp
      call random(A(1,1,j),3)
      A(:,1,j)=fg(:)+sqrt(vars%ini)*A(:,1,j)
   enddo


! Generating measurements
   call obs(ref,d,ih,v,dt,vars,deltaobs) 


! Time stepping between measurement times
   do m=1,nrmes

! Ensemble integration
      na=nint(deltaobs*float(m-1)/dt+1.0)
      nb=nint(deltaobs*float(m)/dt+1.0)
      print '(2(a,i4))','integrating over indexes: na=',na,' nb=',nb
      do j=1,nrsamp
         call lorenz(na,nb,dt,A(1,1,j),vars%dyn) 
      enddo

! EnKF (lag=0) or EnKS update with given lag
      if ((mode == 1).or.(mode == 3)) then 
         call random(Eenkf,3*nrsamp)
         Eenkf=sqrt(vars%mes)*Eenkf
         do j=1,nrsamp
            Senkf(:,j)=A(:,nb,j)
            Denkf(:,j)=d(:,m)+Eenkf(:,j)-Senkf(:,j)
         enddo
         call ensmean(Senkf,aveenkf,3,nrsamp)
         do j=1,nrsamp
            Senkf(:,j)=Senkf(:,j)-aveenkf(:)
         enddo

         call ensmean(Denkf,innovenkf,3,nrsamp)
         Renkf=0.0

         na=max(1,nb-lag)
         nc=nb-na+1
         call analysis(A(:,na:nb,:),Renkf,Eenkf,Senkf,Denkf,innovenkf,3*nc,nrsamp,3,.true.,0.99,13,.false.,.false.,.true.,0,1.0,1)
      endif

   enddo

   if (mode == 2) then   ! Ensemble Smoother analysis
! Random measurement perturbations
         call random(Ees,nrtobs*nrsamp)
         Ees=sqrt(vars%mes)*Ees
! Assemble measurments into one vector
         do i=1,nrmes
            avees((i-1)*nrobs+1:i*nrobs)=d(1:nrobs,i)
         enddo
! Observed ensemble prediction in Ses and ensemble of innovations in Des
         do j=1,nrsamp
            do i=1,nrmes
               nb=nint(deltaobs*real(i)/dt+1.0)
               Ses((i-1)*3+1:i*3,j)=A(:,nb,j)
            enddo
            Des(:,j)=avees(:)+Ees(:,j)-Ses(:,j)
         enddo
! Observed ensemble anomalies
         call ensmean(Ses,avees,nrtobs,nrsamp)
         do j=1,nrsamp
            ses(:,j)=ses(:,j)-avees(:)
         enddo
! Measurement error covariance
         Res=0.0
         do i=1,nrmes
            Res(i,i)=vars%mes
         enddo

         call analysis(A,Res,Ees,Ses,Des,innoves,3*ndim,nrsamp,nrtobs,.true.,0.9999,13,.false.,.false.,.true.,0,1.0,1)
   endif


   call ensmean(A,ave,3*ndim,nrsamp)
   call ensvar(A,ave,var,3*ndim,nrsamp)

   
   open(10,file=trim(fname)//'.dat')
      write(10,'(a)')'TITLE = "Lorenz output"'
      write(10,'(a)')'VARIABLES = "time" "Ref x" "Ref y" "Ref z" "Ave x" "Ave y" "Ave z" &
                                        &"Var x" "Var y" "Var z" "Err x" "Err y" "Err z"'
      write(10,*)'ZONE T="',trim(fname),'"  F=POINT, I=',ndim,', J=1, K=1'
      do i=1,ndim
         write(10,'(13g12.4)')float(i-1)*dt,ref(1:3,i),&
                                            ave(1:3,i),&
                                            sqrt(var(1:3,i)),&
                                            sqrt((ave(1:3,i)-ref(1:3,i))**2)
      enddo
      write(10,*)'ZONE T="Measurements" F=POINT, I=',nrmes,', J=1, K=1'
      do i=1,nrmes
         write(10,'(13g12.4)')float(i)*deltaobs, d(1:3,i), d(1:3,i), d(1:3,i), d(1:3,i) 
      enddo
   close(10)

end
