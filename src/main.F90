program main 
   use mod_dimensions       ! Parameter and type  declarations
   use m_ensmean
   use m_ensvar
   use m_random
   use m_lorenz
   use m_measurements
   use m_readinfile
   use m_set_random_seed2

   implicit none
   type(variances) vars
   integer na,nb,nc
   integer i,j,m

   real, allocatable :: A(:,:,:)    ! Stores the whole ensemble as function of time
   real, allocatable :: obs(:,:)    ! vector of measurements as function of time
   real, allocatable :: ave(:,:)    ! Ensemble average
   real, allocatable :: var(:,:)    ! Ensemble variance
   real, allocatable :: ref(:,:)    ! Reference solution
   real, allocatable :: fg(:)       ! first guess initial condition
   integer, allocatable :: ih(:)    ! Measurement time-index pointer

   real, allocatable :: R(:,:)      ! Measurement error covariance matrix
   real, allocatable :: D(:,:)      ! Perturbed measurement innovations D=obs+E-HA
   real, allocatable :: E(:,:)      ! Measurement perturbations
   real, allocatable :: S(:,:)      ! Predicted measurement anomalies
   real, allocatable :: aves(:)     ! work array for average of S
   real, allocatable :: innov(:)    ! mean of innovations D

   call set_random_seed2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading initial parameters
   call readinfile(vars,dt)

   allocate(A(neq,ndim,nrsamp), obs(nrobs,nrmes), ave(neq,ndim), var(neq,ndim), ref(neq,ndim), fg(neq), ih(0:nrmes))
   select case (mode)
   case(1,3)
      allocate(R(nrobs,nrobs),D(nrobs,nrsamp),E(nrobs,nrsamp),S(nrobs,nrsamp),aves(nrobs),innov(nrobs) )
   case(2)
      allocate(R(nrtobs,nrtobs), D(nrtobs,nrsamp), E(nrtobs,nrsamp), S(nrtobs,nrsamp), aves(nrtobs), innov(nrtobs))
   end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulating reference case
   ref(1,1)= 1.508870
   ref(2,1)=-1.531271
   ref(3,1)= 25.46091
   na=1
   nb=ndim
   call lorenz(na,nb,dt,ref,0.0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Defining first guess initial cond
   call random(fg,3)
   fg(:)=sqrt(vars%ini)*fg(:)+ref(:,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generating initial ensemble
   do j=1,nrsamp
      call random(A(1,1,j),3)
      A(:,1,j)=fg(:)+sqrt(vars%ini)*A(:,1,j)
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generating measurements from reference solution
   call measurements(ref,obs,ih,dt,vars,obsdt) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Time stepping over measurement times
   do m=1,nrmes

! Ensemble integration
      na=ih(m-1)
      nb=ih(m)
      print '(2(a,i4))','integrating over indexes: na=',na,' nb=',nb
      do j=1,nrsamp
         call lorenz(na,nb,dt,A(:,:,j),vars%dyn) 
      enddo
      print *,'Lorenz ensemble done'

! EnKF (lag=0) or EnKS update with given lag
      if ((mode == 1).or.(mode == 3)) then 
         call random(E,3*nrsamp)
         E=sqrt(vars%mes)*E
         do j=1,nrsamp
            S(:,j)=A(:,nb,j)
            D(:,j)=obs(:,m)+E(:,j)-S(:,j)
         enddo
         call ensmean(S,aves,3,nrsamp)
         do j=1,nrsamp
            S(:,j)=S(:,j)-aves(:)
         enddo

         call ensmean(D,innov,3,nrsamp)
         R=0.0


         na=max(1,nb-lag)
         nc=nb-na+1
         print '(3(a,i4))','analysis over indexes: ',na,'-',nb,', nc=',nc
         call analysis(A(:,na:nb,:),R,E,S,D,innov,neq*nc,nrsamp,nrobs,&
                       &.true.,0.99,13,.false.,.false.,.true.,0,1.0,1)
      endif

   enddo

   if (mode == 2) then   ! Ensemble Smoother analysis
! Random measurement perturbations
         call random(E,nrtobs*nrsamp)
         E=sqrt(vars%mes)*E

! Assemble measurments into one vector
         do m=1,nrmes
            aves((m-1)*nrobs+1:m*nrobs)=obs(1:nrobs,m)
         enddo

! Observed ensemble prediction in Ses and ensemble of innovations in Des
         do j=1,nrsamp
            do m=1,nrmes
               S((m-1)*nrobs+1:m*nrobs,j)=A(:,ih(m),j)
            enddo
            D(:,j)=aves(:)+E(:,j)-S(:,j)
         enddo

! Observed ensemble anomalies
         call ensmean(S,aves,nrtobs,nrsamp)
         do j=1,nrsamp
            S(:,j)=S(:,j)-aves(:)
         enddo

! Measurement error covariance
         R=0.0
         do i=1,nrtobs
            R(i,i)=vars%mes
         enddo

         call analysis(A,R,E,S,D,innov,neq*ndim,nrsamp,nrtobs,.true.,0.999,13,.false.,.false.,.true.,0,1.0,1)
   endif


   call ensmean(A,ave,3*ndim,nrsamp)
   call ensvar(A,ave,var,3*ndim,nrsamp)

   
   open(10,file=trim(fname)//'.dat')
      write(10,'(a)')'TITLE = "Lorenz output"'
      write(10,'(a)')'VARIABLES = "time" "Ref x" "Ref y" "Ref z" "Ave x" "Ave y" "Ave z" &
                                        &"Var x" "Var y" "Var z" "Err x" "Err y" "Err z"'
      write(10,*)'ZONE T="',trim(fname),'"  F=POINT, I=',ndim,', J=1, K=1'
      do i=1,ndim
         write(10,'(13g12.4)')float(i-1)*dt,ref(1:neq,i),&
                                            ave(1:neq,i),&
                                            sqrt(var(1:neq,i)),&
                                            sqrt((ave(1:neq,i)-ref(1:neq,i))**2)
      enddo
      write(10,*)'ZONE T="Measurements" F=POINT, I=',nrmes,', J=1, K=1'
      do i=1,nrmes
         write(10,'(13g12.4)')float(i)*obsdt, obs(1:nrobs,i), obs(1:nrobs,i), obs(1:nrobs,i), obs(1:nrobs,i) 
      enddo
   close(10)

end
