! using varying plate convergence rate. pay attention to theta !!!!
! yt(2*i) should not be negative!!!!
!!! with two buffer zones both sides of 20 km wide each
!! this is a restart code with locked seismogenic zone where slip==0
!------------------------------------------------
! D. Li
! Last modified Aug. 2018
! 
! SEAS benchmark BP3
! July 1st, 2021
!------------------------------------------------

program main
  USE mpi 
 
  USE phy3d_module_non
  implicit none
  integer, parameter :: DP=kind(1.d0)

  logical cyclecont
  integer ::Nt_all,Nt,ndt,ndtnext,kk,ii,n,l,ndt_v,ndt_inter,Itout,&  
       i,j,k,nrec,ihrestart,Ifileout, &
       Iperb,Isnapshot,&
       record_cor
  integer,dimension(:) :: s1(12)
  real (8) :: Vint,accuracy, epsv,dt_try, dt,dtmin,dt_did,dt_next, &
       t,tprint_inter, tint_out,tout,&
       tmin_out,tint_cos,tint_sse,&   
       tslip_ave,tslipend,tslip_aveint, tmax, &
       tslipsse,tslipcos, &
       tssestart,tsseend, &
       help

  real (8) ::  tmbegin,tmrun


  real (8), DIMENSION(:), ALLOCATABLE :: x,z,xi,yt,yt0,dydt,yt_scale, &
       slip,slipinc,slipds,sr,zzfric,zzfric2

  !Arrays only defined at master cpu
  real (8), DIMENSION(:), ALLOCATABLE :: x_all,xi_all,z_all,&
       yt_all,yt0_all,dydt_all,yt_scale_all,tau1_all,tau2_all, &
       slip_all,slipinc_all,slipds_all
  real(8), dimension(:),allocatable :: cca_all,ccb_all,xLf_all,seff_all,vi_all,&
        phy1_all,phy2_all

  !output related parameters
  integer :: imv,ias,icos,isse,Ioutput,inul,i_nul,n_nul_int
  real (8) :: vcos,vsse1,vsse2
  real (8), DIMENSION(:), ALLOCATABLE :: maxv,tmv,tas,tcos,tsse
  real (8),dimension(:,:,:),allocatable :: outs1
  real (8), DIMENSION(:,:), ALLOCATABLE :: slipz1_inter, &
       slipz1_cos, &
       slipz1_tau,slipz1_sse

  integer,DIMENSION(:),ALLOCATABLE :: ssetime
  real(8),DIMENSION(:,:),ALLOCATABLE :: slipz1_v

  character(len=40) :: cTemp,filename,ct

  integer:: n_obv, ndp
  real(8) :: dipangle
  integer,dimension(:),allocatable :: pdp
  real(8),dimension(:,:,:),allocatable :: obvs,obvdp
  real(8),dimension(:,:),allocatable :: surf1,surf2,surf3
  real(8) :: vel1,vel2,disp1,disp2

  !MPI RELATED DEFINITIONS
  integer :: ierr,size,myid,master

  call MPI_Init(ierr)
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, size, ierr )

  master = 0 

  !read in intiliazation parameters

  open(12,file='./parameter1.txt',form='formatted',status='old')

  read(12,'(a)')jobname
  read(12,'(a)')foldername
  read(12,'(a)')stiffname
  read(12,'(a)')restartname
  read(12,*)Nab,Nt_all,Nt,Lratio,nprocs,n_obv,ndp,dipangle
  read(12,*)Idin,Idout,Iprofile,Iperb,Isnapshot 
  read(12,*)Vpl
  read(12,*)tmax
  read(12,*)tslip_ave,tslipend,tslip_aveint
  read(12,*)tint_out,tmin_out,tint_cos,tint_sse
  read(12,*)vcos,vsse1,vsse2
  read(12,*)nmv,nas,ncos,nnul,nsse,n_nul_int
  read(12,*)s1(1),s1(2),s1(3),s1(4),s1(5),s1(6),s1(7),s1(8),s1(9),s1(10),s1(11),s1(12)
!!! modified data read in
110 format(A)
  close(12)

  !num. of rows each proc. 

  if(mod(Nt_all,nprocs)/=0)then
     write(*,*)'Nd_all must be integer*nprocs. Change nprocs!'
     STOP
  else
     write(*,*)'Each cpu calculates',Nt_all/nprocs,'cells'
  end if


  if(myid == master)then
     ALLOCATE(xi_all(Nt_all),&
          tau2_all(Nt_all),&
         slipds_all(Nt_all),&
        dydt_all(3*Nt_all),yt_scale_all(3*Nt_all))
     allocate(phy1_all(Nt_all),phy2_all(Nt_all))


     ALLOCATE (outs1(nmv,5,12),maxv(nmv), &
          tmv(nmv),tas(nas),tcos(ncos),tsse(nsse))

!!! modify output number
     ALLOCATE (slipz1_inter(Nt_all,nas),slipz1_cos(Nt_all,ncos), &
          slipz1_tau(Nt_all,nsse),slipz1_sse(Nt_all,nsse) )
     ALLOCATE(slipz1_v(Nt_all,ncos),ssetime(nsse)  )
    
!     allocate(seff_all(Nt_all),slipinc_all(Nt_all),tau1_all(Nt_all),yt0_all(3*Nt_all),yt_all(3*Nt_all),slip_all(Nt_all))
     allocate(obvs(nmv,4,n_obv),obvdp(nmv,3,ndp),pdp(ndp),surf1(n_obv,Nt_all),surf2(n_obv,Nt_all),surf3(n_obv,Nt_all))

 end if


     allocate(slipinc_all(Nt_all),tau1_all(Nt_all),yt0_all(3*Nt_all),yt_all(3*Nt_all),slip_all(Nt_all))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(phy1(Nt),phy2(Nt),zzfric(Nt),zzfric2(Nt),&
       x(Nt),z(Nt),z_all(Nt_all),xi(Nt),&
       tau1(Nt),tau2(Nt),tau0(Nt),slipds(Nt),slip(Nt),slipinc(Nt), &
       yt(3*Nt),dydt(3*Nt),yt_scale(3*Nt),yt0(3*Nt),sr(Nt))

  ALLOCATE (x_all(Nt_all),stiff(Nt,Nt_all),stiff2(Nt,Nt_all))   !!! stiffness of Stuart green calculation

  !Read in stiffness matrix, in nprocs segments

  write(cTemp,*) myid
  write(*,*) cTemp


  open(58, file=trim(stiffname)//'ssGreen_'//trim(adjustl(cTemp))//'.bin',form='unformatted',access='stream')

if(myid==master)then
  open(51,file=trim(stiffname)//'surfGreen'//'.bin',form='unformatted',access='stream')
  open(55,file=trim(stiffname)//'position.bin',form='unformatted',access='stream')
  open(56,file='prfdp'//jobname,form='formatted',status='old')
  do i=1,ndp
    read(56,*) pdp(i)
  end do
 close(56)
end if

  record_cor=Nt*(myid)

  !-------------------------------------------------------------------------------------------
  ! read stiffness from Stuart green calculation.
  !-----------------------------------------------------------------------------------------
  if(myid==master)then
     do k=1,Nt_all
        read(55) xi_all(k),x_all(k),z_all(k) !xi is along the fault-normal  while x is along the strike
        xi_all(k)=xi_all(k)/1000 ! meter to km
        x_all(k)=x_all(k)/1000
        z_all(k)=z_all(k)/1000

     end do

     do j = 1,n_obv ! kernel of surface displacement in x, y,z (downward)
      do k=1,Nt_all
        read(51) surf1(j,k)
        surf1(j,k)=-surf1(j,k)
      end do
      do k=1,Nt_all
        read(51) surf2(j,k)
        surf2(j,k)=-surf2(j,k)
      end do
      do k=1,Nt_all
        read(51) surf3(j,k)
      end do
    end do

  end if

  do i=1,Nt !! observe
     do j=1,Nt_all !! source
        read(58) stiff(i,j) ! kernel of shear traction
        stiff(i,j) = - stiff(i,j)
     end do
     do j=1,Nt_all !! source
        read(58) stiff2(i,j) ! kernel of normal traction
     end do
  end do
  close(58)
if(myid==0)then
  close(51)
  close(55)
end if
  !!-----------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------

 call MPI_Barrier(MPI_COMM_WORLD,ierr)

!  if(myid == master)then
      
     allocate(ccb_all(Nt_all),cca_all(Nt_all),xLf_all(Nt_all),seff_all(Nt_all),vi_all(Nt_all))

     CALL resdep(cca_all,ccb_all,xLf_all,seff_all,vi_all,Nt_all,jobname,foldername)
     write(*,*) ccb_all(Nt_all)
 ! end if


  allocate(ccb(Nt),cca(Nt),xLf(Nt),seff(Nt),vi(Nt))

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Bcast(z_all,Nt_all,MPI_Real8,master,MPI_COMM_WORLD,ierr)


  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  call MPI_Scatter(ccb_all,Nt,MPI_Real8,ccb,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)


  call MPI_Scatter(cca_all,Nt,MPI_Real8,cca,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
  call MPI_Scatter(xLf_all,Nt,MPI_Real8,xLf,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
  call MPI_Scatter(seff_all,Nt,MPI_Real8,seff,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
  call MPI_Scatter(vi_all,Nt,MPI_Real8,vi,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(x_all,Nt_all,MPI_Real8,master,MPI_COMM_WORLD,ierr) 
  call MPI_Scatter(x_all,Nt,MPI_Real8,x,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)


  call CPU_TIME(tmbegin)

  tm1=tmbegin
  tmday=86400.0
  tmelse=0.0
  tmmult=0.0
  tmmidn=0.0

  t=0.d0
  tprint_inter = 0.d0
  tslipcos = 0.d0
  tslipsse = 0.d0
  tout =0.d0

  ndtnext=0
  ndt_v=0
  ndt_inter=0
  ndt =0
  nrec=0
  i_nul = 0 
  inul = 0 
  isse = 0 

  Ioutput = 0    !initial always 0 (output)

  imv=0   !counter for maxv, sliptot, slipthresh1, slipthresh2,slipthresh3 output 
  ias=0 !counter for slip at iz3 and s.z. average slip output 
  icos = 0 

  accuracy = 1.d-4
  epsv = 1.0d-3
  dtmin = 1.d-10
  dt_try= dtmin
  Vint = Vpl


  if(myid==master)then
      open(311,file=trim(foldername)//'fltst_dp000'//jobname,access='append',status='unknown')
      open(312,file=trim(foldername)//'fltst_dp025'//jobname,access='append',status='unknown')
      open(313,file=trim(foldername)//'fltst_dp050'//jobname,access='append',status='unknown')
      open(314,file=trim(foldername)//'fltst_dp075'//jobname,access='append',status='unknown')
      open(315,file=trim(foldername)//'fltst_dp100'//jobname,access='append',status='unknown')
      open(316,file=trim(foldername)//'fltst_dp125'//jobname,access='append',status='unknown')
      open(317,file=trim(foldername)//'fltst_dp150'//jobname,access='append',status='unknown')
      open(318,file=trim(foldername)//'fltst_dp175'//jobname,access='append',status='unknown')
      open(319,file=trim(foldername)//'fltst_dp200'//jobname,access='append',status='unknown')
      open(320,file=trim(foldername)//'fltst_dp250'//jobname,access='append',status='unknown')
      open(321,file=trim(foldername)//'fltst_dp300'//jobname,access='append',status='unknown')
      open(322,file=trim(foldername)//'fltst_dp350'//jobname,access='append',status='unknown')

       do i=311,322
   		write(i,100)'# This is the file header'
		write(i,100)'# problem=SEAS Benchmark No.3; dip=30'
		write(i,100)'# author=D.Li'
                write(i,100)'# code=TriBIE'
		write(i,100)'# date=2021/5/11'
		write(i,100)'# element_size = 100 m'
		write(i,100)'# minimum_time_step = 1e-3'
		write(i,100)'# maximum_time_step = 2e+7'
		write(i,100)'# location = on fault: file name'
		write(i,100)'# Column #1 = Time (s)'
		write(i,100)'# Column #2 = slip (m)'
		write(i,100)'# Column #3 = Slip_rate (log10 m/s)'
		write(i,100)'# Column #4 = Shear stress  (MPa)'
		write(i,100)'# Column #5 = normal stress  (MPa)'
		write(i,100)'# Column #6 = State (log10 s)'
		write(i,100)'# '
		write(i,100)'# The line below lists the names of the data fields:'
		write(i,'(A,1x,A,1x,A,1x,A,1x,A,1x,A,1x)')'t','slip','slip_rate','shear_stress','normal_stress','state'
		write(i,100)'# Below is the time-series data.'		
	end do 
 100	format(A)

       open(30,file=trim(foldername)//'maxvall'//jobname,access='append',status='unknown')
        i=30
        write(i,100)'# This is the file header'
        write(i,100)'# problem=SEAS Benchmark No.3'
        write(i,100)'# author=D.Li'
        write(i,100)'# code=TriBIE'
        write(i,100)'# date=2021/5/11'
        write(i,100)'# element_size = 100 m'
        write(i,100)'# minimum_time_step = 1e-3'
        write(i,100)'# maximum_time_step = 2e+7'
        write(i,100)'# Column #1 = Time (s)'
        write(i,100)'# Column #2 = max_slip_rate (logV m/s)'
        write(i,100)'# Column #3 = slp (m)'
        write(i,100)'# The line below lists the names of the data fields:'
        write(i,'(A,1x,A,1x,A)')'t','max_slip_rate','slip'
        write(i,100)'# Below is the time-series data.'  

      close(30)
      do j=311,322
         close(j)
      end do

      open(401,file=trim(foldername)//'srfst_fn-32'//jobname,access='append',status='unknown')
      open(402,file=trim(foldername)//'srfst_fn-16'//jobname,access='append',status='unknown')
      open(403,file=trim(foldername)//'srfst_fn-08'//jobname,access='append',status='unknown')
      open(404,file=trim(foldername)//'srfst_fn-00'//jobname,access='append',status='unknown')
      open(405,file=trim(foldername)//'srfst_fn+08'//jobname,access='append',status='unknown')
      open(406,file=trim(foldername)//'srfst_fn+16'//jobname,access='append',status='unknown')
      open(407,file=trim(foldername)//'srfst_fn+32'//jobname,access='append',status='unknown')
      !open(408,file=trim(foldername)//'srfst_fn+00'//jobname,access='append',status='unknown')

    do i = 401,407
        write(i,100)'# This is the file header'
        write(i,100)'# problem=SEAS Benchmark No.3'
        write(i,100)'# author=D.Li'
        write(i,100)'# code=TriBIE'
        write(i,100)'# date=2021/7/01'
        write(i,100)'# element_size = 100 m'
        write(i,100)'# minimum_time_step = 1e-3'
        write(i,100)'# maximum_time_step = 2e+7'
        write(i,100)'# location = off fault: file name'
        write(i,100)'# Column #1 = Time (s)'
        write(i,100)'# Column #2 = displacement 1 (m)'
        write(i,100)'# Column #3 = displacement 2 (m)'
        write(i,100)'# Column #4 = velocity_1 (m/s)'
        write(i,100)'# Column #5 = velocity_2 (m/s)'
        write(i,100)'# The line below lists the names of the data fields:'
        write(i,'(A,1x,A,1x,A,1x,A,1x,A)')'t','disp_1','disp_2','vel_1','vel_2'
        write(i,100)'# Below is the time-series data.'
    end do

    do i = 401,407
      close(i)
    end do

   open(501,file=trim(foldername)//'slip'//jobname,access='append',status='unknown')
   open(502,file=trim(foldername)//'shear_stress'//jobname,access='append',status='unknown')
   open(503,file=trim(foldername)//'normal_stress'//jobname,access='append',status='unknown')

  do i=501,503
        write(i,100)'# This is the file header'
        write(i,100)'# problem=SEAS Benchmark No.3'
        write(i,100)'# author=D.Li'
        write(i,100)'# code=TriBIE'
        write(i,100)'# date=2021/7/1'
        write(i,100)'# element_size = 100 m'
        write(i,100)'# Row #1 downdip distance with two zeros'
        write(i,100)'# Column #1 = Time (s)'
        write(i,100)'# Column #2 = max slip rate (log10 m/s)'
        write(i,100)'# Column #3-83 = slip (m)or shear(or normal) stress (Mpa)'
        write(i,100)'# The line below lists the names of the data fields:'
        write(i,100)'x_d'
        write(i,100)'t'
        write(i,100)'max_slip_rate'
        write(i,100)'slip'
        write(i,100)'# Below is the time-series data.'
  end do
122 format(83(E15.7,1X))
      write(501,122) 0d0,0d0, z_all(pdp(:))/dsin(dipangle/180*pi)*1d3
      write(502,122) 0d0,0d0, z_all(pdp(:))/dsin(dipangle/180*pi)*1d3
      write(503,122) 0d0,0d0, z_all(pdp(:))/dsin(dipangle/180*pi)*1d3

    do i=501,503
      close(i)
    end do

     open(1,file=trim(foldername)//'phypara'//jobname, status='unknown')
     write(1,*)'procs num = ', nprocs

  end if

  Ifileout = 60   !file index, after 47

  !----Initial values of velocity, state variable, shear stress and slip--
  !--SET INITIAL VPL FOR THE LOCKED PART TO BE 0 
  ! ! set plate convergence
  if(IDin.eq.0)then 
     disp1=0.0
     disp2=0.0
     do j=1,Nt

        yt(3*j-2)=Vint

        help=(Vint/(2.0*V0))*dexp((f0+ccb(j)*dlog(V0/Vint))/amax)
        tau0(j)=seff(j)*amax*dlog(help+dsqrt(1+help**2))+ eta*Vint
	    help = dlog((2.0*V0/Vint)*dsinh((tau0(j) - eta*Vint)/(cca(j)*seff(j))))
        yt(3*j-1) = seff(j)
        yt(3*j) = (xLf(j)/V0)*exp((cca(j)/ccb(j))*help - f0/ccb(j))
        slip(j)=0.d0
        yt0(3*j-2)=yt(3*j-2)
        yt0(3*j) = yt(3*j)
        tau1(j) = tau0(j)
     end do
  end if

  !------------------------------------------------------------------
  if(IDin.eq.1) then               !if this is a restart job
     if(myid==master)then
        call restart(0,'out',4,Nt_all,t,dt,dt_try,ndt,nrec,yt_all,slip_all)
        write(1,*)'This is a restart job. Start time ',t,' yr'
     end if
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_Bcast(t,1,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(dt,1,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(dt_try,1,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ndt,1,MPI_integer,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(nrec,1,MPI_integer,master,MPI_COMM_WORLD,ierr)
     call MPI_Scatter(yt_all,3*Nt,MPI_Real8,yt,3*Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Scatter(slip_all,Nt,MPI_Real8,slip,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)


     ndtnext = ndt
     tprint_inter = t
     tslip_ave=t 	
     tout = t

  else
     if(myid==master)then
        write(1,*)'Start time ',t,' yr'
     end if
  end if
  if(myid==master)then
     close(1)
  end if

  !----------------------------------------------
  !      Start of Basic Cycle:
  !----------------------------------------------
  cyclecont=.true.

  do while(cyclecont) 

     call derivs(myid,dydt,3*Nt,Nt_all,Nt,t,yt,z_all,x) 

     do j=1,3*Nt
        yt_scale(j)=dabs(yt(j))+dabs(dt_try*dydt(j))
        yt0(j) = yt(j)
     end do
     CALL rkqs(myid,yt,dydt,3*Nt,Nt_all,Nt,t,dt_try,accuracy,yt_scale, &
          dt_did,dt_next,z_all,x)

     dt = dt_did
     dt_try = dt_next

     do i=1,Nt
       if(z_all(Nt*myid+i).ge.40.0*dsin(dipangle/180*PI))then
         yt0(3*i-2)=Vpl
         yt(3*i-2)=Vpl
        end if
        help=(yt(3*i-2)/(2*V0))*dexp((f0+ccb(i)*dlog(V0*yt(3*i)/xLf(i)))/cca(i))
        seff(i) = yt(3*i-1)
        slipinc(i) = 0.5*(yt0(3*i-2)+yt(3*i-2))*dt
        slip(i) = slip(i) + slipinc(i)
        tau1(i) = seff(i)*cca(i)*dlog(help+dsqrt(1+help**2))

     end do

     ndt = ndt + 1



     call MPI_Barrier(MPI_COMM_WORLD,ierr)

     call MPI_Gather(yt0, 3*Nt,MPI_real8,yt0_all,3*Nt,MPI_real8,0,MPI_COMM_WORLD,ierr)
     call MPI_Gather(yt,3*Nt,MPI_real8,yt_all,3*Nt,MPI_real8,0,MPI_COMM_WORLD,ierr)

     call MPI_Gather(slipinc,Nt,MPI_Real8,slipinc_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Gather(slip,Nt,MPI_Real8,slip_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Gather(seff,Nt,MPI_Real8,seff_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Gather(tau1,Nt,MPI_Real8,tau1_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)

     !-------------------
     !      Output:     (a single thread will do the writing while others 
     !      proceed 
     !-------------------

     if(myid==master)then
        imv=imv+1
        tmv(imv)=t*yrs
        maxv(imv) = 0.d0
        do i=1,Nt_all      !!!!!!!!!!!!!!!!! find max velocity
           if(yt_all(3*i-2).ge.maxv(imv))then
              maxv(imv)=yt_all(3*i-2)
              ! maxnum(imv)=i
           end if
        end do
!!!!!  SEAS output variables
       
       do i = 1,12
        outs1(imv,1,i) = slip_all(s1(i))*1.d-3 ! meter
        outs1(imv,2,i) =  dlog10(yt_all(3*s1(i)-2)*1.d-3/yrs) ! log10(V) m/s
        outs1(imv,3,i) = tau1_all(s1(i))/10 ! MPa
        outs1(imv,4,i) = seff_all(s1(i))/10
        outs1(imv,5,i) = dlog10(yt_all(3*s1(i))*yrs) ! log10(theta)
       end do

       do i=1,ndp
        obvdp(imv,1,i)=slip_all(pdp(i))*1d-3
        obvdp(imv,2,i)=tau1_all(pdp(i))/10
        obvdp(imv,3,i)=seff_all(pdp(i))/10
       end do

 do i = 1,n_obv
          vel1=0d0
          vel2=0d0
          disp1=0d0
          disp2=0d0
         do j=1,Nt_all
            vel1 = vel1 + surf1(i,j)*(yt0_all(3*j-2)+yt_all(3*j-2))*0.5
            vel2 = vel2 + surf3(i,j)*(yt0_all(3*j-2)+yt_all(3*j-2))*0.5
          disp1=disp1+surf1(i,j)*slip_all(j)
          disp2=disp2+surf3(i,j)*slip_all(j)
         end do
       obvs(imv,3,i) = vel1/1d3/yrs
       obvs(imv,4,i) = vel2/1d3/yrs
       obvs(imv,1,i) = disp1/1d3
       obvs(imv,2,i) = disp2/1d3
     end do
        !-----Interseismic slip every ? years----

        if (t.ge.tslip_ave)then
           ias = ias + 1 
           tas(ias)=t

           do i=1,Nt_all
              slipz1_inter(i,ias) = slip_all(i)*1.d-3           
           end do

            tslip_ave = tslip_ave + tslip_aveint
        end if

        !------velocity and slip of eq nucleation process -- 


        !----SSE slip ------
        if(t.ge.tssestart.and.t.le.tsseend)then
        end if

        !------coseismic Slip   -------------------

        if((maxv(imv)/yrs).ge.vcos)then
           tslipcos = tslipcos+dt
           if(tslipcos.ge.tint_cos)then    !!!!tint_cos: every 5s output 
		write(*,*) t,dlog10(maxv(imv)*1d-3/yrs)

130 format(E20.13,1X,E15.7)

              icos = icos +1
              tcos(icos) = t 


              do i=1,Nt_all
                 slipz1_cos(i,icos) = slip_all(i)*1.d-3
                 slipz1_v(i,icos) = dlog10(yt_all(3*i-2)*1.d-3/yrs) 
              end do

              tslipcos = 0.d0
           end if
        end if
     end if

     !----Output restart files -------------
     if(myid==master)then
        if(mod(ndt,1000).eq.0)ihrestart=1
        if(IDout.eq.1.and.ihrestart.eq.1)then
	   filename='out0'
           call restart(1,filename,4,Nt_all,t,dt,dt_try,ndt,nrec,yt_all,slip_all)  
           ihrestart=0
        end if
        if(abs(t-tout).le.tmin_out)then
           ihrestart = 1
           Itout=int(tout)
           write(ct,*)Itout
           ct=adjustl(ct)
           filename='out'//trim(ct)
           call restart(1,filename,Ifileout,Nt_all,t,dt,dt_try,ndt,nrec,yt_all,slip_all)
           ihrestart = 0
           tout = tout+tint_out 
        end if
     end if

     !----Output velocity and slip records --- 
     !----velocity in mm/yr ; slip in meter, moment in 10^{14} Nm -- 
     if(myid==master)then 
        Ioutput = 0 
        call output(Ioutput,Isnapshot,Nt_all,Nt,inul,imv,ias,icos,isse,x,&
             tmv,tas,tcos,tsse,maxv,outs1, &
             slipz1_inter,slipz1_tau,slipz1_sse, &
             slipz1_cos,&
             slipz1_v,obvs,n_obv,obvdp,ndp)         

     end if

     !      end of output 
     !----------------------------------------------------
     !   End of output conmmands.
     !----------------------------------------------------  

     !----------------------------------------------------
     !   Go to next cycle if t < tmax or ndt > ndtmax 
     !---------------------------------------------------

     if (t>tmax)cyclecont = .false.

  end do

  !--- Final output ------- 

  if(myid==master)then 
     filename='outlast'
     call restart(1,filename,Ifileout,Nt_all,t,dt,dt_try,ndt,nrec,yt_all,slip_all)
     Ioutput = 1
     call output(Ioutput,Isnapshot,Nt_all,Nt,inul,imv,ias,icos,isse,x,&
          tmv,tas,tcos,tsse,maxv,outs1,&
          slipz1_inter,slipz1_tau,slipz1_sse, &
          slipz1_cos,&
         slipz1_v,obvs,n_obv,obvdp,ndp) 

  end if

  if(myid==master)then
        i=410  
        open(i,file=trim(foldername)//'rupture'//jobname,status='unknown')
        write(i,100)'# This is the file header'
        write(i,100)'# problem=SEAS Benchmark No.5'
        write(i,100)'# author=D.Li '
        write(i,100)'# code=TriBIE'
        write(i,100)'# date=2021/5/11'
        write(i,100)'# element_size = 100 m'
        write(i,100)'# Column #1 = x2 (m)'
        write(i,100)'# Column #2 = x3 (m)'
        write(i,100)'# Column #3 = t (s)'
        write(i,100)'# '
        write(i,100)'# The line below lists the names of the data fields:'
        write(i,'(A,1x,A,1x,A)')'x2','x3','t'
        write(i,100)'# Below is the time-series data.'      

        close(i)
111 format(3(E22.14,1x))
end if

  !---End of final output ----

  call CPU_TIME(tmrun)
  tmrun = tmmidn*tmday + tmrun - tmbegin
  tmmult = tmmult/tmrun*100.
  tmelse = tmelse/tmrun*100.
  if(myid==master)then
     open(10,file=trim(foldername)//'summary'//jobname,status='unknown')

     write(10,*)'processor',myid
     write(10,*)'      PARTS OF RUNNING TIME     '
     write(10,'(T9,A,T40,F10.5)')'message passing (percent)', tmmult
     write(10,'(T9,A,T40,F10.5)')'everything else (percent)', tmelse
     write(10,*)
     write(10,'(T9,A,T45,F20.1)')'total running time(min)',tmrun*0.016667
     write(10,'(T9,A,T45,F20.1)')'total running time(hr)',tmrun*0.016667*0.016667
     write(10,*)
     write(10,*)'       INFORMATION ABOUT THE END OF THE RUN      '
     write(10,'(T9,A,I7)')'ndtend = ', ndt
     write(10,'(T9,A,D20.13)')'tend = ',t
     write(10,'(T9,A,I7)')'nrec = ',nrec
     write(10,*)       
     write(10,*)'Nprocs=', nprocs
     close(10)
  end if


  if(myid==master)then 
     DEALLOCATE (x_all,xi_all,dydt_all,yt_scale_all,&
                phy1_all,phy2_all,tau2_all, &
          slipds_all,&
          cca_all,ccb_all,vi_all,xLf_all,seff_all, &
          maxv,outs1,tmv,tcos,tas,tsse)

     DEALLOCATE (slipz1_inter,slipz1_tau,slipz1_sse,slipz1_cos)
     DEALLOCATE (ssetime,slipz1_v)
     deallocate (pdp,obvs,obvdp,surf1,surf2,surf3)
  end if


  DEALLOCATE (stiff2,stiff)
  deallocate (slipinc_all, yt_all,  tau1_all, slip_all)
  DEALLOCATE (x,z_all,xi,yt,dydt,yt_scale)
  deallocate (phy1,phy2,tau1,tau2,tau0,slip,slipinc,slipds,yt0,zzfric,zzfric2)
  DEALLOCATE (cca,ccb,xLf,seff,vi,sr)
  call MPI_finalize(ierr)

contains

!-----------------------------------------------------------------------------
!    read parameters: sigma_effective, a,b,D_c
!----------------------------------------------------------------------------

    subroutine resdep(cca_all,ccb_all,xLf_all,seff_all,vi_all,Nt_all,jobname,foldername)
      
      real(8),intent(inout) :: cca_all(:),ccb_all(:),xLf_all(:),vi_all(:),seff_all(:)
      character(len=80),intent(in) :: jobname, foldername
      integer,intent(in) :: Nt_all
      integer :: i

      ! allocate(ccb_all(Nt_all),cca_all(Nt_all),xLf_all(Nt_all),seff_all(Nt_all),vi_all(Nt_all))

      open(444,file='var'//jobname,status='old', action='read',form='formatted')
       do i=1,Nt_all
        read(444,*) seff_all(i),xLf_all(i),cca_all(i),ccb_all(i),vi_all(i)
        vi_all(i) = vi_all(i)*yrs*1d3
       end do
      close(444)

      RETURN
    END subroutine resdep


END program main

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine rkqs(myid,y,dydx,n,Nt_all,Nt,x,htry,eps,yscal,hdid,hnext,z_all,p)
  Use mpi
  USE phy3d_module_non, only : nprocs,phi
  implicit none
  integer, parameter :: DP = kind(1.0d0)   
  integer :: n,i,j,k,NMAX,Nt,Nt_all
  real (8) :: eps,hdid,hnext,htry,x
  real (8) :: dydx(n),y(n),yscal(n),z_all(Nt_all),p(Nt) !p is position
  external derivs
  real (8) :: errmax,errmax1,h,htemp,xnew,errmax_all(nprocs)
  real (8), dimension(:), allocatable :: yerr,ytemp
  real (8), parameter :: SAFETY=0.9, PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4

  !MPI RELATED DEFINITIONS
  integer :: ierr,myid,master
  master = 0 

  nmax=n
  h=htry
  allocate (yerr(nmax),ytemp(nmax))
  !        write(*,*)'in rkqs, bf rkck, myid=',myid,y(n-1)
1 call rkck(myid,dydx,h,n,Nt_all,Nt,y,yerr,ytemp,x,derivs,z_all,p)
  !        write(*,*)'in rkqs, af rkck, myid=',myid, ytemp(n-1)
  errmax=0.
  do  i=1,nmax
     j = int(ceiling(real(i)/3)) ! position within central part
     !   if(p(j).gt.-700.and.p(j).lt.180)then
     errmax = max(errmax,dabs(yerr(i)/yscal(i)))
  end do
  errmax=errmax/eps
!  write(*,*) errmax
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  call MPI_Gather(errmax,1,MPI_Real8,errmax_all,1,MPI_Real8,master,MPI_COMM_WORLD,ierr)

  if(myid==master)then
     errmax1=maxval(errmax_all)
  end if
  CALL MPI_BCAST(errmax1,1,MPI_REAL8,master,MPI_COMM_WORLD, ierr)

  if(errmax1.gt.1.)then
     htemp = SAFETY*h*(errmax1**PSHRNK)
     h = dsign(max(dabs(htemp),0.1*dabs(h)),h)
     xnew = x+h
     if(xnew.eq.x) stop !pause 'stepsize underflow in rkqs'
     goto 1
  else
     if(errmax1.gt.ERRCON)then
        hnext=SAFETY*h*(errmax1**PGROW)
     else
        hnext=5.*h
     end if
     hdid=h
     x=x+h 
     do i=1,nmax
        y(i)=ytemp(i)
     end do
  end if

  deallocate (yerr,ytemp)

  RETURN
end subroutine rkqs
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
     subroutine rkck(myid,dydx,h,n,Nt_all,Nt,y,yerr,yout,x,derivs,z_all,p)
       USE phy3d_module_non, only :nprocs
       implicit none
       integer, parameter :: DP = kind(1.0d0)   
       integer :: n,i,NMAX,myid,Nt_all,Nt
       external derivs
       real (8) :: h,x,dydx(n),y(n),yerr(n),yout(n),z_all(Nt_all),p(Nt)
       real (8), dimension(:), ALLOCATABLE :: ak2,ak3,ak4,ak5,ak6,ytemp
       REAL (8),  parameter :: A2=.2,A3=.3,A4=.6,A5=1.,A6=.875, &
            B21=.2,B31=3./40.,B32=9./40.,B41=.3,&
            B42=-.9,B43=1.2,B51=-11./54.,B52=2.5, &
            B53=-70./27.,B54=35./27., B61=1631./55296., &
            B62=175./512.,B63=575./13824.,B64=44275./110592., &
            B65=253./4096.,C1=37./378., C3=250./621.,  &
            C4=125./594.,C6=512./1771.,DC1=C1-2825./27648., &
            DC3=C3-18575./48384.,DC4=C4-13525./55296.,  &
            DC5=-277./14336.,DC6=C6-.25


       nmax = n
       ALLOCATE (ak2(nmax),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),ytemp(NMAX))

       do  i=1,n
          ytemp(i)=y(i)+B21*h*dydx(i)
       end do
       call derivs(myid,ak2,n,Nt_all,Nt,x+A2*h,ytemp,z_all,p)
       do  i=1,n
          ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
       end do
       call derivs(myid,ak3,n,Nt_all,Nt,x+A3*h,ytemp,z_all,p)
       do  i=1,n
          ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
       end do
       call derivs(myid,ak4,n,Nt_all,Nt,x+A4*h,ytemp,z_all,p)
       do  i=1,n
          ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+  &
               B54*ak4(i))
       end do
       call derivs(myid,ak5,n,Nt_all,Nt,x+A5*h,ytemp,z_all,p)
       do  i=1,n
          ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+  &
               B64*ak4(i)+B65*ak5(i))
       end do
       call derivs(myid,ak6,n,Nt_all,Nt,x+A6*h,ytemp,z_all,p)
       do  i=1,n
          yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+   &
               C6*ak6(i))
       end do
       do  i=1,n
          yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+     &
               DC5*ak5(i)+DC6*ak6(i))
       end do
       DEALLOCATE (ak2,ak3,ak4,ak5,ak6,ytemp)
       return
     end subroutine rkck
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
     subroutine derivs(myid,dydt,nv,Nt_all,Nt,t,yt,z_all,x)
       USE mpi
       USE phy3d_module_non, only: phi,phy1,phy2,tau1,tau2, stiff2,stiff,cca,ccb,seff,xLf,eta,f0,Vpl,V0,Lratio,nprocs,pi,&
            tm1,tm2,tmday,tmelse,tmmidn,tmmult
       implicit none
       integer, parameter :: DP = kind(1.0d0)
       integer :: nv,n,i,j,k,kk,l,ii,Nt,Nt_all
       real (8) :: t,yt(nv),dydt(nv),frc 
       real (8) :: deriv3,deriv2,deriv1,small
       real (8) :: psi,help1,help2,help
       real (8) :: SECNDS
       real (8) :: sr(Nt),z_all(Nt_all),x(Nt),zz(Nt),zzfric(Nt),zz_all(Nt_all),zzfric2(Nt)
       intrinsic imag,real

       !MPI RELATED DEFINITIONS
       integer :: ierr,myid,master
       master = 0 

       small=1.d-6
        
       do i=1,Nt
          zz(i)=yt(3*i-2)-Vpl
         if(z_all(Nt*myid+i).ge.40.0*dsin(90.0/180.0*pi))then
             zz(i)=0.0
             yt(3*i-2)=Vpl
          end if
       end do

       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_Gather(zz,Nt,MPI_Real8,zz_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(zz_all,Nt_all,MPI_Real8,master,MPI_COMM_WORLD,ierr)
       !----------------------------------------------------------------------
       !    summation of stiffness of all elements in slab
       !----------------------------------------------------------------------
       ! initilize zzfric
       call CPU_TIME(tm2)
       tmelse=tmelse+tm2-tm1
       tm1=tm2

!!!!!!!!!!!!!!!!!!!!

       ! calculate stiffness*(Vkl-Vpl) of one proccessor.
       do i=1,Nt
          zzfric(i)=0.d0
          zzfric2(i)=0.d0
     
          do j=1, Nt_all
             ! number in all element, total Nt_all
             zzfric(i)=zzfric(i) + stiff(i,j)*zz_all(j) ! shear
             zzfric2(i)=zzfric2(i)+stiff2(i,j)*zz_all(j) ! normal-
          end do
          
       end do

       call CPU_TIME(tm2)
       if ((tm2-tm1) .lt. 0.03)then
          tmmult=tmmult+tm2-tm1
       else
          tmmidn=tmmidn+1
          tmmult=tmmult+tm2-tm1+tmday
       end if
       tm1=tm2

       do i=1,Nt
          frc = f0+cca(i)*dlog(yt(3*i-2)/V0)+ccb(i)*dlog(V0*yt(3*i)/xLf(i))
          psi = dlog(V0*yt(3*i)/xLf(i))
          help1 = yt(3*i-2)/(2*V0)
          help2 = (f0+ccb(i)*psi)/cca(i)
          help = dsqrt(1+(help1*dexp(help2))**2)

          deriv1 = (yt(3*i-1)*ccb(i)/yt(3*i))*help1*dexp(help2)/help
          deriv2 = (yt(3*i-1)*cca(i)/(2*V0))*dexp(help2)/help
!aging        
	      deriv3 = 1-yt(3*i-2)*yt(3*i)/xLf(i)
          dydt(3*i-1) = zzfric2(i)
!slip law	     deriv3 = -yt(2*i-1)*yt(2*i)/xLf(i)*dlog(yt(2*i-1)*yt(2*i)/xLf(i))
          dydt(3*i-2) = -(zzfric(i)+deriv1*deriv3+frc*dydt(3*i-1))/(eta+deriv2) ! total shear traction
          !dydt(3*i-2) = dydtinc*yt(3*i-2)/sr(i)
          !dydt(3*i-1) = dydtinc*yt(3*i-1)/sr(i)
          !dydt(3*i-1) = -(zzfric2(i)+deriv1*deriv3)/(eta+deriv2) 
          dydt(3*i)=deriv3     
!aging       
!slip law	     deriv3 = -yt(2*i-1)*yt(2*i)/xLf(i)*dlog(yt(2*i-1)*yt(2*i)/xLf(i))
       end do

!       if(myid==0) write(*,*) t,dydt(1),dydt(2)

       RETURN
     END subroutine derivs

!-----------------------------------------------------------------------------
!    read parameters: sigma_effective, a,b,D_c
!----------------------------------------------------------------------------

    subroutine resdep_back(cca_all,ccb_all,xLf_all,seff_all,vi_all,Nt_all)
      USE phy3d_module_non, only: yrs,jobname,foldername
      implicit none
      real(8),allocatable :: cca_all(:),ccb_all(:),xLf_all(:),vi_all(:),seff_all(:)
 
      integer, intent(in) :: Nt_all
      integer :: i
 
      allocate(ccb_all(Nt_all),cca_all(Nt_all),xLf_all(Nt_all),seff_all(Nt_all),vi_all(Nt_all))

      !allocate(cca_all(Nt_all),ccb_all(Nt_all),seff_all(Nt_all),xLf_all(Nt_all),vi_all(Nt_all))
      !----------------------------------------------------------------------------
      !     iseff defines what eff. normal stress down-dip profiles
      !     1:     linearly increase to sigmadiff and keep constant
      !     2:     linearly increase to sigmadiff, followed by a sudden drop
      !               to a much lower level of Lffix
      !     3:     linearly increase to sigmadiff, followed by a sudden drop
      !               to Lffix over certain range, then resume sigmadiff at downdip
      !     4:     other profiles to be defined (?)
      !-----------------------------------------------------------------------------

!!! set SSE depth effective normal stress and Dc
!!!! add perturbation and buffer zone at both ends.

 !need to address when j=1 and j=Nd_all!! same in the old openmp f90 file!

      open(444,file='var'//jobname,status='old')
       do i=1,Nt_all
        read(444,*) seff_all(i),xLf_all(i),cca_all(i),ccb_all(i),vi_all(i)
        vi_all(i) = vi_all(i)*yrs*1d3
       end do
      close(444)

      !     To save info about some of the quantities
      RETURN
    END subroutine resdep_back
!       
!------------------------------------------------------------------------------
! restart file
!------------------------------------------------------------------------------

      subroutine restart(inout,filename,Ifileout,Nt_all,t,dt,dt_try,ndt,nrec,yt,slip)
USE phy3d_module_non, ONLY : jobname,foldername,restartname, &
                        tm1,tm2,tmday,tmelse,tmmidn,tmmult
      implicit none
      integer, parameter :: DP = kind(1.0d0)
      integer :: inout,i,ndt,nrec,Ifileout,Nt,Nt_all
      real (8) :: t,dt,dt_try
      real (8) ::  yt(3*Nt_all),slip(Nt_all)
character(len=40) :: filename

      if(inout.eq.0) then
         open(Ifileout,file=trim(restartname),status='old')
          read(Ifileout,*)t,ndt,nrec
          read(Ifileout,*)dt,dt_try
          do i=1,3*Nt_all
             read(Ifileout,*)yt(i)
          end do

          do i=1,Nt_all
             read(Ifileout,*)slip(i)
          end do
	  close(Ifileout)
      else
         open(Ifileout,file=trim(foldername)//trim(filename)//jobname,status='unknown')
         write(Ifileout,*)t,ndt,nrec
         write(Ifileout,*)dt,dt_try
         do i=1,3*Nt_all
              write(Ifileout,*)yt(i)
         end do
         do i=1,Nt_all
              write(Ifileout,*)slip(i)
         end do
         close(Ifileout)
        end if

      RETURN
      END

!------Output -------------------------------------------
!--------------------------------------------------------
subroutine output(Ioutput,Isnapshot,Nt_all,Nt,inul,imv,ias,icos,isse,x,&
    tmv,tas,tcos,tsse,maxv,outs1,&
     slipz1_inter,slipz1_tau,slipz1_sse,&
     slipz1_cos,&
    slipz1_v,obvs,n_obv,obvdp,ndp) 


USE phy3d_module_non, only: xmu,nmv,nas,ncos,nnul,nsse,yrs,Vpl,Nl, &
		foldername,jobname
implicit none
integer, parameter :: DP = kind(1.0d0)
integer :: Nt,Nt_all,i,j,k,l,kk,inul,imv,ias,icos,isse,Ioutput,Isnapshot

real (8) :: x(Nt),maxv(nmv),outs1(nmv,5,12),&
	tmv(nmv),tas(nas),tcos(ncos),tsse(nsse)

real (8) :: slipz1_inter(Nt_all,nas),slipz1_cos(Nt_all,ncos),&
        slipz1_tau(Nt_all,nsse),slipz1_sse(Nt_all,nsse), &
      slipz1_v(Nt_all,ncos)
real(8):: obvs(nmv,4,n_obv),obvdp(nmv,3,ndp)
integer :: n_obv, ndp

if(Ioutput == 0)then    !output during run 


   if(imv==nmv)then
      open(30,file=trim(foldername)//'maxvall'//jobname,access='append',status='unknown')
      open(311,file=trim(foldername)//'fltst_dp000'//jobname,access='append',status='unknown')
      open(312,file=trim(foldername)//'fltst_dp025'//jobname,access='append',status='unknown')
      open(313,file=trim(foldername)//'fltst_dp050'//jobname,access='append',status='unknown')
      open(314,file=trim(foldername)//'fltst_dp075'//jobname,access='append',status='unknown')
      open(315,file=trim(foldername)//'fltst_dp100'//jobname,access='append',status='unknown')
      open(316,file=trim(foldername)//'fltst_dp125'//jobname,access='append',status='unknown')
      open(317,file=trim(foldername)//'fltst_dp150'//jobname,access='append',status='unknown')
      open(318,file=trim(foldername)//'fltst_dp175'//jobname,access='append',status='unknown')
      open(319,file=trim(foldername)//'fltst_dp200'//jobname,access='append',status='unknown')
      open(320,file=trim(foldername)//'fltst_dp250'//jobname,access='append',status='unknown')
      open(321,file=trim(foldername)//'fltst_dp300'//jobname,access='append',status='unknown')
      open(322,file=trim(foldername)//'fltst_dp350'//jobname,access='append',status='unknown')

      open(401,file=trim(foldername)//'srfst_fn-32'//jobname,access='append',status='unknown')
      open(402,file=trim(foldername)//'srfst_fn-16'//jobname,access='append',status='unknown')
      open(403,file=trim(foldername)//'srfst_fn-08'//jobname,access='append',status='unknown')
      open(404,file=trim(foldername)//'srfst_fn-00'//jobname,access='append',status='unknown')
      open(405,file=trim(foldername)//'srfst_fn+08'//jobname,access='append',status='unknown')
      open(406,file=trim(foldername)//'srfst_fn+16'//jobname,access='append',status='unknown')
      open(407,file=trim(foldername)//'srfst_fn+32'//jobname,access='append',status='unknown')
      !open(408,file=trim(foldername)//'srfst_fn+00'//jobname,access='append',status='unknown')
 
   open(501,file=trim(foldername)//'slip'//jobname,access='append',status='unknown')
   open(502,file=trim(foldername)//'shear_stress'//jobname,access='append',status='unknown')
   open(503,file=trim(foldername)//'normal_stress'//jobname,access='append',status='unknown')

      do i=1,nmv
         write(30,130)tmv(i),dlog10(maxv(i)*1d-3/yrs)
        do j=311,322
         write(j,110) tmv(i),outs1(i,1,j-310),outs1(i,2,j-310),outs1(i,3,j-310),outs1(i,4,j-310), &
           outs1(i,5,j-310)
        end do

        do j=401,407
         write(j,110) tmv(i),obvs(i,1,j-400),obvs(i,2,j-400),obvs(i,3,j-400),obvs(i,4,j-400)
        end do
        write(501,144) tmv(i),dlog10(maxv(i)*1d-3/yrs),obvdp(i,1,:)
        write(502,144) tmv(i),dlog10(maxv(i)*1d-3/yrs),obvdp(i,2,:)
        write(503,144) tmv(i),dlog10(maxv(i)*1d-3/yrs),obvdp(i,3,:)

       end do
      close(30)
      
      do j=311,322
         close(j)
      end do
      
      do j=401,407
        close(j)
      end do
      close(501)
      close(502)
      close(503) 
 
      imv = 0
   end if

    if(ias==nas)then
       open(31,form='unformatted',file=trim(foldername)//'slipz1-inter'//jobname,access='append',status='unknown')
       open(34,file=trim(foldername)//'t-inter'//jobname,access='append',status='unknown')

       do j=1,nas
          do i=1,Nt_all
             write(31) slipz1_inter(i,j)
          end do
       end do
       do i=1,nas
          write(34,*) tas(i)
       end do
       close(31)
       close(34)
       ias = 0 
end if

    if(icos==ncos)then
       open(42,form='unformatted',file=trim(foldername)//'slipz1-v'//jobname,access='append',status='unknown')
       open(45,form='unformatted',file=trim(foldername)//'slipz1-cos'//jobname,access='append',status='unknown')
       open(48,file=trim(foldername)//'t-cos'//jobname,access='append',status='unknown')

       do j=1,ncos
          do i=1,Nt_all
             write(42) slipz1_v(i,j)
             write(45) slipz1_cos(i,j)
          end do
       end do

       do i=1,ncos
          write(48,*) tcos(i)
       end do
       close(42)
       close(45)
       close(48)
       icos = 0 
end if


	if(inul == nnul)then
	open(52,file=trim(foldername)//'vs-nul'//jobname,access='append',status='unknown')
        open(53,file=trim(foldername)//'nul-time'//jobname,access='append',status='unknown')
        do j=1,nnul
           do kk=1,Nt_all
              !              write(52,160)v_nul(kk,j),slip_nul(kk,j)
           end do
        end do
        do i=52,53
           close(i)
        end do
        inul = 0 
	end if
 

   if(isse==nsse)then
      open(25,form='unformatted',file=trim(foldername)//'slipz1_sse'//jobname,access='append',status='unknown')
      open(26,form='unformatted',file=trim(foldername)//'slipz1_tau'//jobname,access='append',status='unknown')
      open(28,file=trim(foldername)//'t_sse'//jobname, access='append', status='unknown')
      do j = 1,nsse 
         do i=1,Nt_all
            write(25) slipz1_sse(i,j)
            write(26) slipz1_tau(i,j)
         end do
      end do

      do i=1,nsse
         write(28,*) tsse(i)
      end do

      isse = 0
      close(25)
      close(26)
      close(28)
  end if


else

   if((imv>0).and.(imv<nmv))then
      open(30,file=trim(foldername)//'maxvall'//jobname,access='append',status='unknown')
      open(311,file=trim(foldername)//'fltst_dp000'//jobname,access='append',status='unknown')
      open(312,file=trim(foldername)//'fltst_dp025'//jobname,access='append',status='unknown')
      open(313,file=trim(foldername)//'fltst_dp050'//jobname,access='append',status='unknown')
      open(314,file=trim(foldername)//'fltst_dp075'//jobname,access='append',status='unknown')
      open(315,file=trim(foldername)//'fltst_dp100'//jobname,access='append',status='unknown')
      open(316,file=trim(foldername)//'fltst_dp125'//jobname,access='append',status='unknown')
      open(317,file=trim(foldername)//'fltst_dp150'//jobname,access='append',status='unknown')
      open(318,file=trim(foldername)//'fltst_dp175'//jobname,access='append',status='unknown')
      open(319,file=trim(foldername)//'fltst_dp200'//jobname,access='append',status='unknown')
      open(320,file=trim(foldername)//'fltst_dp250'//jobname,access='append',status='unknown')
      open(321,file=trim(foldername)//'fltst_dp300'//jobname,access='append',status='unknown')
      open(322,file=trim(foldername)//'fltst_dp350'//jobname,access='append',status='unknown')

     open(401,file=trim(foldername)//'srfst_fn-32'//jobname,access='append',status='unknown')
      open(402,file=trim(foldername)//'srfst_fn-16'//jobname,access='append',status='unknown')
      open(403,file=trim(foldername)//'srfst_fn-08'//jobname,access='append',status='unknown')
      open(404,file=trim(foldername)//'srfst_fn-00'//jobname,access='append',status='unknown')
      open(405,file=trim(foldername)//'srfst_fn+08'//jobname,access='append',status='unknown')
      open(406,file=trim(foldername)//'srfst_fn+16'//jobname,access='append',status='unknown')
      open(407,file=trim(foldername)//'srfst_fn+32'//jobname,access='append',status='unknown')
      !open(408,file=trim(foldername)//'srfst_fn+00'//jobname,access='append',status='unknown')

   open(501,file=trim(foldername)//'slip'//jobname,access='append',status='unknown')
   open(502,file=trim(foldername)//'shear_stress'//jobname,access='append',status='unknown')
   open(503,file=trim(foldername)//'normal_stress'//jobname,access='append',status='unknown')

      do i=1,imv
         write(30,130)tmv(i),dlog10(maxv(i)*1d-3/yrs)
        do j=311,322
         write(j,110) tmv(i),outs1(i,1,j-310),outs1(i,2,j-310),outs1(i,3,j-310),outs1(i,4,j-310), &
           outs1(i,5,j-310)
        end do

       do j=401,407
         write(j,110) tmv(i),obvs(i,1,j-400),obvs(i,2,j-400),obvs(i,3,j-400),obvs(i,4,j-400)
        end do
        write(501,144) tmv(i),dlog10(maxv(i)*1d-3/yrs),obvdp(i,1,:)
        write(502,144) tmv(i),dlog10(maxv(i)*1d-3/yrs),obvdp(i,2,:)
        write(503,144) tmv(i),dlog10(maxv(i)*1d-3/yrs),obvdp(i,3,:)

      end do
       close(30)
       do j=311,322
         close(j)
       end do 
      do j=401,407
        close(j)
      end do
      close(501)
      close(502)
      close(503)

      imv = 0
    end if

    if((ias>0).and.(ias<nas))then
       open(31,form='unformatted',file=trim(foldername)//'slipz1-inter'//jobname,access='append',status='unknown')
       open(34,file=trim(foldername)//'t-inter'//jobname,access='append',status='unknown')

       do j=1,ias
          do i=1,Nt_all
             write(31) slipz1_inter(i,j)
          end do
       end do
       do i=1,ias
          write(34,*) tas(i)
       end do

       close(31)
       close(34)

       ias = 0 
     end if


   if(isse<nsse.and.isse>0)then
      open(25,form='unformatted',file=trim(foldername)//'slipz1_sse'//jobname,access='append',status='unknown')
      open(26,form='unformatted',file=trim(foldername)//'slipz1_tau'//jobname,access='append',status='unknown')
      open(28,file=trim(foldername)//'t_sse'//jobname, access='append',status='unknown')
      do j=1,isse 
         do i=1,Nt_all
            write(25) slipz1_sse(i,j)
            write(26) slipz1_tau(i,j)
         end do
      end do
      do i=1,isse
         write(28,*) tsse(i)
      end do

      isse = 0
      close(25)
      close(26)
      close(28)
  end if


     if((icos>0).and.(icos<ncos))then
       open(45,form='unformatted',file=trim(foldername)//'slipz1-cos'//jobname,access='append',status='unknown')
      open(42,form='unformatted',file=trim(foldername)//'slipz1-v'//jobname,access='append',status='unknown')
       open(48,file=trim(foldername)//'t-cos'//jobname,access='append',status='unknown')
       
       do j=1,icos
          do i=1,Nt_all
             write(45) slipz1_cos(i,j)
             write(42) slipz1_v(i,j)
          end do
       end do
       do i=1,icos
          write(48,*) tcos(i)
       end do
       close(42)     
       close(45)
       close(48)
       icos = 0 
      end if

                 
	if((inul>0).and.(inul<nnul))then
	open(52,file=trim(foldername)//'vs-nul'//jobname,access='append',status='unknown')
        open(53,file=trim(foldername)//'nul-time'//jobname, &
             access='append',status='unknown')
		do j=1,inul
              do kk=1,Nt_all
     !            write(52,160)v_nul(kk,j),slip_nul(kk,j)
            end do 
      !      write(53,140)tnul(j)
		end do 
		do i=52,53
			close(i)
		end do 
		inul = 0 
end if

end if 

 110    format(E22.15,7(1X,E15.7))
 120    format(E20.13,4X,E20.13,4X,I6)
 130    format(E22.15,2(1X,E15.7))
 140    format(E20.13)
 150    format(E22.13,3(1X,E15.7))
 160    format(E20.13,1x,E20.13)
 500    format(E15.8,1X,E20.13)
 600    format(E15.4,1X,E13.6,1X,E15.8)
 700    format(E13.6)
 900    format(E15.8)
144 format(E22.14,166(1X,E15.7))

RETURN
END subroutine output


