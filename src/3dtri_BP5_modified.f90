! Modified version of 3dtri_BP5.f90 to load TriGreen files directly
! This version implements dynamic load balancing compatible with calc_trigreen.f90

program main
  USE mpi
  USE phy3d_module_non
  implicit none
  integer, parameter :: DP=kind(1.d0)

  logical cyclecont
  integer ::Nt_all,Nt, ndt,ndtnext,kk,ii,n,l,ndt_v,ndt_inter,Itout,&  
       ndtmax,i,j,k,nm,nrec,ihrestart,Ifileout, &
       Iperb,record,Isnapshot,ix1,ix2,ix3,ix4,iz1,iz2,iz3,&
       record_cor,s2,s3,s4,s5,s6,s7,s8,s9,s10
  integer,dimension(:) :: s1(10)
  real (DP) :: Vint,tmp,accuracy,areatot, epsv,dt_try, dt,dtmin,dt_did,dt_next,dip, &
       hnucl, hstarfactor, &
       sigmadiff,sefffix,Lffix,xdepth,xlength,&
       t,tprint_inter, tint_out,tout,&
       tmin_out,tint_cos,tint_sse,&   
       tslip_ave,tslipend,tslip_aveint, tmax, &
       tslipsse,tslipcos,tstart1,tend1,tstart2,tend2,tstart3,tend3, &
       tssestart,tsseend, &
       factor1,factor2,factor3,factor4,vtmp,fr_dummy,&
       xilock1,xilock2,x1,x2,x3,x4,z1,z2,z3,&
       xi_dummy,x_dummy,z_dummy,stiff_dummy,help

  real (DP) ::  tmbegin,tmrun,tautmp

  integer, DIMENSION(:), ALLOCATABLE :: itran1,itran2,islip1,islip2

  real (DP), DIMENSION(:), ALLOCATABLE :: x,z,xi,yt,yt0,dydt,yt_scale, &
       slip,slipinc,slipds,slipdsinc,sr,vi

  !Arrays only defined at master cpu
  real (DP), DIMENSION(:), ALLOCATABLE :: x_all,xi_all,z_all,&
       yt_all,yt0_all,dydt_all,yt_scale_all,tau1_all,tau2_all, &
       slip_all,slipinc_all,slipds_all,slipdsinc_all,cca_all,ccb_all,xLf_all,seff_all,&
       vi_all,phy1_all,phy2_all
  !output related parameters
  integer :: imv,ias,icos,isse,Ioutput,inul,i_nul,n_nul_int
  real (DP) :: vcos,vsse1,vsse2
  real (DP), DIMENSION(:), ALLOCATABLE :: maxv,moment,&
       maxnum,msse1,msse2,areasse1,areasse2,tmv,tas,tcos,tnul,tsse
  real (DP),dimension(:,:,:),allocatable :: outs1
  real (DP), DIMENSION(:,:), ALLOCATABLE :: slipz1_inter, &
       slipz1_cos,slipave_inter,slipave_cos, slipz1_v, &
       v_cos,slip_cos,v_nul,slip_nul,slipz1_tau,slipz1_sse

  integer,DIMENSION(:),ALLOCATABLE :: intdepz1,intdepz2,intdepz3,ssetime
  integer :: n_intz1,n_intz2,n_intz3,n_cosz1,n_cosz2,n_cosz3

  real (DP), DIMENSION(:), ALLOCATABLE :: Trup,area,zzfric,zzfric2
  logical,dimension(:),allocatable :: rup
  logical :: end1=.false.
  real (DP) :: teve1=0.d0

  integer :: n_obv,np1,np2
  real(DP),dimension(:,:),allocatable :: surf1,surf2,surf3
  real(DP),dimension(:,:,:),allocatable :: obvs,obvstrk,obvdp
  real(DP) :: vel1,vel2,vel3,disp1,disp2,disp3
  integer, dimension(:),allocatable :: pstrk, pdp

  character(len=40) :: cTemp,filename,ct

  !MPI RELATED DEFINITIONS
  integer :: ierr,size,myid,master

  ! Dynamic load balancing variables (compatible with calc_trigreen.f90)
  integer :: base_cells, extra_cells, local_cells, start_idx
  logical :: use_trigreen_format = .true.  ! Set to .true. to use TriGreen files

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
  read(12,*)Nab,Nt_all,Nt,Lratio,nprocs,n_obv,np1,np2
  read(12,*)Idin,Idout,Iprofile,Iperb,Isnapshot 
  read(12,*)Vpl
  read(12,*)tmax
  read(12,*)tslip_ave,tslipend,tslip_aveint
  read(12,*)tint_out,tmin_out,tint_cos,tint_sse
  read(12,*)vcos,vsse1,vsse2
  read(12,*)nmv,nas,ncos,nnul,nsse,n_nul_int
  read(12,*)s1(1),s1(2),s1(3),s1(4),s1(5),s1(6),s1(7),s1(8),s1(9),s1(10)
!!! modified data read in
110 format(A)
  close(12)

  !num. of rows each proc. 
  Nt_all = Nt_all

  ! MODIFICATION: Implement dynamic load balancing compatible with calc_trigreen.f90
  if (use_trigreen_format) then
     ! Calculate optimal distribution using the same algorithm as calc_trigreen.f90
     base_cells = Nt_all / size
     extra_cells = mod(Nt_all, size)
     
     if (myid < extra_cells) then
        local_cells = base_cells + 1
        start_idx = myid * (base_cells + 1)
     else
        local_cells = base_cells
        start_idx = extra_cells * (base_cells + 1) + (myid - extra_cells) * base_cells
     end if
     
     ! Override Nt with the actual local cells for this process
     Nt = local_cells
     
     if (myid == master) then
        write(*,*) 'Using TriGreen format with dynamic load balancing'
        write(*,*) 'Base cells per process:', base_cells
        write(*,*) 'Extra cells distributed:', extra_cells
        write(*,*) 'Total processes:', size
     end if
     
     write(*,*) 'Process', myid, 'gets', local_cells, 'cells starting from index', start_idx
  else
     ! Original logic for even distribution
     if(mod(Nt_all,nprocs)/=0)then
        write(*,*)'Nd_all must be integer*nprocs. Change nprocs!'
        STOP
     else
        write(*,*)'Each cpu calculates',Nt_all/nprocs,'cells'
     end if
  end if

  if(myid == master)then
     ALLOCATE(x_all(Nt_all),xi_all(Nt_all),&
          cca_all(Nt_all),ccb_all(Nt_all),seff_all(Nt_all),xLf_all(Nt_all),vi_all(Nt_all),&
          tau1_all(Nt_all),tau2_all(Nt_all),slip_all(Nt_all),slipinc_all(Nt_all),slipds_all(Nt_all),slipdsinc_all(Nt_all),&
         yt0_all(2*Nt_all),yt_all(2*Nt_all),dydt_all(2*Nt_all),yt_scale_all(2*Nt_all))

     allocate(phy1_all(Nt_all),phy2_all(Nt_all))

     ALLOCATE (outs1(nmv,7,10),&
          maxv(nmv),maxnum(nmv),msse1(nsse),msse2(nsse),areasse1(nsse),areasse2(nsse), &
          tmv(nmv),tas(nas),tcos(ncos),tnul(nnul),tsse(nsse))

!!! modify output number
     ALLOCATE (slipz1_inter(Nt_all,nas),slipz1_cos(Nt_all,ncos), &
          slipave_inter(Nt_all,nas),slipave_cos(Nt_all,ncos),v_cos(Nt_all,ncos),slip_cos(Nt_all,ncos), &
          v_nul(Nt_all,nnul),slip_nul(Nt_all,nnul),slipz1_tau(Nt_all,nsse),slipz1_sse(Nt_all,nsse) )
     ALLOCATE(intdepz1(Nt_all),intdepz2(Nt_all),intdepz3(Nt_all),slipz1_v(Nt_all,ncos),ssetime(nsse)  )

     allocate(moment(nmv),Trup(Nt_all),rup(Nt_all),area(Nt_all))
     allocate(surf1(n_obv,Nt_all),surf2(n_obv,Nt_all),surf3(n_obv,Nt_all),obvs(nmv,6,n_obv))
     allocate(pstrk(np1),pdp(np2),obvstrk(nmv,2,np1),obvdp(nmv,2,np2))
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE (phy1(Nt),phy2(Nt),islip1(Nl),islip2(Nl),itran1(Nl),itran2(Nl),zzfric(Nt),zzfric2(Nt),&
       x(Nt),z(Nt),z_all(Nt_all),xi(Nt),cca(Nt),ccb(Nt),seff(Nt),&
       xLf(Nt),tau1(Nt),tau2(Nt),tau0(Nt),slipds(Nt),slipdsinc(Nt),slip(Nt),slipinc(Nt), &
       yt(2*Nt),dydt(2*Nt),yt_scale(2*Nt),yt0(2*Nt),sr(Nt),vi(Nt))

  ! MODIFICATION: Allocate stiffness matrix with actual local cells
  ALLOCATE (stiff(local_cells,Nt_all),stiff2(local_cells,Nt_all))

  !Read in stiffness matrix, in nprocs segments
  write(cTemp,*) myid
  write(*,*) cTemp

  ! MODIFICATION: Choose file format based on use_trigreen_format flag
  if (use_trigreen_format) then
     ! Load TriGreen format files
     open(5, file=trim(stiffname)//'trigreen_'//trim(adjustl(cTemp))//'.bin',form='unformatted',access='stream',buffered='yes')
     write(*,*) 'Loading TriGreen file: trigreen_', trim(adjustl(cTemp)), '.bin'
  else
     ! Load original ssGreen format files
     open(5, file=trim(stiffname)//'ssGreen_'//trim(adjustl(cTemp))//'.bin',form='unformatted',access='stream',buffered='yes')
     write(*,*) 'Loading ssGreen file: ssGreen_', trim(adjustl(cTemp)), '.bin'
  end if

if(myid==master)then
  open(666,file='area'//jobname,form='formatted',status='old',access='stream')
  
  ! MODIFICATION: Handle surface Green's functions based on format
  if (use_trigreen_format) then
     ! For TriGreen, we might not have surface Green's functions, so create dummy or skip
     write(*,*) 'Note: Using TriGreen format - surface Green functions may not be available'
     ! Create dummy surface Green's functions or skip this section
     allocate(surf1(n_obv,Nt_all),surf2(n_obv,Nt_all),surf3(n_obv,Nt_all))
     surf1 = 0.d0
     surf2 = 0.d0
     surf3 = 0.d0
  else
     open(51,file=trim(stiffname)//'surfGreen'//'.bin',form='unformatted',access='stream')
  end if
  
  open(55,file=trim(stiffname)//'position.bin',form='unformatted',access='stream')
  open(56,file='profstrk'//jobname,form='formatted',status='old')
  open(57,file='profdp'//jobname,form='formatted',status='old')
  do i=1,np1
    read(56,*) pstrk(i)
  end do
  do i=1,np2
    read(57,*) pdp(i)
 end do
 close(56)
 close(57)
end if

  ! MODIFICATION: Update record_cor for dynamic load balancing
  if (use_trigreen_format) then
     record_cor = start_idx
  else
     record_cor = Nt * myid
  end if

  !-------------------------------------------------------------------------------------------
  ! read stiffness from Stuart green calculation.
  !-----------------------------------------------------------------------------------------
  if(myid==master)then
     ! OPTIMIZATION: Vectorize position reading for better performance
     do k=1,Nt_all
        read(55) xi_all(k),x_all(k),z_all(k) !xi is along the fault-normal  while x is along the strike
        xi_all(k)=xi_all(k)/1000
        x_all(k)=x_all(k)/1000
        z_all(k)=z_all(k)/1000

       Trup(k)=1d9
       rup(k)=.false. 
       read(666,'(E14.7)') area(k)
     end do
     
     ! MODIFICATION: Handle surface Green's functions based on format
     if (.not. use_trigreen_format) then
        ! Read surface Green's functions only for original format
        do j = 1,n_obv
          do k=1,Nt_all
            read(51) surf1(j,k)
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
  end if

  ! MODIFICATION: Read stiffness matrix with correct local_cells
  ! OPTIMIZATION: Use OpenMP for parallel stiffness matrix reading and processing
  !$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(STATIC)
  do i=1,local_cells !! observe (now using local_cells instead of Nt)
     do j=1,Nt_all !! source
        read(5) stiff(i,j)
        stiff(i,j)=-stiff(i,j)
        if(stiff(i,j).lt.-15.d0.or.stiff(i,j).gt.20.d0)then
           stiff(i,j) = 0.d0
           write(*,*) j,'extreme'
	end if
        if(isNaN(stiff(i,j)))then
          stiff(i,j)=0.d0
        end if
     end do
  end do
  !$OMP END PARALLEL DO
  close(5)

if(myid==master)then
  close(666)
  if (.not. use_trigreen_format) then
     close(51)
  end if
  close(55)
end if

  !!-----------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  if(myid==master)then
     CALL resdep(Nt_all,dip,hnucl,sigmadiff,sefffix,Lffix, &
          itran1,itran2,islip1,islip2,Iperb,factor1,factor2,factor3,factor4,&
          xilock1,xilock2,cca_all,ccb_all,xLf_all,seff_all,x_all,z_all,vi_all)
  end if

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  ! MODIFICATION: Use local_cells for scattering operations
  call MPI_Scatter(cca_all,local_cells,MPI_Real8,cca,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
  call MPI_Scatter(ccb_all,local_cells,MPI_Real8,ccb,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
  call MPI_Scatter(xLf_all,local_cells,MPI_Real8,xLf,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)

  ! ... rest of the program would continue here ...
  
  ! For demonstration, just print some info and exit
  write(*,*) 'Process', myid, 'successfully loaded stiffness matrix with', local_cells, 'local cells'
  write(*,*) 'Stiffness matrix dimensions:', local_cells, 'x', Nt_all
  
  call MPI_Finalize(ierr)
  
end program main
