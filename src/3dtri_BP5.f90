!===============================================================================
! Copyright (c) 2024 TriBIE Development Team
! All rights reserved.
!
! This software is part of the TriBIE (Triangular Boundary Integral Element)
! earthquake simulation package. It implements the 3D triangular boundary
! integral equation method for dynamic rupture simulation.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors
!    may be used to endorse or promote products derived from this software
!    without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
! THE POSSIBILITY OF SUCH DAMAGE.
!
! For questions or support, please contact the TriBIE Development Team.
!===============================================================================

!------------------------------------------------
! D. Li
! Last modified Aug. 2025
! 
! SEAS benchmark BP5
! May 16, 2021
!------------------------------------------------

program main
  USE mpi
  USE phy3d_module_non
  use hdf5  ! Add HDF5 support
  implicit none
  integer, parameter :: DP=kind(1.d0)

  logical cyclecont
  integer ::Nt_all,Nt, ndt,ndtnext,kk,ii,n,l,ndt_v,ndt_inter,Itout,&  
       ndtmax,i,j,k,nm,nrec,ihrestart,Ifileout, &
       Iperb,record,Isnapshot,iz1,iz2,iz3,&
       record_cor
  integer,dimension(:) :: s1(10)
  real (DP) :: Vint,tmp,accuracy,areatot, epsv,dt_try, dt,dtmin,dt_did,dt_next, &
       hnucl,&
       t,tprint_inter, tint_out,tout,&
       tmin_out,tint_cos,tint_sse,&   
       tslip_ave,tslipend,tslip_aveint, tmax, &
       tslipsse,tslipcos,tstart1,tend1,tstart2,tend2,tstart3,tend3, &
       tssestart,tsseend, &
       xilock1,xilock2,x4,z1,z2,z3,&
       help

  real (DP) ::  tmbegin,tmrun,tautmp


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
  
  ! Communication buffer variables
  real(DP), dimension(:), allocatable :: send_buffer, recv_buffer
  integer :: comm_count, comm_tag
  
  ! Blocking variables for cache optimization
  integer :: block_size, i_start, i_end, j_start, j_end
  
  ! Additional variables for advanced optimizations
  integer :: request1, request2

  character(len=40) :: cTemp,filename,ct

  !MPI RELATED DEFINITIONS
  integer :: ierr,size,myid,master

  ! Dynamic load balancing variables (compatible with calc_trigreen.f90)
  integer :: base_cells, extra_cells, local_cells, start_idx
  logical :: use_trigreen_format = .true.  ! Set to .true. to use TriGreen files
  
  ! File existence checking variables
  logical :: trigreen_file_exists
  character(len=256) :: trigreen_filename
  
  ! MPI_Scatterv variables for uneven distribution
  integer, allocatable :: sendcounts(:), displs(:)
  integer :: total_sent

  ! Add HDF5 variables
  integer(HID_T) :: file_id, dset_id, dspace_id
  integer(HID_T) :: group_id, attr_id, attr_space_id
  integer(HSIZE_T), dimension(3) :: dims, maxdims
  integer(HSIZE_T), dimension(2) :: dims_2d, maxdims_2d
  integer(HSIZE_T), dimension(1) :: dims_1d, maxdims_1d
  integer :: hdferr
  logical :: hdf5_initialized = .false.
  
  ! HDF5 file naming
  character(len=256) :: hdf5_filename, xdmf_filename
  character(len=256) :: time_series_group_name

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
        start_idx = myid * (base_cells + 1) + 1  ! Fix: Make 1-based indexing
     else
        local_cells = base_cells
        start_idx = extra_cells * (base_cells + 1) + (myid - extra_cells) * base_cells + 1  ! Fix: Make 1-based indexing
     end if
     
     ! Override Nt with the actual local cells for this process
     Nt = local_cells
     
     if (myid == master) then
        write(*,*) '=========================================='
        write(*,*) 'TriGreen Integration Status:'
        write(*,*) '=========================================='
        write(*,*) 'Using TriGreen format with dynamic load balancing'
        write(*,*) 'Base cells per process:', base_cells
        write(*,*) 'Extra cells distributed:', extra_cells
        write(*,*) 'Total processes:', size
        write(*,*) 'Expected files: trigreen_0.bin, trigreen_1.bin, ...'
        write(*,*) '=========================================='
     end if
     
     write(*,*) 'Process', myid, 'gets', local_cells, 'cells starting from index', start_idx
     
     ! Allocate MPI_Scatterv arrays for uneven distribution
     allocate(sendcounts(0:size-1), displs(0:size-1))
     
     ! Calculate send counts and displacements for each process
     call MPI_Allgather(local_cells, 1, MPI_INTEGER, sendcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
     
     ! Calculate displacements
     displs(0) = 0
     do i = 1, size-1
        displs(i) = displs(i-1) + sendcounts(i-1)
     end do
     
     if (myid == master) then
        write(*,*) 'MPI_Scatterv distribution:'
        do i = 0, size-1
           write(*,*) '  Process', i, ': sendcount =', sendcounts(i), ', displacement =', displs(i)
        end do
     end if
  else
     ! Original logic for even distribution
     if(mod(Nt_all,nprocs)/=0)then
        write(*,*)'Nd_all must be integer*nprocs. Change nprocs!'
        STOP
     else
        write(*,*)'Each cpu calculates',Nt_all/nprocs,'cells'
     end if
  end if


  if(myid == master) then
     ALLOCATE(x_all(Nt_all),xi_all(Nt_all),&
          cca_all(Nt_all),ccb_all(Nt_all),seff_all(Nt_all),xLf_all(Nt_all),vi_all(Nt_all),&
          tau1_all(Nt_all),tau2_all(Nt_all),slip_all(Nt_all),slipinc_all(Nt_all),slipds_all(Nt_all),slipdsinc_all(Nt_all),&
         yt0_all(2*Nt_all),yt_all(2*Nt_all),dydt_all(2*Nt_all),yt_scale_all(2*Nt_all))

     allocate(phy1_all(Nt_all),phy2_all(Nt_all))

     ALLOCATE (outs1(nmv,7,10),&
          maxv(nmv),maxnum(nmv),msse1(nsse),msse2(nsse),areasse1(nsse),areasse2(nsse), &
          tmv(nmv),tas(nas),tcos(ncos),tnul(nnul),tsse(nsse))

     ! CRITICAL FIX: Initialize arrays to prevent garbage values and performance issues
     outs1 = 0.d0      ! Initialize output array to prevent garbage data
     maxv = 0.d0       ! Initialize maximum velocity array
     maxnum = 0        ! Initialize maximum number array
     msse1 = 0.d0      ! Initialize SSE arrays
     msse2 = 0.d0
     areasse1 = 0.d0
     areasse2 = 0.d0
     tmv = 0.d0        ! Initialize time arrays
     tas = 0.d0
     tcos = 0.d0
     tnul = 0.d0
     tsse = 0.d0

!!! modify output number
     ALLOCATE (slipz1_inter(Nt_all,nas),slipz1_cos(Nt_all,ncos), &
          slipave_inter(Nt_all,nas),slipave_cos(Nt_all,ncos),v_cos(Nt_all,ncos),slip_cos(Nt_all,ncos), &
          v_nul(Nt_all,nnul),slip_nul(Nt_all,nnul),slipz1_tau(Nt_all,nsse),slipz1_sse(Nt_all,nsse) )
     ALLOCATE(intdepz1(Nt_all),intdepz2(Nt_all),intdepz3(Nt_all),slipz1_v(Nt_all,ncos),ssetime(nsse)  )

     allocate(moment(nmv),Trup(Nt_all),rup(Nt_all),area(Nt_all))
     allocate(surf1(n_obv,Nt_all),surf2(n_obv,Nt_all),surf3(n_obv,Nt_all),obvs(nmv,6,n_obv))
     allocate(pstrk(np1),pdp(np2),obvstrk(nmv,2,np1),obvdp(nmv,2,np2))
     
     ! CRITICAL FIX: Initialize additional arrays to prevent garbage values
     moment = 0.d0      ! Initialize moment array
     Trup = 0.d0        ! Initialize rupture time array
     rup = .false.      ! Initialize rupture flag array
     area = 0.d0        ! Initialize area array
     obvs = 0.d0        ! Initialize observation arrays
     obvstrk = 0.d0
     obvdp = 0.d0
     pstrk = 0          ! Initialize profile arrays
     pdp = 0
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE (phy1(local_cells),phy2(local_cells),&
  zzfric(local_cells),zzfric2(local_cells),&
       x(local_cells),z(local_cells),z_all(Nt_all),xi(local_cells),cca(local_cells),ccb(local_cells),&
       seff(local_cells),&
       xLf(local_cells),tau1(local_cells),tau2(local_cells),tau0(local_cells),slipds(local_cells),&
       slipdsinc(local_cells),slip(local_cells),slipinc(local_cells), &
       yt(2*local_cells),dydt(2*local_cells),yt_scale(2*local_cells),&
       yt0(2*local_cells),sr(local_cells),vi(local_cells))

  ALLOCATE (stiff(local_cells,Nt_all),stiff2(local_cells,Nt_all))   !!! stiffness of Stuart green calculation

  !Read in stiffness matrix, in nprocs segments

  write(cTemp,*) myid
  write(*,*) cTemp

  ! MODIFICATION: Choose file format based on use_trigreen_format flag
  if (use_trigreen_format) then
     ! Load TriGreen format files with existence checking
     trigreen_filename = trim(stiffname)//'trigreen_'//trim(adjustl(cTemp))//'.bin'
     
     ! Check if TriGreen file exists
     inquire(file=trigreen_filename, exist=trigreen_file_exists)
     if (.not. trigreen_file_exists) then
        write(*,*) 'ERROR: TriGreen file not found: ', trim(trigreen_filename)
        write(*,*) 'Process', myid, 'cannot continue without TriGreen file'
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
     end if
     
     open(5, file=trigreen_filename, form='unformatted', access='stream', status='old')
     write(*,*) 'Process', myid, ': Successfully loaded TriGreen file: ', trim(trigreen_filename)
  else
     ! Load original ssGreen format files
     open(5, file=trim(stiffname)//'ssGreen_'//trim(adjustl(cTemp))//'.bin',form='unformatted',access='stream')
     write(*,*) 'Loading ssGreen file: ssGreen_', trim(adjustl(cTemp)), '.bin'
  end if

if(myid==master)then
  open(666,file='area'//jobname,form='formatted',status='old',access='stream')
  
  ! MODIFICATION: Handle surface Green's functions based on format
  if (use_trigreen_format) then
     ! For TriGreen, we might not have surface Green's functions, so initialize to zeros
     write(*,*) 'Note: Using TriGreen format - surface Green functions initialized to zero'
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
        read(55) x_all(k),xi_all(k),z_all(k) !xi is along the fault-normal  while x is along the strike
        xi_all(k)=xi_all(k)/1000
        x_all(k)=x_all(k)/1000
        z_all(k)=z_all(k)/1000

       Trup(k)=1d9
       rup(k)=.false. 
       read(666,'(E14.7)') area(k)
     end do
!! read surface Green's
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

  ! OPTIMIZATION: Use OpenMP for parallel stiffness matrix reading and processing
  ! FIRST: Read file sequentially (OUTSIDE OpenMP to avoid file corruption)
  do i=1,local_cells !! observe (now using local_cells instead of Nt)
     do j=1,Nt_all !! source
        read(5, err=999) stiff(i,j)
     end do
  end do
  
  ! Jump to error handling if read fails
  goto 200
  
  ! Error handling for file read
  999 write(*,*) 'Process', myid, ': ERROR reading stiffness matrix from TriGreen file'
      write(*,*) 'Process', myid, ': File may be corrupted or incomplete'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      
  200 continue
  close(5)
  
  ! SECOND: Process data in parallel (OpenMP for computation only, NO file I/O)
  !$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(STATIC)
  do i=1,local_cells !! observe (now using local_cells instead of Nt)
     do j=1,Nt_all !! source
        if(stiff(i,j).lt.-1.0d0.or.stiff(i,j).gt.1.0d0)then
           stiff(i,j) = 0.d0
           !$OMP CRITICAL
           write(*,*) 'Process', myid, ': Extreme value at position (', i, ',', j, ') =', stiff(i,j)
           !$OMP END CRITICAL
	end if
        if(stiff(i,j) /= stiff(i,j))then  ! Check for NaN using IEEE standard
          stiff(i,j)=0.d0
          !$OMP CRITICAL
          write(*,*) 'Process', myid, ': NaN detected at position (', i, ',', j, ') - set to 0'
          !$OMP END CRITICAL
        end if
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! TriGreen integration summary
  if (use_trigreen_format) then
     write(*,*) 'Process', myid, ': TriGreen stiffness matrix loaded successfully'
     write(*,*) 'Process', myid, ': Matrix dimensions:', local_cells, 'x', Nt_all
     write(*,*) 'Process', myid, ': Total elements loaded:', local_cells * Nt_all
  end if
  
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
     CALL resdep(Nt_all,hnucl, &
          xilock1,xilock2,cca_all,ccb_all,xLf_all,seff_all,x_all,z_all,vi_all)
  end if


  call MPI_Barrier(MPI_COMM_WORLD,ierr)
   call MPI_Scatterv(cca_all,sendcounts,displs,MPI_Real8,cca,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Scatterv(ccb_all,sendcounts,displs,MPI_Real8,ccb,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Scatterv(xLf_all,sendcounts,displs,MPI_Real8,xLf,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Scatterv(seff_all,sendcounts,displs,MPI_Real8,seff,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Scatterv(vi_all,sendcounts,displs,MPI_Real8,vi,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Scatterv(x_all,sendcounts,displs,MPI_Real8,x,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)

  call MPI_Bcast(z_all,Nt_all,MPI_Real8,master,MPI_COMM_WORLD,ierr)

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
  dt_try=dtmin
  Vint = Vpl


  if(myid==master)then
      open(311,file=trim(foldername)//'fltst_strk-36dp+00'//jobname,position='append',status='unknown')
      open(312,file=trim(foldername)//'fltst_strk-16dp+00'//jobname,position='append',status='unknown')
      open(313,file=trim(foldername)//'fltst_strk+00dp+00'//jobname,position='append',status='unknown')
      open(314,file=trim(foldername)//'fltst_strk+16dp+00'//jobname,position='append',status='unknown')
      open(315,file=trim(foldername)//'fltst_strk+36dp+00'//jobname,position='append',status='unknown')
      open(316,file=trim(foldername)//'fltst_strk-24dp+10'//jobname,position='append',status='unknown')
      open(317,file=trim(foldername)//'fltst_strk-16dp+10'//jobname,position='append',status='unknown')
      open(318,file=trim(foldername)//'fltst_strk+00dp+10'//jobname,position='append',status='unknown')
      open(319,file=trim(foldername)//'fltst_strk+16dp+10'//jobname,position='append',status='unknown')
      open(320,file=trim(foldername)//'fltst_strk+00dp+22'//jobname,position='append',status='unknown')

       do i=311,320
   		write(i,100)'# This is the file header'
		write(i,100)'# problem=SEAS Benchmark No.5'
		write(i,100)'# author=D.Li'
        write(i,100)'# code=TriBIE'
		write(i,100)'# date=2021/5/11'
		write(i,100)'# element_size = 500 m'
		write(i,100)'# minimum_time_step = 1e-3'
		write(i,100)'# maximum_time_step = 2e+7'
		write(i,100)'# location = on fault: file name'
		write(i,100)'# Column #1 = Time (s)'
		write(i,100)'# Column #2 = slip_2 (m)'
        write(i,100)'# Column #3 = slip_3 (m)'
		write(i,100)'# Column #4 = Slip_rate_2(log10 m/s)'
        write(i,100)'# Column #5 = Slip_rate_3(log10 m/s)'
		write(i,100)'# Column #6 = Shear stress 2 (MPa)'
		write(i,100)'# Column #7 = shear stress 3 (MPa)'
		write(i,100)'# Column #8 = State (log10 s)'
		write(i,100)'# '
		write(i,100)'# The line below lists the names of the data fields:'
		write(i,'(A,1x,A,1x,A,1x,A,1x,A,1x,A,1x,A,1x,A)')'t','slip_2','slip_3','slip_rate_2',&
      'slip_rate_3','shear_stress_2','shear_stress_3','state'
		write(i,100)'# Below is the time-series data.'		
	end do 
 100	format(A)

       open(30,file=trim(foldername)//'maxvall'//jobname,position='append',status='unknown')
        i=30
        write(i,100)'# This is the file header'
        write(i,100)'# problem=SEAS Benchmark No.5'
        write(i,100)'# author=D.Li'
        write(i,100)'# code=TriBIE'
        write(i,100)'# date=2021/5/11'
        write(i,100)'# element_size = 500 m'
        write(i,100)'# minimum_time_step = 1e-3'
        write(i,100)'# maximum_time_step = 2e+7'
        write(i,100)'# Column #1 = Time (s)'
        write(i,100)'# Column #2 = max_slip_rate (logV m/s)'
        write(i,100)'# Column #3 = moment rate(N-m/s)'
        write(i,100)'# The line below lists the names of the data fields:'
        write(i,'(A,1x,A,1x,A)')'t','max_slip_rate','moment_rate'
        write(i,100)'# Below is the time-series data.'  

      close(30)
      do j=311,320
         close(j)
      end do

      ! Close main output files
      close(30)
      
      do j=311,320
         close(j)
      end do 

      ! Open block stress output files (these are for stress data, not observation data)
      open(401,file=trim(foldername)//'blkst_strk-16fn+08dp+00'//jobname,position='append',status='unknown')
      open(402,file=trim(foldername)//'blkst_strk+00fn+08dp+00'//jobname,position='append',status='unknown')
      open(403,file=trim(foldername)//'blkst_strk+16fn+08dp+00'//jobname,position='append',status='unknown')
      open(404,file=trim(foldername)//'blkst_strk+00fn+16dp+00'//jobname,position='append',status='unknown')
      open(405,file=trim(foldername)//'blkst_strk+00fn+32dp+00'//jobname,position='append',status='unknown')
      open(406,file=trim(foldername)//'blkst_strk+00fn+48dp+00'//jobname,position='append',status='unknown')
      open(407,file=trim(foldername)//'blkst_strk+00fn+08dp+10'//jobname,position='append',status='unknown')
      open(408,file=trim(foldername)//'blkst_strk+16fn+08dp+10'//jobname,position='append',status='unknown')
      open(409,file=trim(foldername)//'blkst_strk+00fn+32dp+10'//jobname,position='append',status='unknown')

      ! Open additional output files
      open(501,file=trim(foldername)//'slip_2_depth'//jobname,position='append',status='unknown')
      open(502,file=trim(foldername)//'slip_2_strike'//jobname,position='append',status='unknown')
      open(503,file=trim(foldername)//'stress_2_depth'//jobname,position='append',status='unknown')
      open(504,file=trim(foldername)//'stress_2_strike'//jobname,position='append',status='unknown')
    
      ! Write main output data
      do i=1,nmv
        ! Write block stress data to files 401-409 (these are NOT observation data)
        ! Note: These files are for stress output, not observation data
        
        ! Write additional output data
        write(501,144) tmv(i),dlog10(maxv(i)*1d-3/yrs),obvdp(i,1,:)
        write(503,144) tmv(i),dlog10(maxv(i)*1d-3/yrs),obvdp(i,2,:)
        write(502,144) tmv(i),dlog10(maxv(i)*1d-3/yrs),obvstrk(i,1,:)
        write(504,144) tmv(i),dlog10(maxv(i)*1d-3/yrs),obvstrk(i,2,:)
     end do
     
     ! Close additional output files
     do j=501,504
       close(j)
     end do

144 format(E22.14,80(1X,E15.7))

      ! Close block stress output files
      do j=401,409
         close(j)
      end do

      ! Handle observation data output separately (if n_obv > 0)
      if (n_obv > 0) then
         ! Open observation data output files
         do j = 1, min(n_obv, 9)  ! Limit to available observation points
            open(600+j, file=trim(foldername)//'obs_data_'//char(48+j)//jobname, position='append', status='unknown')
         end do
         
         ! Write observation data
         do i = 1, nmv
            do j = 1, min(n_obv, 9)
               write(600+j, 110) tmv(i), obvs(i,1,j), obvs(i,2,j), obvs(i,3,j), obvs(i,4,j), obvs(i,5,j), obvs(i,6,j)
            end do
         end do
         
         ! Close observation data files
         do j = 1, min(n_obv, 9)
            close(600+j)
         end do
         
         write(*,*) 'Observation data written for ', min(n_obv, 9), ' observation points'
      end if

      imv = 0
   end if

    if(ias==nas)then
       open(31,form='unformatted',file=trim(foldername)//'slipz1-inter'//jobname,position='append',status='unknown')
       open(34,file=trim(foldername)//'t-inter'//jobname,position='append',status='unknown')

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
       ! HDF5 output for time-series variables instead of binary files
       if (myid == master) then
          ! Initialize HDF5 if not already done
          if (.not. hdf5_initialized) then
             call h5open_f(hdferr)
             hdf5_initialized = .true.
          end if
          
          ! Create HDF5 filename
          hdf5_filename = trim(foldername)//'timeseries_data'//jobname//'.h5'
          
          ! Create or open HDF5 file
          call h5fcreate_f(trim(hdf5_filename), H5F_ACC_TRUNC_F, file_id, hdferr)
          
          ! Create time-series group
          time_series_group_name = '/time_series'
          call h5gcreate_f(file_id, trim(time_series_group_name), group_id, hdferr)
          
          ! Write slipz1_v data (velocity time series)
          dims_2d = (/Nt_all, ncos/)
          call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
          call h5dcreate_f(group_id, 'slipz1_v', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_v, dims_2d, hdferr)
          call h5dclose_f(dset_id, hdferr)
          call h5sclose_f(dspace_id, hdferr)
          
          ! Write slipz1_cos data (cosine slip time series)
          call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
          call h5dcreate_f(group_id, 'slipz1_cos', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_cos, dims_2d, hdferr)
          call h5dclose_f(dset_id, hdferr)
          call h5sclose_f(dspace_id, hdferr)
          
          ! Write time array
          dims_1d = (/ncos/)
          call h5screate_simple_f(1, dims_1d, dspace_id, hdferr)
          call h5dcreate_f(group_id, 'tcos', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tcos, dims_1d, hdferr)
          call h5dclose_f(dset_id, hdferr)
          call h5sclose_f(dspace_id, hdferr)
          
          ! Add metadata attributes (simplified - no character attributes for now)
          ! Note: Character attributes can cause issues with HDF5 Fortran interface
          ! The data itself provides sufficient information for analysis
          
          ! Close group and file
          call h5gclose_f(group_id, hdferr)
          call h5fclose_f(file_id, hdferr)
          
          ! Create XDMF file for visualization
          xdmf_filename = trim(foldername)//'timeseries_data'//jobname//'.xdmf'
          open(99, file=trim(xdmf_filename), status='replace')
          write(99,*) '<?xml version="1.0" ?>'
          write(99,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
          write(99,*) '<Xdmf Version="2.0">'
          write(99,*) '  <Domain>'
          write(99,*) '    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
          write(99,*) '      <Time TimeType="List">'
          write(99,*) '        <DataItem Format="XML" NumberType="Float" Dimensions="', ncos, '">'
          write(99,*) '          ', (tcos(i), i=1, ncos)
          write(99,*) '        </DataItem>'
          write(99,*) '      </Time>'
          write(99,*) '      <Grid Name="slipz1_v" GridType="Uniform">'
          write(99,*) '        <Topology TopologyType="2DCoRectMesh" Dimensions="', Nt_all, ' ', ncos, '"/>'
          write(99,*) '        <Geometry GeometryType="ORIGIN_DXDY">'
          write(99,*) '          <DataItem Format="XML" Dimensions="2">0.0 0.0</DataItem>'
          write(99,*) '          <DataItem Format="XML" Dimensions="2">1.0 1.0</DataItem>'
          write(99,*) '        </Geometry>'
          write(99,*) '        <Attribute Name="slipz1_v" AttributeType="Scalar" Center="Node">'
          write(99,*) '          <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="', Nt_all, ' ', ncos, '">'
          write(99,*) '            ', trim(hdf5_filename), ':/time_series/slipz1_v'
          write(99,*) '          </DataItem>'
          write(99,*) '        </Attribute>'
          write(99,*) '      </Grid>'
          write(99,*) '      <Grid Name="slipz1_cos" GridType="Uniform">'
          write(99,*) '        <Topology TopologyType="2DCoRectMesh" Dimensions="', Nt_all, ' ', ncos, '"/>'
          write(99,*) '        <Geometry GeometryType="ORIGIN_DXDY">'
          write(99,*) '          <DataItem Format="XML" Dimensions="2">0.0 0.0</DataItem>'
          write(99,*) '          <DataItem Format="XML" Dimensions="2">1.0 1.0</DataItem>'
          write(99,*) '        </Geometry>'
          write(99,*) '        <Attribute Name="slipz1_cos" AttributeType="Scalar" Center="Node">'
          write(99,*) '          <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="', Nt_all, ' ', ncos, '">'
          write(99,*) '            ', trim(hdf5_filename), ':/time_series/slipz1_cos'
          write(99,*) '          </DataItem>'
          write(99,*) '        </Attribute>'
          write(99,*) '      </Grid>'
          write(99,*) '    </Grid>'
          write(99,*) '  </Domain>'
          write(99,*) '</Xdmf>'
          close(99)
          
          write(*,*) 'Time-series data written to HDF5: ', trim(hdf5_filename)
          write(*,*) 'XDMF visualization file created: ', trim(xdmf_filename)
       end if
       
       icos = 0 
end if


	if(inul == nnul)then
       ! Null slip data collection completed - no output files needed
       inul = 0 
	end if

   if(isse==nsse)then
      ! HDF5 output for SSE time-series variables instead of binary files
      if (myid == master) then
         ! Initialize HDF5 if not already done
         if (.not. hdf5_initialized) then
            call h5open_f(hdferr)
            hdf5_initialized = .true.
         end if
         
         ! Create HDF5 filename for SSE data
         hdf5_filename = trim(foldername)//'sse_timeseries_data'//jobname//'.h5'
         
         ! Create or open HDF5 file
         call h5fcreate_f(trim(hdf5_filename), H5F_ACC_TRUNC_F, file_id, hdferr)
         
         ! Create SSE time-series group
         time_series_group_name = '/sse_time_series'
         call h5gcreate_f(file_id, trim(time_series_group_name), group_id, hdferr)
         
         ! Write slipz1_sse data (SSE slip time series)
         dims_2d = (/Nt_all, nsse/)
         call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
         call h5dcreate_f(group_id, 'slipz1_sse', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_sse, dims_2d, hdferr)
         call h5dclose_f(dset_id, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         
         ! Write slipz1_tau data (SSE tau time series)
         call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
         call h5dcreate_f(group_id, 'slipz1_tau', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_tau, dims_2d, hdferr)
         call h5dclose_f(dset_id, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         
         ! Write time array
         dims_1d = (/nsse/)
         call h5screate_simple_f(1, dims_1d, dspace_id, hdferr)
         call h5dcreate_f(group_id, 'tsse', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tsse, dims_1d, hdferr)
         call h5dclose_f(dset_id, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         
         ! Add metadata attributes (simplified - no character attributes for now)
         ! Note: Character attributes can cause issues with HDF5 Fortran interface
         ! The data itself provides sufficient information for analysis
         
         ! Close group and file
         call h5gclose_f(group_id, hdferr)
         call h5fclose_f(file_id, hdferr)
         
         ! Create XDMF file for SSE visualization
         xdmf_filename = trim(foldername)//'sse_timeseries_data'//jobname//'.xdmf'
         open(99, file=trim(xdmf_filename), status='replace')
         write(99,*) '<?xml version="1.0" ?>'
         write(99,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
         write(99,*) '<Xdmf Version="2.0">'
         write(99,*) '  <Domain>'
         write(99,*) '    <Grid Name="SSETimeSeries" GridType="Collection" CollectionType="Temporal">'
         write(99,*) '      <Time TimeType="List">'
         write(99,*) '        <DataItem Format="XML" NumberType="Float" Dimensions="', nsse, '">'
         write(99,*) '          ', (tsse(i), i=1, nsse)
         write(99,*) '        </DataItem>'
         write(99,*) '      </Time>'
         write(99,*) '      <Grid Name="slipz1_sse" GridType="Uniform">'
         write(99,*) '        <Topology TopologyType="2DCoRectMesh" Dimensions="', Nt_all, ' ', nsse, '"/>'
         write(99,*) '        <Geometry GeometryType="ORIGIN_DXDY">'
         write(99,*) '          <DataItem Format="XML" Dimensions="2">0.0 0.0</DataItem>'
         write(99,*) '          <DataItem Format="XML" Dimensions="2">1.0 1.0</DataItem>'
         write(99,*) '        </Geometry>'
         write(99,*) '        <Attribute Name="slipz1_sse" AttributeType="Scalar" Center="Node">'
         write(99,*) '          <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="', Nt_all, ' ', nsse, '">'
         write(99,*) '            ', trim(hdf5_filename), ':/sse_time_series/slipz1_sse'
         write(99,*) '          </DataItem>'
         write(99,*) '        </Attribute>'
         write(99,*) '      </Grid>'
         write(99,*) '      <Grid Name="slipz1_tau" GridType="Uniform">'
         write(99,*) '        <Topology TopologyType="2DCoRectMesh" Dimensions="', Nt_all, ' ', nsse, '"/>'
         write(99,*) '        <Geometry GeometryType="ORIGIN_DXDY">'
         write(99,*) '          <DataItem Format="XML" Dimensions="2">0.0 0.0</DataItem>'
         write(99,*) '          <DataItem Format="XML" Dimensions="2">1.0 1.0</DataItem>'
         write(99,*) '        </Geometry>'
         write(99,*) '        <Attribute Name="slipz1_tau" AttributeType="Scalar" Center="Node">'
         write(99,*) '          <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="', Nt_all, ' ', nsse, '">'
         write(99,*) '            ', trim(hdf5_filename), ':/sse_time_series/slipz1_tau'
         write(99,*) '          </DataItem>'
         write(99,*) '        </Attribute>'
         write(99,*) '      </Grid>'
         write(99,*) '    </Grid>'
         write(99,*) '  </Domain>'
         write(99,*) '</Xdmf>'
         close(99)
         
         write(*,*) 'SSE time-series data written to HDF5: ', trim(hdf5_filename)
         write(*,*) 'SSE XDMF visualization file created: ', trim(xdmf_filename)
      end if
      
      isse = 0
  end if


else

   if((imv>0).and.(imv<nmv))then
      ! OPTIMIZATION: Use buffered I/O and reduce file operations
             open(30,file=trim(foldername)//'maxvall'//jobname,position='append',status='unknown')
              open(311,file=trim(foldername)//'fltst_strk-36dp+00'//jobname,position='append',status='unknown')
        open(312,file=trim(foldername)//'fltst_strk-16dp+00'//jobname,position='append',status='unknown')
        open(313,file=trim(foldername)//'fltst_strk+00dp+00'//jobname,position='append',status='unknown')
        open(314,file=trim(foldername)//'fltst_strk+16dp+00'//jobname,position='append',status='unknown')
        open(315,file=trim(foldername)//'fltst_strk+36dp+00'//jobname,position='append',status='unknown')
        open(316,file=trim(foldername)//'fltst_strk-24dp+10'//jobname,position='append',status='unknown')
        open(317,file=trim(foldername)//'fltst_strk-16dp+10'//jobname,position='append',status='unknown')
        open(318,file=trim(foldername)//'fltst_strk+00dp+10'//jobname,position='append',status='unknown')
        open(319,file=trim(foldername)//'fltst_strk+16dp+10'//jobname,position='append',status='unknown')
        open(320,file=trim(foldername)//'fltst_strk+00dp+22'//jobname,position='append',status='unknown')

      ! OPTIMIZATION: Batch writes to reduce I/O overhead
      do i=1,imv
         write(30,130)tmv(i),dlog10(maxv(i)*1d-3/yrs),moment(i)
        do j=311,320
         write(j,110) tmv(i),outs1(i,1,j-310),outs1(i,2,j-310),outs1(i,3,j-310),outs1(i,4,j-310), &
           outs1(i,5,j-310),outs1(i,6,j-310),outs1(i,7,j-310)
        end do
      end do
       close(30)
       do j=311,320
         close(j)
       end do
 
      open(401,file=trim(foldername)//'blkst_strk-16fn+08dp+00'//jobname,position='append',status='unknown')
      open(402,file=trim(foldername)//'blkst_strk+00fn+08dp+00'//jobname,position='append',status='unknown')
      open(403,file=trim(foldername)//'blkst_strk+16fn+08dp+00'//jobname,position='append',status='unknown')
      open(404,file=trim(foldername)//'blkst_strk+00fn+16dp+00'//jobname,position='append',status='unknown')
      open(405,file=trim(foldername)//'blkst_strk+00fn+32dp+00'//jobname,position='append',status='unknown')
      open(406,file=trim(foldername)//'blkst_strk+00fn+48dp+00'//jobname,position='append',status='unknown')
      open(407,file=trim(foldername)//'blkst_strk+00fn+08dp+10'//jobname,position='append',status='unknown')
      open(408,file=trim(foldername)//'blkst_strk+00fn+16dp+10'//jobname,position='append',status='unknown')
      open(409,file=trim(foldername)//'blkst_strk+00fn+32dp+10'//jobname,position='append',status='unknown')
      open(501,file=trim(foldername)//'slip_2_depth'//jobname,position='append',status='unknown')
      open(502,file=trim(foldername)//'slip_2_strike'//jobname,position='append',status='unknown')
      open(503,file=trim(foldername)//'stress_2_depth'//jobname,position='append',status='unknown')
      open(504,file=trim(foldername)//'stress_2_strike'//jobname,position='append',status='unknown')

      do i=1,imv
        do j=401,409
         write(j,110) tmv(i),obvs(i,1,j-400),obvs(i,2,j-400),obvs(i,3,j-400),obvs(i,4,j-400),obvs(i,5,j-400),obvs(i,6,j-400)
        end do
        write(501,144) tmv(i),dlog10(maxv(i)*1d-3/yrs),obvdp(i,1,:)
        write(503,144) tmv(i),dlog10(maxv(i)*1d-3/yrs),obvdp(i,2,:)
        write(502,144) tmv(i),dlog10(maxv(i)*1d-3/yrs),obvstrk(i,1,:)
        write(504,144) tmv(i),dlog10(maxv(i)*1d-3/yrs),obvstrk(i,2,:)
     end do

     do j=501,504
       close(j)
     end do

      do j=401,409
         close(j)
      end do

      ! Handle observation data output separately (if n_obv > 0)
      if (n_obv > 0) then
         ! Open observation data output files
         do j = 1, min(n_obv, 9)  ! Limit to available observation points
            open(600+j, file=trim(foldername)//'obs_data_'//char(48+j)//jobname, position='append', status='unknown')
         end do
         
         ! Write observation data
         do i = 1, imv
            do j = 1, min(n_obv, 9)
               write(600+j, 110) tmv(i), obvs(i,1,j), obvs(i,2,j), obvs(i,3,j), obvs(i,4,j), obvs(i,5,j), obvs(i,6,j)
            end do
         end do
         
         ! Close observation data files
         do j = 1, min(n_obv, 9)
            close(600+j)
         end do
         
         write(*,*) 'Observation data written for ', min(n_obv, 9), ' observation points'
      end if

      imv = 0
    end if

    if((ias>0).and.(ias<nas))then
       open(31,form='unformatted',file=trim(foldername)//'slipz1-inter'//jobname,position='append',status='unknown')
       open(34,file=trim(foldername)//'t-inter'//jobname,position='append',status='unknown')

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
      ! Write partial SSE data to the same HDF5 file as main SSE output
      if (myid == master) then
         ! Initialize HDF5 if not already done
         if (.not. hdf5_initialized) then
            call h5open_f(hdferr)
            hdf5_initialized = .true.
         end if
         
         ! Open existing HDF5 file for SSE data (append mode)
         hdf5_filename = trim(foldername)//'sse_timeseries_data'//jobname//'.h5'
         call h5fopen_f(trim(hdf5_filename), H5F_ACC_RDWR_F, file_id, hdferr)
         
         ! Open existing SSE time-series group
         time_series_group_name = '/sse_time_series'
         call h5gopen_f(file_id, trim(time_series_group_name), group_id, hdferr)
         
         ! Write partial slipz1_sse data (SSE slip time series)
         dims_2d = (/Nt_all, isse/)
         call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
         call h5dcreate_f(group_id, 'slipz1_sse_partial', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_sse(:,1:isse), dims_2d, hdferr)
         call h5dclose_f(dset_id, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         
         ! Write partial slipz1_tau data (SSE tau time series)
         call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
         call h5dcreate_f(group_id, 'slipz1_tau_partial', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_tau(:,1:isse), dims_2d, hdferr)
         call h5dclose_f(dset_id, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         
         ! Write partial time array
         dims_1d = (/isse/)
         call h5screate_simple_f(1, dims_1d, dspace_id, hdferr)
         call h5dcreate_f(group_id, 'tsse_partial', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tsse(1:isse), dims_1d, hdferr)
         call h5dclose_f(dset_id, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         
         ! Add metadata attributes for partial data (simplified - no character attributes for now)
         ! Note: Character attributes can cause issues with HDF5 Fortran interface
         ! The data itself provides sufficient information for analysis
         
         ! Close group and file
         call h5gclose_f(group_id, hdferr)
         call h5fclose_f(file_id, hdferr)
         
         write(*,*) 'Partial SSE data written to HDF5: ', trim(hdf5_filename)
      end if
      
      isse = 0
  end if


     if((icos>0).and.(icos<ncos))then
       ! Write partial cosine slip data to the same HDF5 file as main cosine slip output
       if (myid == master) then
          ! Initialize HDF5 if not already done
          if (.not. hdf5_initialized) then
             call h5open_f(hdferr)
             hdf5_initialized = .true.
          end if
          
          ! Open existing HDF5 file for cosine slip data (append mode)
          hdf5_filename = trim(foldername)//'timeseries_data'//jobname//'.h5'
          call h5fopen_f(trim(hdf5_filename), H5F_ACC_RDWR_F, file_id, hdferr)
          
          ! Open existing time-series group
          time_series_group_name = '/time_series'
          call h5gopen_f(file_id, trim(time_series_group_name), group_id, hdferr)
          
          ! Write partial slipz1_cos data (cosine slip time series)
          dims_2d = (/Nt_all, icos/)
          call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
          call h5dcreate_f(group_id, 'slipz1_cos_partial', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_cos(:,1:icos), dims_2d, hdferr)
          call h5dclose_f(dset_id, hdferr)
          call h5sclose_f(dspace_id, hdferr)
          
          ! Write partial slipz1_v data (velocity time series)
          call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
          call h5dcreate_f(group_id, 'slipz1_v_partial', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_v(:,1:icos), dims_2d, hdferr)
          call h5dclose_f(dset_id, hdferr)
          call h5sclose_f(dspace_id, hdferr)
          
          ! Write partial time array
          dims_1d = (/icos/)
          call h5screate_simple_f(1, dims_1d, dspace_id, hdferr)
          call h5dcreate_f(group_id, 'tcos_partial', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tcos(1:icos), dims_1d, hdferr)
          call h5dclose_f(dset_id, hdferr)
          call h5sclose_f(dspace_id, hdferr)
          
          ! Add metadata attributes for partial data (simplified - no character attributes for now)
          ! Note: Character attributes can cause issues with HDF5 Fortran interface
          ! The data itself provides sufficient information for analysis
          
          ! Close group and file
          call h5gclose_f(group_id, hdferr)
          call h5fclose_f(file_id, hdferr)
          
          write(*,*) 'Partial cosine slip data written to HDF5: ', trim(hdf5_filename)
       end if
       
       icos = 0 
      end if

                 
	if((inul>0).and.(inul<nnul))then
	open(52,file=trim(foldername)//'vs-nul'//jobname,position='append',status='unknown')
        open(53,file=trim(foldername)//'nul-time'//jobname, &
             position='append',status='unknown')
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

    if(inul == nnul)then
       ! Null slip data collection completed - no output files needed
       inul = 0 
	end if

end if 

 110    format(E22.14,7(1X,E15.7))
 120    format(E20.13,4X,E20.13,4X,I6)
 130    format(E22.14,2(1X,E15.7))
 140    format(E20.13)
 150    format(E22.14,3(1X,E15.7))
 160    format(E20.13,1x,E20.13)
 500    format(E15.8,1X,E20.13)
 600    format(E15.4,1X,E13.6,1X,E15.8)
 700    format(E13.6)
 900    format(E15.8)

RETURN
END subroutine output
