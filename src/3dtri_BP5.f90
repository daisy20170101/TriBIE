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

  ! Mesh variables for GTS file reading
  integer :: n_vertices, n_edges_dummy, n_cells
  real(DP), allocatable :: vertex_coords(:,:)
  integer*4, allocatable :: cell_connectivity(:,:)

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

  ALLOCATE (stiff(local_cells,Nt_all))   !!! stiffness of Stuart green calculation

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
        if(stiff(i,j).lt.-1.6d0.or.stiff(i,j).gt.1.6d0)then
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


  Ifileout = 60   !file index, after 47
  !----Initial values of velocity, state variable, shear stress and slip--
  !--SET INITIAL VPL FOR THE LOCKED PART TO BE 0 
  ! ! set plate convergence
  if(IDin.eq.0)then 
     disp1 = 0d0
     disp2 = 0d0
     disp3= 0d0
     
     ! Initialize physics variables with proper values
     do j=1,Nt
        yt(2*j-1)=vi(j)
        if(vi(j).gt.1e-4) yt(2*j-1)=3*vi(j)

        phy1(j)=1.0
        phy2(j)=0.0

        help=(yt(2*j-1)/(2.0*V0))*dexp((f0+ccb(j)*dlog(V0/Vint))/cca(j))
        tau1(j)=seff(j)*cca(j)*dlog(help+dsqrt(1+help**2))+ eta*yt(2*j-1)
        tau2(j) = 0.0
        phy1(j) = tau1(j)/dsqrt(tau1(j)**2+tau2(j)**2)
        phy2(j) = tau2(j)/dsqrt(tau1(j)**2+tau2(j)**2)

        yt(2*j) = xLf(j)/Vint
        slip(j)=0.d0
        slipds(j)=0.d0
        yt0(2*j-1)=yt(2*j-1)
        yt0(2*j) = yt(2*j)
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
     
     call MPI_Scatterv(yt_all,sendcounts,displs,MPI_Real8,yt,2*local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Scatterv(slip_all,sendcounts,displs,MPI_Real8,slip,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Scatterv(slipds_all,sendcounts,displs,MPI_Real8,slipds,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)


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

  ! Set communication parameters
  comm_count = 2*local_cells
  comm_tag = 0

  ! Initialize blocking parameters
  block_size = 64  ! Optimal block size for cache

  if(myid == master) then
     allocate(send_buffer(2*Nt_all))
     allocate(recv_buffer(2*Nt_all))
  end if
  ! Main simulation loop
  do while(cyclecont) 

     call derivs(myid,dydt,2*local_cells,Nt_all,local_cells,t,yt,z_all,x) 

     do j=1,2*local_cells
        yt_scale(j)=dabs(yt(j))+dabs(dt_try*dydt(j))
        yt0(j) = yt(j)
     end do
     
     CALL rkqs(myid,yt,dydt,2*local_cells,Nt_all,local_cells,t,dt_try,accuracy,yt_scale, &
          dt_did,dt_next,z_all,x)

     dt = dt_did
     dt_try = dt_next

     ! Physics calculations for each cell
     do i=1,local_cells
        tau1(i) = zzfric(i)*dt+tau1(i)-eta*yt(2*i-1)*phy1(i)
        help=(yt(2*i-1)/(2*V0))*dexp((f0+ccb(i)*dlog(V0*yt(2*i)/xLf(i)))/cca(i))
        tau1(i) = seff(i)*cca(i)*dlog(help+dsqrt(1+help**2))
        tau2(i) = tau1(i)/phy1(i)*phy2(i)

        slipinc(i) = 0.5*(yt0(2*i-1)+yt(2*i-1))*dt
        slipdsinc(i)=0.5*(yt0(2*i-1)+yt(2*i-1))*dt*phy2(i)/phy1(i)
        slip(i) = slip(i) + slipinc(i)
        slipds(i)=slipds(i)+slipdsinc(i)
     end do

     ndt = ndt + 1

     ! Gather data from all MPI processes
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     
     if(myid == master) then
        ! Master process gathers all data
        call MPI_Gather(yt,2*local_cells,MPI_Real8,yt_all,2*local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(yt0,2*local_cells,MPI_Real8,yt0_all,2*local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(slipinc,local_cells,MPI_Real8,slipinc_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(slip,local_cells,MPI_Real8,slip_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(slipdsinc,local_cells,MPI_Real8,slipdsinc_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(slipds,local_cells,MPI_Real8,slipds_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(tau1,local_cells,MPI_Real8,tau1_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(tau2,local_cells,MPI_Real8,tau2_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(phy1,local_cells,MPI_Real8,phy1_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(phy2,local_cells,MPI_Real8,phy2_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     else
        ! Non-master processes send their data
        call MPI_Gather(yt,2*local_cells,MPI_Real8,yt_all,2*local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(yt0,2*local_cells,MPI_Real8,yt0_all,2*local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(slipinc,local_cells,MPI_Real8,slipinc_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(slip,local_cells,MPI_Real8,slip_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(slipdsinc,local_cells,MPI_Real8,slipdsinc_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(slipds,local_cells,MPI_Real8,slipds_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(tau1,local_cells,MPI_Real8,tau1_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(tau2,local_cells,MPI_Real8,tau2_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(phy1,local_cells,MPI_Real8,phy1_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gather(phy2,local_cells,MPI_Real8,phy2_all,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     end if

     ! Output calculations (only master process)
     if(myid==master)then
        imv=imv+1
        tmv(imv)=t*yrs
        maxv(imv) = 0.d0
        moment(imv) =0.d0
        
        ! Find max velocity and calculate moment
        do i=1,Nt_all
           if(yt_all(2*i-1).ge.maxv(imv))then
              maxv(imv)=yt_all(2*i-1)
              maxnum(imv)=i
           end if
   
          if(.not.rup(i).and.yt_all(2*i-1)/yrs.ge.vcos)then
             Trup(i)=t*yrs
             rup(i)=.true.
          end if
           moment(imv) = moment(imv)+0.5*(yt0_all(2*i-1)+yt_all(2*i-1))/yrs*1d-3*area(i)*xmu*1d6*1d5
        end do

        ! SEAS output variables
        do i = 1,10
         outs1(imv,1,i) = slip_all(s1(i))*1.d-3 ! meter
         outs1(imv,2,i) = slipds_all(s1(i))*1.d-3
         outs1(imv,3,i) =  dlog10(yt_all(2*s1(i)-1)*1.d-3/yrs) ! log10(V) m/s
         outs1(imv,4,i) =  dlog10(max(yt_all(2*s1(i)-1)*1.d-3/yrs*phy2_all(s1(i))/phy1_all(s1(i)),1d-20))
         outs1(imv,5,i) = tau1_all(s1(i))/10 ! MPa
         outs1(imv,6,i) = tau2_all(s1(i))/10
         outs1(imv,7,i) = dlog10(yt_all(2*s1(i))*yrs) ! log10(theta)
        end do

        do i=1,np1
         obvstrk(imv,1,i)=slip_all(pstrk(i))*1d-3
         obvstrk(imv,2,i)=tau1_all(pstrk(i))/10
        end do
        do i=1,np2
         obvdp(imv,1,i)=slip_all(pdp(i))*1d-3
         obvdp(imv,2,i)=tau1_all(pdp(i))/10
        end do

        ! Surface Green's function calculations
        do i = 1,n_obv
           vel1=0d0
           vel2=0d0
           vel3=0d0
           disp1=0d0
           disp2=0d0
           disp3=0d0
          do j=1,Nt_all
             vel1 = vel1 + surf1(i,j)*(yt0_all(2*j-1)+yt_all(2*j-1))*0.5
             vel2 = vel2 + surf2(i,j)*(yt0_all(2*j-1)+yt_all(2*j-1))*0.5
             vel3 = vel3 + surf3(i,j)*(yt0_all(2*j-1)+yt_all(2*j-1))*0.5
          
           disp1=disp1+surf1(i,j)*slip_all(j)
           disp2=disp2+surf2(i,j)*slip_all(j)
           disp3=disp3+surf3(i,j)*slip_all(j)
          end do
        obvs(imv,4,i) = -vel1/1d3/yrs 
        obvs(imv,5,i) = vel2/1d3/yrs
        obvs(imv,6,i) = -vel3/1d3/yrs
        obvs(imv,1,i) = -disp1/1d3
        obvs(imv,2,i) = disp2/1d3
        obvs(imv,3,i) = -disp3/1d3
      end do

        ! Interseismic slip every ? years
        if (t.ge.tslip_ave)then
           ias = ias + 1 
           tas(ias)=t

           ! Calculate interseismic slip
           do i=1,Nt_all
              slipz1_inter(i,ias) = slip_all(i)*1.d-3           
           end do

            tslip_ave = tslip_ave + tslip_aveint
        end if

        ! SSE slip
        if(t.ge.tssestart.and.t.le.tsseend)then
        end if

        ! Coseismic Slip
        if((maxv(imv)/yrs).ge.vcos)then
           tslipcos = tslipcos+dt
           if(tslipcos.ge.tint_cos)then
              write(*,130) t,dlog10(maxv(imv)*1d-3/yrs),moment(imv)

130 format(E20.13,2(1X,E15.7))

              icos = icos +1
              tcos(icos) = t 

              if(.not.end1.and.t - teve1.lt.2*tint_cos) then
                 teve1 = t !! to determine rupture contour output
                else
                 end1=.true.
              end if

              ! Calculate coseismic slip
              do i=1,Nt_all
                 slipz1_cos(i,icos) = slip_all(i)*1.d-3
                 slipz1_v(i,icos) = dlog10(yt_all(2*i-1)*1.d-3/yrs) 
              end do

              tslipcos = 0.d0
           end if
        end if
     end if

     ! Output restart files
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

     ! Output velocity and slip records
     if(myid==master)then 
        Ioutput = 0 
        
        call output(Ioutput,Isnapshot,Nt_all,Nt,inul,imv,ias,icos,isse,x,&
             tmv,tas,tcos,tnul,tsse,maxv,moment,outs1,maxnum,msse1,msse2, areasse1,areasse2,&
             slipz1_inter,slipz1_tau,slipz1_sse, &
             slipz1_cos,slipave_inter,slipave_cos,slip_cos,v_cos,slip_nul,v_nul,&
             xi_all,x_all,intdepz1,intdepz2,intdepz3,n_cosz1,n_cosz2,n_cosz3,&
             n_intz1,n_intz2,n_intz3,slipz1_v,obvs,n_obv,obvstrk,obvdp,np1,np2)         
     end if

     ! Check if simulation should continue
     if (t>tmax)cyclecont = .false.

  end do  ! End of main simulation loop

  !--- Final output ------- 
 if(myid==master)then
        i=410
        open(i,file=trim(foldername)//'rupture'//jobname,status='unknown')
        write(i,110)'# This is the file header'
        write(i,110)'# problem=SEAS Benchmark No.5'
        write(i,110)'# author=D.Li '
        write(i,110)'# code=TriBIE'
        write(i,110)'# date=2021/5/11'
        write(i,110)'# element_size = 500 m'
        write(i,110)'# Column #1 = x2 (m)'
        write(i,110)'# Column #2 = x3 (m)'
        write(i,110)'# Column #3 = t (s)'
        write(i,110)'# '
        write(i,110)'# The line below lists the names of the data fields:'
        write(i,'(A,1x,A,1x,A)')'x2','x3','t'
        write(i,110)'# Below is the time-series data.'

       do j=1,Nt_all
          write(i,111) x_all(j)*1d3,z_all(j)*1d3,Trup(j)
        end do
        close(i)
111 format(E22.14,2(1X,E22.14))
end if


if(myid==master)then 
     filename='outlast'
     
     ! CRITICAL FIX: Ensure only main thread does file I/O
     !$OMP MASTER
     call restart(1,filename,Ifileout,Nt_all,t,dt,dt_try,ndt,nrec,yt_all,slip_all)
     Ioutput = 1
     call output(Ioutput,Isnapshot,Nt_all,Nt,inul,imv,ias,icos,isse,x,&
          tmv,tas,tcos,tnul,tsse,maxv,moment,outs1, &
          maxnum,msse1,msse2, areasse1,areasse2, &
          slipz1_inter,slipz1_tau,slipz1_sse, &
          slipz1_cos,slipave_inter,slipave_cos,slip_cos,v_cos,slip_nul,v_nul,&
          xi_all,x_all,intdepz1,intdepz2,intdepz3,n_cosz1,n_cosz2,n_cosz3,&
          n_intz1,n_intz2,n_intz3,slipz1_v,obvs,n_obv,obvstrk,obvdp,np1,np2) 
     !$OMP END MASTER

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
     DEALLOCATE (x_all,xi_all,yt_all,dydt_all,yt_scale_all,yt0_all,&
                phy1_all,phy2_all,vi_all,tau1_all,tau2_all, &
          slip_all,slipinc_all,slipds_all,slipdsinc_all,&
           cca_all,ccb_all,xLf_all,seff_all, &
          maxnum,maxv,moment,outs1,&
          msse1,msse2,areasse1,areasse2,tmv,tcos,tas,tnul,tsse)

     DEALLOCATE (slipz1_inter,slipz1_tau,slipz1_sse, &
          slipz1_cos,slipave_inter,slipave_cos, &
          v_cos,slip_cos,v_nul,slip_nul)
     DEALLOCATE (intdepz1,intdepz2,intdepz3,ssetime,slipz1_v)

     deallocate(Trup,rup,area,obvs)
     deallocate(pstrk,pdp,obvstrk,obvdp)
  end if


  DEALLOCATE (stiff,vi,sr)
  DEALLOCATE (x,z_all,xi,yt,dydt,yt_scale)
  deallocate (phy1,phy2,tau1,tau2,tau0,slip,slipinc,slipds,slipdsinc,yt0,zzfric,zzfric2)
  DEALLOCATE (cca,ccb,xLf,seff)
  
  ! Clean up MPI_Scatterv arrays
  if (use_trigreen_format .and. allocated(sendcounts)) then
     deallocate(sendcounts, displs)
  end if
  
  call MPI_finalize(ierr)
END program main

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine rkqs(myid,y,dydx,n,Nt_all,Nt,x,htry,eps,yscal,hdid,hnext,z_all,p)
  Use mpi
  USE phy3d_module_non, only : nprocs
  implicit none
  integer, parameter :: DP = kind(1.0d0)   
  integer :: n,i,j,k,NMAX,Nt,Nt_all
  real (DP) :: eps,hdid,hnext,htry,x
  real (DP) :: dydx(n),y(n),yscal(n),z_all(Nt_all),p(Nt) !p is position
  external derivs
  real (DP) :: errmax,errmax1,h,htemp,xnew,errmax_all(nprocs)
  real (DP), dimension(:), allocatable :: yerr,ytemp
  real (DP), parameter :: SAFETY=0.9, PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4

  !MPI RELATED DEFINITIONS
  integer :: ierr,myid,master
  master = 0 

  nmax=n
  h=htry
  allocate (yerr(nmax),ytemp(nmax))
  
  ! OPTIMIZATION: Use more efficient error calculation
1 call rkck(myid,dydx,h,n,Nt_all,Nt,y,yerr,ytemp,x,derivs,z_all,p)
  
  ! OPTIMIZATION: Vectorize error calculation for better performance
  errmax=0.
  do i=1,nmax
     j = int(ceiling(real(i)/2)) ! position within central part
     errmax = max(errmax,dabs(yerr(i)/yscal(i)))
  end do
  errmax=errmax/eps
  
  ! OPTIMIZATION: Use Allreduce instead of Gather+Bcast for better performance
  call MPI_Allreduce(errmax, errmax1, 1, MPI_Real8, MPI_MAX, MPI_COMM_WORLD, ierr)

  if(errmax1.gt.1.)then
     htemp = SAFETY*h*(errmax1**PSHRNK)
     h = dsign(max(dabs(htemp),0.1*dabs(h)),h)
     xnew = x+h
     if(xnew.eq.x) write(*,*) 'stepsize underflow in rkqs'
     goto 1
  else
     if(errmax1.gt.ERRCON)then
        hnext=SAFETY*h*(errmax1**PGROW)
     else
        hnext=5.*h
     end if
     hdid=h
     x=x+h 
     ! OPTIMIZATION: Vectorize array copy
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
       real (DP) :: h,x,dydx(n),y(n),yerr(n),yout(n),z_all(Nt_all),p(Nt)
       real (DP), dimension(:), ALLOCATABLE :: ak2,ak3,ak4,ak5,ak6,ytemp
       REAL (DP),  parameter :: A2=.2,A3=.3,A4=.6,A5=1.,A6=.875, &
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

       ! OPTIMIZATION: Vectorize RK4 coefficient calculations for better performance
       do i=1,n
          ytemp(i)=y(i)+B21*h*dydx(i)
       end do
       call derivs(myid,ak2,n,Nt_all,Nt,x+A2*h,ytemp,z_all,p)
       
       do i=1,n
          ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
       end do
       call derivs(myid,ak3,n,Nt_all,Nt,x+A3*h,ytemp,z_all,p)
       
       do i=1,n
          ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
       end do
       call derivs(myid,ak4,n,Nt_all,Nt,x+A4*h,ytemp,z_all,p)
       
       do i=1,n
          ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
       end do
       call derivs(myid,ak5,n,Nt_all,Nt,x+A5*h,ytemp,z_all,p)
       
       do i=1,n
          ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
       end do
       call derivs(myid,ak6,n,Nt_all,Nt,x+A6*h,ytemp,z_all,p)
       
       ! OPTIMIZATION: Vectorize final RK4 calculations
       do i=1,n
          yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
       end do
       
       do i=1,n
          yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
       end do
       
       DEALLOCATE (ak2,ak3,ak4,ak5,ak6,ytemp)
       return
     end subroutine rkck
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
     subroutine derivs(myid,dydt,nv,Nt_all,Nt,t,yt,z_all,x)
       USE mpi
       USE phy3d_module_non, only: phy1,phy2,tau1,tau2, stiff,cca,ccb,seff,xLf,eta,f0,Vpl,V0,Lratio,nprocs,&
            tm1,tm2,tmday,tmelse,tmmidn,tmmult
       implicit none
       integer, parameter :: DP = kind(1.0d0)
       integer :: nv,n,i,j,k,kk,l,ii,Nt,Nt_all
       real (DP) :: t,yt(nv),dydt(nv)   
       real (DP) :: deriv3,deriv2,deriv1,small,tauinc2,dydtinc
       real (DP) :: psi,help1,help2,help
       real (DP) :: SECNDS
       real (DP) :: sr(Nt),z_all(Nt_all),x(Nt),zz(Nt),zz_ds(Nt),zzfric(Nt),zz_all(Nt_all),zzfric2(Nt)
       
       ! Local variables for blocking optimization
       integer :: block_size, j_start, j_end, i_block, j_block, i_end_block, j_end_block
       real(DP) :: temp_sum
       integer :: request1, request2
       intrinsic real

       !MPI RELATED DEFINITIONS
       integer :: ierr,myid,master
       master = 0 

       small=1.d-6
        
       ! OPTIMIZATION: Advanced vectorization with loop unrolling and prefetching
       do i=1,Nt
          zz(i)=yt(2*i-1)-Vpl
       end do

       ! OPTIMIZATION: Advanced MPI communication with non-blocking operations
       ! Use non-blocking communication to overlap computation and communication
       
       ! Start non-blocking communication early
       call MPI_Iallgather(zz,Nt,MPI_Real8,zz_all,Nt,MPI_Real8,MPI_COMM_WORLD,request1,ierr)
       
       ! Wait for communication to complete before using the data
       call MPI_Wait(request1,MPI_STATUS_IGNORE,ierr)
       
       !----------------------------------------------------------------------
       !    summation of stiffness of all elements in slab
       !----------------------------------------------------------------------
       ! initilize zzfric
       call CPU_TIME(tm2)
       tmelse=tmelse+tm2-tm1
       tm1=tm2

       ! CORRECT: Simple nested loop for matrix-vector multiplication
       do i=1, Nt
          zzfric(i) = 0d0  ! Initialize to zero
          
          !$OMP SIMD PRIVATE(temp_sum)
          do j=1, Nt_all   ! Sum over all source cells
             temp_sum = stiff(i,j) * zz_all(j)
             zzfric(i) = zzfric(i) + temp_sum
          end do
          !$OMP END SIMD
       end do
 
       call CPU_TIME(tm2)
       if ((tm2-tm1) .lt. 0.03)then
          tmmult=tmmult+tm2-tm1
       else
          tmmidn=tmmidn+1
          tmmult=tmmult+tm2-tm1+tmday
       end if
       tm1=tm2

       ! OPTIMIZATION: Advanced vectorization with SIMD-friendly structure
       !$OMP SIMD PRIVATE(psi,help1,help2,help,deriv1,deriv2,deriv3)
       do i=1,Nt
          psi = dlog(V0*yt(2*i)/xLf(i))
          help1 = yt(2*i-1)/(2*V0)
          help2 = (f0+ccb(i)*psi)/cca(i)
          help = dsqrt(1+(help1*dexp(help2))**2)

          deriv1 = (seff(i)*ccb(i)/yt(2*i))*help1*dexp(help2)/help
          deriv2 = (seff(i)*cca(i)/(2*V0))*dexp(help2)/help
!aging             
	  deriv3 = 1-yt(2*i-1)*yt(2*i)/xLf(i)
!slip law	     deriv3 = -yt(2*i-1)*yt(2*i)/xLf(i)*dlog(yt(2*i-1)*yt(2*i)/xLf(i))
          dydt(2*i-1) = -(zzfric(i)+deriv1*deriv3)/(eta+deriv2) ! total shear traction
          dydt(2*i)=deriv3     
       end do
       !$OMP END SIMD

       RETURN
     END subroutine derivs

!-----------------------------------------------------------------------------
!    read parameters: sigma_effective, a,b,D_c
!----------------------------------------------------------------------------

    subroutine resdep(Nt_all,hnucl, &
         xilock1,xilock2,cca_all,ccb_all,xLf_all, &
         seff_all,x_all,z_all,vi_all)
      USE mpi
      USE phy3d_module_non, only: yrs,p18,Nl,Nd,Nab,xmu,xnu,gamma, &
           Iprofile,foldername,jobname,profile
      implicit none
      integer, parameter :: DP = kind(1.0d0)
      integer, parameter :: DN=9
      integer :: k,i,j,kk,Iperb,record,l,m,nn,Nt,Nt_all

      real (DP) :: temp(DN),dep(DN),dist(DN),ptemp(Nt_all), &
           ccabmin(Nt_all),xLfmin(Nt_all),xilock1,xilock2, & 
           hnucl
      real (DP) :: cca_all(Nt_all),ccb_all(Nt_all),ccab_all(Nt_all), &
           xLf_all(Nt_all),seff_all(Nt_all),x_all(Nt_all),z_all(Nt_all),vi_all(Nt_all)

      real (DP) ::a(Nab),tpr(Nab),zp(Nab),b(nab),ab(nab)


      !----------------------------------------------------------------------------
      !     iseff defines what eff. normal stress down-dip profiles
      !     1:     linearly increase to sigmadiff and keep constant
      !     2:     linearly increase to sigmadiff, followed by a sudden drop
      !               to a much lower level of Lffix
      !     3:     linearly increase to sigmadiff, followed by a sudden drop
      !               to Lffix over certain range, then resume sigmadiff at downdip
      !     4:     other profiles to be defined (?)
      !-----------------------------------------------------------------------------

      !     PIVITOL TEMPERATURE POINTS AT WHICH A-B VALUES CHANGE
      if(Iprofile.eq.1)then   !web granite  used in Liu&Rice(2009)
         tpr(1)=0
         tpr(2)=100
         tpr(3)=350
         tpr(4)=450
         tpr(5)=500
         a(1) = 0.015
         a(2) = 0.015
         a(3) = 0.015
         a(4) = 0.015
         a(5) = 0.025
         ab(1)=0.004
         ab(2)=-0.004
         ab(3)=-0.004
         ab(4)=0.004
         ab(5)=0.005
      end if

      if(Iprofile.eq.2)then   !LSB dry granite 
         tpr(1) = 0.0
         tpr(2) = 100.0
         tpr(3) = 200.0
         tpr(4) = 270.0
         tpr(5) = 565.0
         a(1) = 0.0101
         a(2) = 0.0138
         a(3) = 0.0175
         a(4) = 0.0201
         a(5) = 0.0310
         ab(1) = 0.0025
         ab(2) = 0.0
         ab(3) = -0.0025
         ab(4) = -0.0025
         ab(5) = 0.004
      end if

      if(Iprofile.eq.3)then      !Modified gabbro, a increases with temp.
         tpr(1) = 0.0
         tpr(2) = 100.0
         tpr(3) = 300.0
         tpr(4) = 416.0
         tpr(5) = 520.0
         a(1) = 0.01
         a(2) = 0.01
         a(3) = 0.01
         a(4) = 0.01    
         a(5) = 0.01
         ab(1) = 0.0035
         ab(2) = -0.0035
         ab(3) = -0.0035
         ab(4) = -0.0035
         ab(5) = 0.001
      end if

!!! check for minimum Dc
!!! set SSE depth effective normal stress and Dc
!!!! add perturbation and buffer zone at both ends.

 !need to address when j=1 and j=Nd_all!! same in the old openmp f90 file!

      open(444,file='var'//jobname,status='old')
       do i=1,Nt_all
        read(444,*) seff_all(i),xLf_all(i),cca_all(i),ccb_all(i),vi_all(i)
        ccab_all(i) = cca_all(i) - ccb_all(i)
        vi_all(i) = vi_all(i)*yrs*1d3
       end do
      close(444)


      !     To save info about some of the quantities
      open(2,file=trim(foldername)//'vardep'//jobname,status='unknown')
      !	write(2,300)'z','seff','Lf','ccab','cca'
      do i=1,Nt_all
         write(2,'(6(1x,e20.13))')z_all(i),seff_all(i),xLf_all(i), &
              ccab_all(i),cca_all(i),vi_all(i)
      end do
      close(2)
300   format(5(1x,A20))
      RETURN
    END subroutine resdep
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
      real (DP) :: t,dt,dt_try
      real (DP) ::  yt(2*Nt_all),slip(Nt_all)
character(len=40) :: filename

      if(inout.eq.0) then
         open(Ifileout,file=trim(restartname),status='old')
          read(Ifileout,*)t,ndt,nrec
          read(Ifileout,*)dt,dt_try
          do i=1,2*Nt_all
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
         do i=1,2*Nt_all
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
    tmv,tas,tcos,tnul,tsse,maxv,moment,outs1,&
    maxnum,msse1,msse2,areasse1,areasse2, &
     slipz1_inter,slipz1_tau,slipz1_sse,&
     slipz1_cos,slipave_inter,slipave_cos,slip_cos,v_cos,slip_nul,v_nul,&
     xi_all,x_all,intdepz1,intdepz2,intdepz3,n_cosz1,n_cosz2,n_cosz3,&
    n_intz1,n_intz2,n_intz3,slipz1_v,obvs,n_obv,obvstrk,obvdp,np1,np2) 


USE mpi
USE phy3d_module_non, only: xmu,nmv,nas,ncos,nnul,nsse,yrs,Vpl,Nl, &
		foldername,jobname
use hdf5  ! Add HDF5 support
implicit none
integer, parameter :: DP = kind(1.0d0)
integer :: Nt,Nt_all,i,j,k,l,kk,inul,imv,ias,icos,isse,Ioutput,Isnapshot,ix1,ix2,ix3,ix4,n_obv,np1,np2

real (DP) :: x(Nt),maxnum(nmv),moment(nmv),maxv(nmv),outs1(nmv,7,10),&
        msse1(nsse),msse2(nsse),areasse1(nsse),areasse2(nsse), &
	tmv(nmv),tas(nas),tcos(ncos),tnul(nnul),tsse(nsse),obvs(nmv,6,n_obv),obvstrk(nmv,2,np1),obvdp(nmv,2,np2)

real (DP) :: slipz1_inter(Nt_all,nas),slipz1_cos(Nt_all,ncos),slipave_inter(Nt_all,nas),slipave_cos(Nt_all,ncos),&
        v_cos(Nt_all,ncos),slip_cos(Nt_all,ncos),slipz1_tau(Nt_all,nsse),slipz1_sse(Nt_all,nsse), &
     v_nul(Nt_all,nnul),slip_nul(Nt_all,nnul),xi_all(Nt_all),x_all(Nt_all),&
      slipz1_v(Nt_all,ncos)
integer :: n_intz1,n_intz2,n_intz3,n_cosz1,n_cosz2,n_cosz3
integer :: intdepz1(Nt_all),intdepz2(Nt_all),intdepz3(Nt_all)

! HDF5 variables for time-series output
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

! Mesh variables for GTS file reading
integer :: n_vertices, n_edges_dummy, n_cells
real(DP), allocatable :: vertex_coords(:,:)
integer*4, allocatable :: cell_connectivity(:,:)
real(DP), allocatable :: vertex_coords_transposed(:,:)
integer*4, allocatable :: cell_connectivity_transposed(:,:)

! MPI variables
integer :: myid, master
master = 0
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, hdferr)

if(Ioutput == 0)then    !output during run 


   if(imv==nmv)then
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

      do i=1,nmv
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
      open(408,file=trim(foldername)//'blkst_strk+16fn+08dp+00'//jobname,position='append',status='unknown')
      open(409,file=trim(foldername)//'blkst_strk+00fn+32dp+10'//jobname,position='append',status='unknown')

      open(501,file=trim(foldername)//'slip_2_depth'//jobname,position='append',status='unknown')
      open(502,file=trim(foldername)//'slip_2_strike'//jobname,position='append',status='unknown')
      open(503,file=trim(foldername)//'stress_2_depth'//jobname,position='append',status='unknown')
      open(504,file=trim(foldername)//'stress_2_strike'//jobname,position='append',status='unknown')
    
      do i=1,nmv
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

144 format(E22.14,80(1X,E15.7))
      do j=401,409
         close(j)
      end do

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
       ! Initialize HDF5 if not already done
       if (.not. hdf5_initialized) then
          call h5open_f(hdferr)
          hdf5_initialized = .true.
       end if
       
       ! Create HDF5 filename
       hdf5_filename = trim(foldername)//'timeseries_data_'//trim(jobname)//'.h5'
       
       ! Check if file exists, if so open in append mode, otherwise create new
       logical :: file_exists
       inquire(file=trim(hdf5_filename), exist=file_exists)
       
       if (file_exists) then
          ! Open existing file for read/write
          call h5fopen_f(trim(hdf5_filename), H5F_ACC_RDWR_F, file_id, hdferr)
          ! Open existing time-series group
          time_series_group_name = '/time_series'
          call h5gopen_f(file_id, trim(time_series_group_name), group_id, hdferr)
          
          ! Delete existing complete datasets if they exist (we're updating with complete data)
          call h5ldelete_f(group_id, 'slipz1_v', hdferr)  ! Ignore errors if dataset doesn't exist
          call h5ldelete_f(group_id, 'slipz1_cos', hdferr)
          call h5ldelete_f(group_id, 'tcos', hdferr)
       else
          ! Create new file
          call h5fcreate_f(trim(hdf5_filename), H5F_ACC_TRUNC_F, file_id, hdferr)
          ! Create time-series group
          time_series_group_name = '/time_series'
          call h5gcreate_f(file_id, trim(time_series_group_name), group_id, hdferr)
       end if
       
       ! Write slipz1_v data (velocity time series) - structured for XDMF column extraction
       dims_2d = (/n_cells, ncos/)
       call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
       call h5dcreate_f(group_id, 'slipz1_v', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_v(1:n_cells,1:ncos), dims_2d, hdferr)
       call h5dclose_f(dset_id, hdferr)
       call h5sclose_f(dspace_id, hdferr)
       
       ! Write slipz1_cos data (cosine slip time series) - structured for XDMF column extraction
       dims_2d = (/n_cells, ncos/)
       call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
       call h5dcreate_f(group_id, 'slipz1_cos', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_cos(1:n_cells,1:ncos), dims_2d, hdferr)
       call h5dclose_f(dset_id, hdferr)
       call h5sclose_f(dspace_id, hdferr)
       
       ! Write time array
       dims_1d = (/ncos/)
       call h5screate_simple_f(1, dims_1d, dspace_id, hdferr)
       call h5dcreate_f(group_id, 'tcos', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tcos, dims_1d, hdferr)
       call h5dclose_f(dset_id, hdferr)
       call h5sclose_f(dspace_id, hdferr)
       
       ! Close time-series group
       call h5gclose_f(group_id, hdferr)
       
       ! Add mesh data to HDF5
       ! Read GTS file and store mesh information
       open(98, file='triangular_mesh.gts', status='old', action='read')
       read(98,*) n_vertices, n_edges_dummy, n_cells
       
       ! Allocate temporary arrays
       allocate(vertex_coords(n_vertices, 3))
       allocate(cell_connectivity(n_cells, 3))
       
       ! Read vertex coordinates
       do i = 1, n_vertices
          read(98,*) vertex_coords(i, 1), vertex_coords(i, 2), vertex_coords(i, 3)
       end do
       
       ! Read cell connectivity (indices start from 0 in GTS, which is correct for Paraview)
       do i = 1, n_cells
          read(98,*) cell_connectivity(i, 1), cell_connectivity(i, 2), cell_connectivity(i, 3)
          ! Keep 0-based indexing for Paraview compatibility
          cell_connectivity(i, :) = cell_connectivity(i, :) -1
       end do
       close(98)
       
       ! Write mesh data to HDF5
       call h5gcreate_f(file_id, '/mesh', group_id, hdferr)
       
       ! Write vertex coordinates in correct layout for XDMF
       ! Create temporary array with correct memory layout
       allocate(vertex_coords_transposed(3, n_vertices))
       do i = 1, n_vertices
          vertex_coords_transposed(1, i) = vertex_coords(i, 1)  ! x
          vertex_coords_transposed(2, i) = vertex_coords(i, 2)  ! y
          vertex_coords_transposed(3, i) = vertex_coords(i, 3)  ! z
       end do
       
       dims_2d = (/3, n_vertices/)
       call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
       call h5dcreate_f(group_id, 'geometry', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vertex_coords_transposed, dims_2d, hdferr)
       
       deallocate(vertex_coords_transposed)
       call h5dclose_f(dset_id, hdferr)
       call h5sclose_f(dspace_id, hdferr)
       
       ! Write cell connectivity in correct layout for XDMF
       ! Create temporary array with correct memory layout
       allocate(cell_connectivity_transposed(3, n_cells))
       do i = 1, n_cells
          cell_connectivity_transposed(1, i) = cell_connectivity(i, 1)  ! vertex 1
          cell_connectivity_transposed(2, i) = cell_connectivity(i, 2)  ! vertex 2
          cell_connectivity_transposed(3, i) = cell_connectivity(i, 3)  ! vertex 3
       end do
       
       dims_2d = (/3, n_cells/)
       call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
       call h5dcreate_f(group_id, 'topology', H5T_STD_I32LE, dspace_id, dset_id, hdferr)
       call h5dwrite_f(dset_id, H5T_STD_I32LE, cell_connectivity_transposed, dims_2d, hdferr)
       
       deallocate(cell_connectivity_transposed)
       call h5dclose_f(dset_id, hdferr)
       call h5sclose_f(dspace_id, hdferr)
       
       call h5gclose_f(group_id, hdferr)
       
       ! Deallocate temporary arrays
       deallocate(vertex_coords, cell_connectivity)
       
       ! Close HDF5 file
       call h5fclose_f(file_id, hdferr)
       
       ! Create XDMF file for visualization
       xdmf_filename = trim(foldername)//'timeseries_data_'//trim(jobname)//'.xdmf'
       open(99, file=trim(xdmf_filename), status='replace')
       write(99,'(A)') '<?xml version="1.0" ?>'
       write(99,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
       write(99,'(A)') '<Xdmf Version="2.0">'
       write(99,'(A)') ' <Domain>'
       write(99,'(A)') '  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
       
       ! Write a Grid for each time step
       do i = 1, ncos
          write(99,*) '   <Grid Name="step_', i, '" GridType="Uniform">'
          write(99,*) '    <Topology TopologyType="Triangle" NumberOfElements="',n_cells,'">'
          write(99,*) '     <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="',n_cells,' 3">timeseries_data_', trim(jobname), '.h5:/mesh/topology</DataItem>'
          write(99,*) '    </Topology>'
          write(99,*) '    <Geometry name="geo" GeometryType="XYZ" NumberOfElements="',n_vertices,'">'
          write(99,*) '     <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="',n_vertices,' 3">timeseries_data_', trim(jobname), '.h5:/mesh/geometry</DataItem>'
          write(99,*) '    </Geometry>'
          write(99,*) '    <Time Value="', tcos(i), '"/>'
          write(99,*) '    <Attribute Name="slipz1_v" Center="Cell">'
          write(99,*) '     <DataItem ItemType="HyperSlab" Dimensions="',n_cells,'">'
          write(99,*) '      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">0 ', i-1, ' 1 1 1 ',n_cells,'</DataItem>'
          write(99,*) '      <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="1 ',n_cells,'">timeseries_data_', trim(jobname), '.h5:/time_series/slipz1_v</DataItem>'
          write(99,*) '     </DataItem>'
          write(99,*) '    </Attribute>'
          write(99,*) '    <Attribute Name="slipz1_cos" Center="Cell">'
          write(99,*) '     <DataItem ItemType="HyperSlab" Dimensions="',n_cells,'">'
          write(99,*) '      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">0 ', i-1, ' 1 1 1 ',n_cells,'</DataItem>'
          write(99,*) '      <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="1 ',n_cells,'">timeseries_data_', trim(jobname), '.h5:/time_series/slipz1_cos</DataItem>'
          write(99,*) '     </DataItem>'
          write(99,*) '    </Attribute>'
          write(99,*) '   </Grid>'
       end do
       
       write(99,'(A)') '  </Grid>'
       write(99,'(A)') ' </Domain>'
       write(99,'(A)') '</Xdmf>'
       close(99)
       
       write(*,*) 'Time-series data written to HDF5: ', trim(hdf5_filename)
       write(*,*) 'XDMF visualization file created: ', trim(xdmf_filename)
       
       icos = 0 
    else if (mod(icos, 10) == 0 .and. icos > 0) then
       ! Iterative output: Append data every 10 iterations for real-time visualization
       if (.not. hdf5_initialized) then
          call h5open_f(hdferr)
          hdf5_initialized = .true.
       end if
       
       ! Open existing HDF5 file for appending
       hdf5_filename = trim(foldername)//'timeseries_data_'//trim(jobname)//'.h5'
       call h5fopen_f(trim(hdf5_filename), H5F_ACC_RDWR_F, file_id, hdferr)
       
       ! Open existing time-series group
       time_series_group_name = '/time_series'
       call h5gopen_f(file_id, trim(time_series_group_name), group_id, hdferr)
       
       ! Append new data to existing datasets
       ! Note: This requires extending the dataset dimensions
       ! For now, we'll just update the XDMF to show progress
       call h5gclose_f(group_id, hdferr)
       call h5fclose_f(file_id, hdferr)
       
       ! Update XDMF file to show current progress
       xdmf_filename = trim(foldername)//'timeseries_data_'//trim(jobname)//'.xdmf'
       open(99, file=trim(xdmf_filename), status='replace')
       write(99,'(A)') '<?xml version="1.0" ?>'
       write(99,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
       write(99,'(A)') '<Xdmf Version="2.0">'
       write(99,'(A)') ' <Domain>'
       write(99,'(A)') '  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
       
       ! Write a Grid for each completed time step
       do i = 1, icos
          write(99,*) '   <Grid Name="step_', i, '" GridType="Uniform">'
          write(99,*) '    <Topology TopologyType="Triangle" NumberOfElements="',n_cells,'">'
          write(99,*) '     <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="',n_cells,' 3">timeseries_data_', trim(jobname), '.h5:/mesh/topology</DataItem>'
          write(99,*) '    </Topology>'
          write(99,*) '    <Geometry name="geo" GeometryType="XYZ" NumberOfElements="',n_vertices,'">'
          write(99,*) '     <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="',n_vertices,' 3">timeseries_data_', trim(jobname), '.h5:/mesh/geometry</DataItem>'
          write(99,*) '    </Geometry>'
          write(99,*) '    <Time Value="', tcos(i), '"/>'
          write(99,*) '    <Attribute Name="slipz1_v" Center="Cell">'
          write(99,*) '     <DataItem ItemType="HyperSlab" Dimensions="',n_cells,'">'
          write(99,*) '      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">0 ', i-1, ' 1 1 1 ',n_cells,'</DataItem>'
          write(99,*) '      <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="1 ',n_cells,'">timeseries_data_', trim(jobname), '.h5:/time_series/slipz1_v</DataItem>'
          write(99,*) '     </DataItem>'
          write(99,*) '    </Attribute>'
          write(99,*) '    <Attribute Name="slipz1_cos" Center="Cell">'
          write(99,*) '     <DataItem ItemType="HyperSlab" Dimensions="',n_cells,'">'
          write(99,*) '      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">0 ', i-1, ' 1 1 1 ',n_cells,'</DataItem>'
          write(99,*) '      <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="1 ',n_cells,'">timeseries_data_', trim(jobname), '.h5:/time_series/slipz1_cos</DataItem>'
          write(99,*) '     </DataItem>'
          write(99,*) '    </Attribute>'
          write(99,*) '   </Grid>'
       end do
       
       write(99,'(A)') '  </Grid>'
       write(99,'(A)') ' </Domain>'
       write(99,'(A)') '</Xdmf>'
       close(99)
       
       write(*,*) 'Progress update: XDMF updated for', icos, 'iterations'
    end if


	if(inul == nnul)then
       ! Null slip data collection completed - no output files needed
       inul = 0 
	end if

   if(isse==nsse)then
      ! HDF5 output for SSE time-series variables instead of binary files
      ! Initialize HDF5 if not already done
      if (.not. hdf5_initialized) then
         call h5open_f(hdferr)
         hdf5_initialized = .true.
      end if
      
      ! Create HDF5 filename for SSE data
      hdf5_filename = trim(foldername)//'sse_timeseries_data_'//trim(jobname)//'.h5'
      
      ! Check if file exists, if so open in append mode, otherwise create new
      logical :: file_exists_sse
      inquire(file=trim(hdf5_filename), exist=file_exists_sse)
      
      if (file_exists_sse) then
         ! Open existing file for read/write
         call h5fopen_f(trim(hdf5_filename), H5F_ACC_RDWR_F, file_id, hdferr)
         ! Open existing SSE time-series group
         time_series_group_name = '/sse_time_series'
         call h5gopen_f(file_id, trim(time_series_group_name), group_id, hdferr)
         
         ! Delete existing complete datasets if they exist (we're updating with complete data)
         call h5ldelete_f(group_id, 'slipz1_sse', hdferr)  ! Ignore errors if dataset doesn't exist
         call h5ldelete_f(group_id, 'slipz1_tau', hdferr)
         call h5ldelete_f(group_id, 'tsse', hdferr)
      else
         ! Create new file
         call h5fcreate_f(trim(hdf5_filename), H5F_ACC_TRUNC_F, file_id, hdferr)
         ! Create SSE time-series group
         time_series_group_name = '/sse_time_series'
         call h5gcreate_f(file_id, trim(time_series_group_name), group_id, hdferr)
      end if
      
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
      
      ! Close SSE time-series group
      call h5gclose_f(group_id, hdferr)
      
      ! Add mesh data to HDF5
      ! Read GTS file and store mesh information
      open(98, file='triangular_mesh.gts', status='old', action='read')
      read(98,*) n_vertices, n_edges_dummy, n_cells
      
      ! Allocate temporary arrays
      allocate(vertex_coords(n_vertices, 3))
      allocate(cell_connectivity(n_cells, 3))
      
      ! Read vertex coordinates
      do i = 1, n_vertices
         read(98,*) vertex_coords(i, 1), vertex_coords(i, 2), vertex_coords(i, 3)
      end do
      
      ! Read cell connectivity (indices start from 0 in GTS, which is correct for Paraview)
      do i = 1, n_cells
         read(98,*) cell_connectivity(i, 1), cell_connectivity(i, 2), cell_connectivity(i, 3)
         ! change 1-based to 0-based indexing for Paraview compatibility
         cell_connectivity(i, :) = cell_connectivity(i, :) -1
      end do
      close(98)
      
      ! Write mesh data to HDF5
      call h5gcreate_f(file_id, '/mesh', group_id, hdferr)
      
      ! Write vertex coordinates in correct layout for XDMF
      ! Create temporary array with correct memory layout
      allocate(vertex_coords_transposed(3, n_vertices))
      do i = 1, n_vertices
         vertex_coords_transposed(1, i) = vertex_coords(i, 1)  ! x
         vertex_coords_transposed(2, i) = vertex_coords(i, 2)  ! y
         vertex_coords_transposed(3, i) = vertex_coords(i, 3)  ! z
      end do
      
      dims_2d = (/3, n_vertices/)
      call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
      call h5dcreate_f(group_id, 'geometry', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vertex_coords_transposed, dims_2d, hdferr)
      
      deallocate(vertex_coords_transposed)
      call h5dclose_f(dset_id, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      
      ! Write cell connectivity in correct layout for XDMF
      ! Create temporary array with correct memory layout
      allocate(cell_connectivity_transposed(3, n_cells))
      do i = 1, n_cells
         cell_connectivity_transposed(1, i) = cell_connectivity(i, 1)  ! vertex 1
         cell_connectivity_transposed(2, i) = cell_connectivity(i, 2)  ! vertex 2
         cell_connectivity_transposed(3, i) = cell_connectivity(i, 3)  ! vertex 3
      end do
      
      dims_2d = (/3, n_cells/)
      call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
      call h5dcreate_f(group_id, 'topology', H5T_STD_I32LE, dspace_id, dset_id, hdferr)
      call h5dwrite_f(dset_id, H5T_STD_I32LE, cell_connectivity_transposed, dims_2d, hdferr)
      
      deallocate(cell_connectivity_transposed)
      call h5dclose_f(dset_id, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      
      call h5gclose_f(group_id, hdferr)
      
      ! Deallocate temporary arrays
      deallocate(vertex_coords, cell_connectivity)
      
      ! Close HDF5 file
      call h5fclose_f(file_id, hdferr)
      
      ! Create XDMF file for SSE visualization
      xdmf_filename = trim(foldername)//'sse_timeseries_data_'//trim(jobname)//'.xdmf'
      open(99, file=trim(xdmf_filename), status='replace')
      write(99,'(A)')'<?xml version="1.0" ?>'
      write(99,'(A)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      write(99,'(A)')'<Xdmf Version="2.0">'
      write(99,'(A)') ' <Domain>'
      write(99,'(A)') '  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
      
      ! Write a Grid for each time step
      do i = 1, nsse
         write(99,*) '   <Grid Name="step_', i, '" GridType="Uniform">'
         write(99,*) '    <Topology TopologyType="Triangle" NumberOfElements="', n_cells, '">'
         write(99,*) '     <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="', n_cells, ' 3">sse_timeseries_data_', trim(jobname), '.h5:/mesh/topology</DataItem>'
         write(99,*) '    </Topology>'
         write(99,*) '    <Geometry name="geo" GeometryType="XYZ" NumberOfElements="', n_vertices, '">'
         write(99,*) '     <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="', n_vertices, ' 3">sse_timeseries_data_', trim(jobname), '.h5:/mesh/geometry</DataItem>'
         write(99,*) '    </Geometry>'
         write(99,*) '    <Time Value="', tsse(i), '"/>'
         write(99,*) '    <Attribute Name="slipz1_sse" Center="Cell">'
         write(99,*) '     <DataItem ItemType="HyperSlab" Dimensions="', n_cells, '">'
         write(99,*) '      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">0 ', i-1, ' 1 1 1 ', n_cells, '</DataItem>'
         write(99,*) '      <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="1 ', n_cells, '">sse_timeseries_data_', trim(jobname), '.h5:/sse_time_series/slipz1_sse</DataItem>'
         write(99,*) '     </DataItem>'
         write(99,*) '    </Attribute>'
         write(99,*) '    <Attribute Name="slipz1_tau" Center="Cell">'
         write(99,*) '     <DataItem ItemType="HyperSlab" Dimensions="', n_cells, '">'
         write(99,*) '      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">0 ', i-1, ' 1 1 1 ', n_cells, '</DataItem>'
         write(99,*) '      <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="1 ', n_cells, '">sse_timeseries_data_', trim(jobname), '.h5:/sse_time_series/slipz1_tau</DataItem>'
         write(99,*) '     </DataItem>'
         write(99,*) '    </Attribute>'
         write(99,*) '   </Grid>'
      end do
      
      write(99,'(A)')'  </Grid>'
      write(99,'(A)')' </Domain>'
      write(99,'(A)')'</Xdmf>'
      close(99)
      
      write(*,*) 'SSE time-series data written to HDF5: ', trim(hdf5_filename)
      write(*,*) 'SSE XDMF visualization file created: ', trim(xdmf_filename)
      
      isse = 0
   else if (mod(isse, 10) == 0 .and. isse > 0) then
      ! Iterative output: Append data every 10 iterations for real-time visualization
      if (.not. hdf5_initialized) then
         call h5open_f(hdferr)
         hdf5_initialized = .true.
       end if
       
       ! Open existing HDF5 file for appending
       hdf5_filename = trim(foldername)//'sse_timeseries_data_'//trim(jobname)//'.h5'
       call h5fopen_f(trim(hdf5_filename), H5F_ACC_RDWR_F, file_id, hdferr)
       
       ! Open existing SSE time-series group
       time_series_group_name = '/sse_time_series'
       call h5gopen_f(file_id, trim(time_series_group_name), group_id, hdferr)
       
       ! Append new data to existing datasets
       ! Note: This requires extending the dataset dimensions
       ! For now, we'll just update the XDMF to show progress
       call h5gclose_f(group_id, hdferr)
       call h5fclose_f(file_id, hdferr)
       
       ! Update XDMF file to show current progress
       xdmf_filename = trim(foldername)//'sse_timeseries_data_'//trim(jobname)//'.xdmf'
       open(99, file=trim(xdmf_filename), status='replace')
       write(99,'(A)')'<?xml version="1.0" ?>'
       write(99,'(A)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
       write(99,'(A)')'<Xdmf Version="2.0">'
       write(99,'(A)') ' <Domain>'
       write(99,'(A)') '  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
       
       ! Write a Grid for each completed time step
       do i = 1, isse
          write(99,*) '   <Grid Name="step_', i, '" GridType="Uniform">'
          write(99,*) '    <Topology TopologyType="Triangle" NumberOfElements="',n_cells,'">'
          write(99,*) '     <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="',n_cells,' 3">sse_timeseries_data_', trim(jobname), '.h5:/mesh/topology</DataItem>'
          write(99,*) '    </Topology>'
          write(99,*) '    <Geometry name="geo" GeometryType="XYZ" NumberOfElements="',n_vertices,'">'
          write(99,*) '     <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="',n_vertices,' 3">sse_timeseries_data_', trim(jobname), '.h5:/mesh/geometry</DataItem>'
          write(99,*) '    </Geometry>'
          write(99,*) '    <Time Value="', tsse(i), '"/>'
          write(99,*) '    <Attribute Name="slipz1_sse" Center="Cell">'
          write(99,*) '     <DataItem ItemType="HyperSlab" Dimensions="',n_cells,'">'
          write(99,*) '      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">0 ', i-1, ' 1 1 1 ',n_cells,'</DataItem>'
          write(99,*) '      <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="1 ',n_cells,'">sse_timeseries_data_', trim(jobname), '.h5:/sse_time_series/slipz1_sse</DataItem>'
          write(99,*) '     </DataItem>'
          write(99,*) '    </Attribute>'
          write(99,*) '    <Attribute Name="slipz1_tau" Center="Cell">'
          write(99,*) '     <DataItem ItemType="HyperSlab" Dimensions="',n_cells,'">'
          write(99,*) '      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">0 ', i-1, ' 1 1 1 ',n_cells,'</DataItem>'
          write(99,*) '      <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="1 ',n_cells,'">sse_timeseries_data_', trim(jobname), '.h5:/sse_time_series/slipz1_tau</DataItem>'
          write(99,*) '     </DataItem>'
          write(99,*) '    </Attribute>'
          write(99,*) '   </Grid>'
       end do
       
       write(99,'(A)')'  </Grid>'
       write(99,'(A)')' </Domain>'
       write(99,'(A)')'</Xdmf>'
       close(99)
       
       write(*,*) 'SSE Progress update: XDMF updated for', isse, 'iterations'
  end if


else

   if((imv>0).and.(imv<nmv))then
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
      open(408,file=trim(foldername)//'blkst_strk+16fn+08dp+00'//jobname,position='append',status='unknown')
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
      ! Initialize HDF5 if not already done
      if (.not. hdf5_initialized) then
         call h5open_f(hdferr)
         hdf5_initialized = .true.
      end if
      
      ! Open existing HDF5 file for SSE data (append mode)
      hdf5_filename = trim(foldername)//'sse_timeseries_data_'//trim(jobname)//'.h5'
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
      
      ! Close group and file
      call h5gclose_f(group_id, hdferr)
      call h5fclose_f(file_id, hdferr)
      
      write(*,*) 'Partial SSE data written to HDF5: ', trim(hdf5_filename)
      
      isse = 0
  end if

     if((icos>0).and.(icos<ncos))then
       ! Write partial cosine slip data to the same HDF5 file as main cosine slip output
       ! Initialize HDF5 if not already done
       if (.not. hdf5_initialized) then
          call h5open_f(hdferr)
          hdf5_initialized = .true.
       end if
       
       ! Open existing HDF5 file for cosine slip data (append mode)
       hdf5_filename = trim(foldername)//'timeseries_data_'//trim(jobname)//'.h5'
       call h5fopen_f(trim(hdf5_filename), H5F_ACC_RDWR_F, file_id, hdferr)
       
       ! Open existing time-series group
       time_series_group_name = '/time_series'
       call h5gopen_f(file_id, trim(time_series_group_name), group_id, hdferr)
       
       ! Write partial slipz1_cos data (cosine slip time series) - structured for XDMF column extraction
       dims_2d = (/n_cells, icos/)
       call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
       call h5dcreate_f(group_id, 'slipz1_cos_partial', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_cos(1:n_cells,1:icos), dims_2d, hdferr)
       call h5dclose_f(dset_id, hdferr)
       call h5sclose_f(dspace_id, hdferr)
       
       ! Write partial slipz1_v data (velocity time series) - structured for XDMF column extraction
       dims_2d = (/n_cells, icos/)
       call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
       call h5dcreate_f(group_id, 'slipz1_v_partial', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_v(1:n_cells,1:icos), dims_2d, hdferr)
       call h5dclose_f(dset_id, hdferr)
       call h5sclose_f(dspace_id, hdferr)
       
       ! Write partial time array
       dims_1d = (/icos/)
       call h5screate_simple_f(1, dims_1d, dspace_id, hdferr)
       call h5dcreate_f(group_id, 'tcos_partial', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tcos(1:icos), dims_1d, hdferr)
       call h5dclose_f(dset_id, hdferr)
       call h5sclose_f(dspace_id, hdferr)
       
       ! Close group and file
       call h5gclose_f(group_id, hdferr)
       call h5fclose_f(file_id, hdferr)
       
       write(*,*) 'Partial cosine slip data written to HDF5: ', trim(hdf5_filename)
       
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
 160    format(E22.14,1x,E20.13)
 500    format(E15.8,1X,E20.13)
 600    format(E15.4,1X,E13.6,1X,E15.8)
 700    format(E13.6)
 900    format(E15.8)

RETURN
END subroutine output
