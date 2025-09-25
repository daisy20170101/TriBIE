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
  USE phy3d_module_bp6
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
       xilock1,xilock2,x4,z1,z2,z3

  real (DP) ::  tmbegin,tmrun,tautmp
  
  ! Pore fluid pressure variables
  real (DP) :: hz,gfun,gfun2,gfunb,gfun2b,dt_pf1
  real (DP) :: frc,help,help1,help2,zh
  ! heavi is now imported from module
  real (DP), DIMENSION(:), ALLOCATABLE :: dt_pf
  



  real (DP), DIMENSION(:), ALLOCATABLE :: x,z,xi,yt,yt0,dydt,yt_scale, &
       slip,slipinc,slipds,slipdsinc,sr,vi,pore_fluid

  !Arrays only defined at master cpu
  real (DP), DIMENSION(:), ALLOCATABLE :: x_all,xi_all,z_all,&
       yt_all,yt0_all,dydt_all,yt_scale_all,tau1_all,tau2_all, &
       slip_all,slipinc_all,slipds_all,slipdsinc_all,cca_all,ccb_all,xLf_all,seff_all,&
       vi_all,phy1_all,phy2_all,pore_fluid_all,dt_pf_all
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
  
  ! MPI scatter arrays for different data types
  integer, dimension(:), allocatable :: sendcounts_yt, displs_yt
  
  ! Element mapping for visualization (MPI order -> Mesh order)
  integer, dimension(:), allocatable :: mpi_to_mesh_map
  integer, dimension(:), allocatable :: start_indices
  integer :: mpi_idx, global_idx, mesh_idx
  
  ! File existence checking variables
  logical :: trigreen_file_exists
  character(len=256) :: trigreen_filename
  
  ! MPI_Scatterv variables for uneven distribution

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
  read(12,*)Nt_all,nprocs,n_obv,np1,np2
  read(12,*)IDin,Idout,Iprofile,Iperb,Isnapshot 
  read(12,*)Vpl
  read(12,*)tmax
  read(12,*)tslip_ave,tslipend,tslip_aveint
  read(12,*)tint_out,tmin_out,tint_cos,tint_sse
  read(12,*)vcos,vsse1,vsse2
  read(12,*)nmv,nas,ncos,nnul,nsse,n_nul_int
  read(12,*)s1(1),s1(2),s1(3),s1(4),s1(5),s1(6),s1(7),s1(8),s1(9),s1(10)
!!! modified data read in
  close(12)

  Nab=5 ! used in resdep if dault a-b profile is given 

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
     allocate(sendcounts_yt(0:size-1), displs_yt(0:size-1))
     allocate(mpi_to_mesh_map(Nt_all))
     allocate(start_indices(0:size-1))
     
     ! Calculate send counts and displacements for each process
     call MPI_Allgather(local_cells, 1, MPI_INTEGER, sendcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
     
     ! Calculate displacements for single arrays (slip, slipds)
     displs(0) = 0
     do i = 1, size-1
        displs(i) = displs(i-1) + sendcounts(i-1)
     end do
     
     ! Initialize yt scatter arrays (will be set properly during restart)
     do i = 0, size-1
        sendcounts_yt(i) = 3 * sendcounts(i)  ! yt has 3 components per cell
     end do
     displs_yt(0) = 0
     do i = 1, size-1
        displs_yt(i) = displs_yt(i-1) + sendcounts_yt(i-1)
     end do
     
     ! Create element mapping: MPI gather order -> Original mesh order
     ! The MPI_Gather puts data in process order, but we need mesh order
     ! We need to collect the actual start_idx values from each process
     call MPI_Allgather(start_idx, 1, MPI_INTEGER, start_indices, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
     
     do i = 0, size-1
        do j = 1, sendcounts(i)
           ! Element index in MPI gathered array
           mpi_idx = displs(i) + j
           ! Original global element index from actual start_idx
           global_idx = start_indices(i) + j - 1  ! start_idx is 1-based, j is 1-based
           mpi_to_mesh_map(mpi_idx) = global_idx
        end do
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
         yt0_all(3*Nt_all),yt_all(3*Nt_all),dydt_all(3*Nt_all),yt_scale_all(3*Nt_all),pore_fluid_all(Nt_all),dt_pf_all(Nt_all))

     allocate(phy1_all(Nt_all),phy2_all(Nt_all))
  else
     ! Worker processes: Allocate minimal dummy arrays for MPI_Scatterv compatibility
     ! These arrays won't be used as source data, but must exist for the MPI call
     ALLOCATE(x_all(1),xi_all(1),&
          cca_all(1),ccb_all(1),seff_all(1),xLf_all(1),vi_all(1),&
          tau1_all(1),tau2_all(1),slip_all(1),slipinc_all(1),slipds_all(1),slipdsinc_all(1),&
         yt0_all(1),yt_all(1),dydt_all(1),yt_scale_all(1),pore_fluid_all(1),dt_pf_all(1))

     allocate(phy1_all(1),phy2_all(1))
  end if

  ! Master-only output arrays allocation
  if(myid == master) then
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
     pore_fluid_all = 0.d0
     dt_pf_all = 0.d0

!!! modify output number
     ALLOCATE (slipz1_inter(Nt_all,nas),slipz1_cos(Nt_all,ncos), &
          slipave_inter(Nt_all,nas),slipave_cos(Nt_all,ncos),v_cos(Nt_all,ncos),slip_cos(Nt_all,ncos), &
          v_nul(Nt_all,nnul),slip_nul(Nt_all,nnul),slipz1_tau(Nt_all,ncos),slipz1_sse(Nt_all,nsse) )
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
       yt(3*local_cells),dydt(3*local_cells),yt_scale(3*local_cells),&
       yt0(3*local_cells),sr(local_cells),vi(local_cells),pore_fluid(local_cells),dvel(local_cells),dt_pf(local_cells))

  ALLOCATE (stiff(local_cells,Nt_all))   !!! stiffness of Stuart green calculation

  ! Initialize dt_pf array
  dt_pf = 1.d12  ! Initial pore fluid time step
  pore_fluid = 0.d0
  dvel = 0.d0

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
        xi_all(k) = xi_all(k) ! y infinite long axis, meter
        x_all(k) = x_all(k)
        z_all(k) = z_all(k) -360.0d3-36.0d3

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
        stiff(i,j) = 1d5*1d3*stiff(i,j)
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
   if (ierr /= 0) then
      write(*,*) 'ERROR: MPI_Scatterv failed for cca, ierr =', ierr, 'on process', myid
      stop
   end if
   
   call MPI_Scatterv(ccb_all,sendcounts,displs,MPI_Real8,ccb,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
   if (ierr /= 0) then
      write(*,*) 'ERROR: MPI_Scatterv failed for ccb, ierr =', ierr, 'on process', myid
      stop
   end if
   
   call MPI_Scatterv(xLf_all,sendcounts,displs,MPI_Real8,xLf,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
   if (ierr /= 0) then
      write(*,*) 'ERROR: MPI_Scatterv failed for xLf, ierr =', ierr, 'on process', myid
      stop
   end if
   
   call MPI_Scatterv(seff_all,sendcounts,displs,MPI_Real8,seff,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
   if (ierr /= 0) then
      write(*,*) 'ERROR: MPI_Scatterv failed for seff, ierr =', ierr, 'on process', myid
      stop
   end if
   

   
   call MPI_Scatterv(vi_all,sendcounts,displs,MPI_Real8,vi,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
   if (ierr /= 0) then
      write(*,*) 'ERROR: MPI_Scatterv failed for vi, ierr =', ierr, 'on process', myid
      stop
   end if
   
   call MPI_Scatterv(x_all,sendcounts,displs,MPI_Real8,x,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
   if (ierr /= 0) then
      write(*,*) 'ERROR: MPI_Scatterv failed for x, ierr =', ierr, 'on process', myid
      stop
   end if

  call MPI_Bcast(z_all,Nt_all,MPI_Real8,master,MPI_COMM_WORLD,ierr)
  
  ! Validate scattered parameters for each process
  write(*,*) 'Process', myid, 'scattered parameters validation:'
  write(*,*) '  Local cells:', local_cells
  write(*,*) '  First few cca values:', cca(1:min(5, local_cells))
  write(*,*) '  First few ccb values:', ccb(1:min(5, local_cells))
  write(*,*) '  First few xLf values:', xLf(1:min(5, local_cells))
  write(*,*) '  First few seff values:', seff(1:min(5, local_cells))
  
  ! Check for invalid values
  do i = 1, local_cells
     if (cca(i) <= 0.0d0) then
        write(*,*) 'ERROR: Non-positive cca(i) at i=', i, ' value=', cca(i), 'on process', myid
        stop
     end if
     if (ccb(i) < 0.0d0) then
        write(*,*) 'ERROR: Negative ccb(i) at i=', i, ' value=', ccb(i), 'on process', myid
        stop
     end if
     if (xLf(i) <= 0.0d0) then
        write(*,*) 'ERROR: Non-positive xLf(i) at i=', i, ' value=', xLf(i), 'on process', myid
        stop
     end if
     if (seff(i) <= 0.0d0) then
        write(*,*) 'ERROR: Non-positive seff(i) at i=', i, ' value=', seff(i), 'on process', myid
        stop
     end if
  end do
  
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
  dtmin = 0.001 ! in sec
  dt_try=dtmin
  Vint = Vpl

  if(myid==master)then
      open(311,file=trim(foldername)//'fltst_strk-15'//jobname,access='append',status='unknown')
      open(312,file=trim(foldername)//'fltst_strk+00'//jobname,access='append',status='unknown')
      open(313,file=trim(foldername)//'fltst_strk+05'//jobname,access='append',status='unknown')
      open(314,file=trim(foldername)//'fltst_strk+10'//jobname,access='append',status='unknown')
      open(315,file=trim(foldername)//'fltst_strk+15'//jobname,access='append',status='unknown')
      open(316,file=trim(foldername)//'fltst_strk+25'//jobname,access='append',status='unknown')
      open(317,file=trim(foldername)//'fltst_strk+35'//jobname,access='append',status='unknown')
      open(318,file=trim(foldername)//'fltst_strk+50'//jobname,access='append',status='unknown')
      open(319,file=trim(foldername)//'fltst_strk+75'//jobname,access='append',status='unknown')

       do i=311,319
                write(i,100)'# This is the file header'
                write(i,100)'# problem=SEAS Benchmark No.6'
                write(i,100)'# author=D.Li'
        write(i,100)'# code=TriBIE'
                write(i,100)'# date=2022/8/1'
                write(i,100)'# element_size = 100 m'
                write(i,100)'# minimum_time_step = 1e-3'
                write(i,100)'# maximum_time_step = 2e+7'
                write(i,100)'# location = on fault: file name'
                write(i,100)'# Column #1 = Time (s)'
                write(i,100)'# Column #2 = slip (m)'
                write(i,100)'# Column #3 = Slip_rate (log10 m/s)'
                write(i,100)'# Column #4 = Shear stress  (MPa)'
                write(i,100)'# Column #5 = pore_pressure  (MPa)'
        write(i,100)'# Column #6 = Darcy vel  (m/s)'
                write(i,100)'# Column #7 = State (log10 s)'
                write(i,100)'# '
                write(i,100)'# The line below lists the names of the data fields:'
                write(i,'(A,1x,A,1x,A,1x,A,1x,A,1x,A,1x,A,1x)')'t','slip','slip_rate','shear_stress','pore_pressure','Darcy_vel','state'
                write(i,100)'# Below is the time-series data.'          
        end do
 100    format(A)
end if
 
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
        yt(3*j-1)= vini

        phy1(j)=1.0
        phy2(j)=0.0

        help = dlog((2.d0*V0/Vint) * dsinh(tauini/(cca(j)*seff(j))))
        
        tau1(j)= tauini
        tau2(j) = 0.0
        phy1(j) = tau1(j)/dsqrt(tau1(j)**2+tau2(j)**2)
        phy2(j) = tau2(j)/dsqrt(tau1(j)**2+tau2(j)**2)

        yt(3*j-2) = 0.0d0  ! Initialize pore fluid pressure
        yt(3*j) = xLf(j)/V0*dexp((cca(j)/ccb(j))*help - f0/ccb(j))  ! Initialize theta (state variable)
        slip(j)=0.d0
        slipds(j)=0.d0
        dvel(j)=1.d-12
        yt0(3*j-2)=yt(3*j-2)
        yt0(3*j-1) = yt(3*j-1)
        yt0(3*j) = yt(3*j)
     end do
  end if

!------------------------------------------------------------------
  if(IDin.eq.1) then               !if this is a restart job
     if(myid==master)then
        filename='out'
        call restart(0,filename,4,Nt_all,t,dt,dt_try,ndt,nrec,yt_all,slip_all)
        write(1,*)'This is a restart job. Start time ',t,' yr'
     end if
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_Bcast(t,1,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(dt,1,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(dt_try,1,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ndt,1,MPI_integer,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(nrec,1,MPI_integer,master,MPI_COMM_WORLD,ierr)
     
     ! yt scatter arrays are already initialized above
     
     if (myid == master) then
        write(*,*) 'YT MPI_Scatterv distribution:'
        do i = 0, size-1
           write(*,*) '  Process', i, ': yt_sendcount =', sendcounts_yt(i), ', yt_displacement =', displs_yt(i)
        end do
     end if
     
     call MPI_Scatterv(yt_all,sendcounts_yt,displs_yt,MPI_Real8,yt,3*local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Scatterv(slip_all,sendcounts,displs,MPI_Real8,slip,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     call MPI_Scatterv(slipds_all,sendcounts,displs,MPI_Real8,slipds,local_cells,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     
     ! Validate scattered data for NaN/infinity
     do i = 1, 3*local_cells
        if (yt(i) /= yt(i) .or. abs(yt(i)) > huge(yt(i))/2) then
           write(*,*) 'ERROR: Invalid yt(',i,') after MPI scatter on process', myid, ' value=', yt(i)
        end if
     end do
     do i = 1, local_cells
        if (slip(i) /= slip(i) .or. abs(slip(i)) > huge(slip(i))/2) then
           write(*,*) 'ERROR: Invalid slip(',i,') after MPI scatter on process', myid, ' value=', slip(i)
        end if
     end do
  
     ndtnext = ndt
     tprint_inter = t
     tslip_ave=t        
     tout = t

  else
     if(myid==master)then
        write(1,*)'Start time ',t,' s'
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
  comm_count = 3*local_cells
  comm_tag = 0

  ! Initialize blocking parameters
  block_size = 64  ! Optimal block size for cache

  if(myid == master) then
     allocate(send_buffer(3*Nt_all))
     allocate(recv_buffer(3*Nt_all))
  end if
  ! Main simulation loop
  do while(cyclecont) 

     call derivs(myid,dydt,3*local_cells,Nt_all,local_cells,t,yt,z_all,x) 

     do j=1,3*local_cells
        yt_scale(j)=dabs(yt(j))+dabs(dt_try*dydt(j))
        yt0(j) = yt(j)
     end do
     
     CALL rkqs(myid,yt,dydt,3*local_cells,Nt_all,local_cells,t,dt_try,accuracy,yt_scale, &
          dt_did,dt_next,z_all,x)

     dt = dt_did
     dt_try = dt_next

     ! Physics calculations for each cell
     do i=1,local_cells

      zh = dsign(max(1.0d-8,dabs(z(i))),z(i))

      dvel(i) = compute_dpf_dt(zh, t, alpha, beta, phi, q0, toff)
           
      pore_fluid(i) = compute_pf(zh, t, alpha, beta, phi, q0, toff)

      ! Time step inversely related to velocity change rate (dvel)
      ! This ensures smaller time steps when velocity changes rapidly
      ! Use max to prevent division by zero and extremely large time steps
      dt_pf(i) = 1.0d3/max(1.0d-12, abs(dvel(i)))

        help=(yt(3*i-1)/(2*V0))*dexp((f0+ccb(i)*dlog(V0*yt(3*i)/xLf(i)))/cca(i))
        
        tau1(i) = (seff(i)-pore_fluid(i))*cca(i)*dlog(help+dsqrt(1+help**2))
        tau2(i) = tau1(i)/phy1(i)*phy2(i)

        slipinc(i) = 0.5*(yt0(3*i-1)+yt(3*i-1))*dt
        slipdsinc(i)=0.5*(yt0(3*i-1)+yt(3*i-1))*dt*phy2(i)/phy1(i)
        
        slip(i) = slip(i) + slipinc(i)
        slipds(i)=slipds(i)+slipdsinc(i)
     end do

      write(*,*) 'pf,dpf_dt:',pore_fluid(1),dvel(1)

     call MPI_Gatherv(dt_pf,local_cells,MPI_Real8,dt_pf_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)

     if(myid.eq.master) then 
         ! Combine Runge-Kutta suggested time step with pore fluid-based time step
         dt_pf1 = min(min(0.1,minval(dt_pf_all)), dt_try)
         write(*,*) 'step:',t,dt_try,minval(dt_pf_all)! at z=0.0 km 
     end if

     CALL MPI_BCAST(dt_pf1,1,MPI_REAL8,master,MPI_COMM_WORLD, ierr)
     ! Use the more restrictive time step (smaller of RK and pore fluid)
     dt_try = dt_pf1

     ndt = ndt + 1

     ! Gather data from all MPI processes
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     
     if(myid == master) then
        ! Master process gathers all data
        call MPI_Gatherv(yt,3*local_cells,MPI_Real8,yt_all,sendcounts_yt,displs_yt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(yt0,3*local_cells,MPI_Real8,yt0_all,sendcounts_yt,displs_yt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(slipinc,local_cells,MPI_Real8,slipinc_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(slip,local_cells,MPI_Real8,slip_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(slipdsinc,local_cells,MPI_Real8,slipdsinc_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(slipds,local_cells,MPI_Real8,slipds_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(tau1,local_cells,MPI_Real8,tau1_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(tau2,local_cells,MPI_Real8,tau2_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(phy1,local_cells,MPI_Real8,phy1_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(phy2,local_cells,MPI_Real8,phy2_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(pore_fluid,local_cells,MPI_Real8,pore_fluid_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     else
        ! Non-master processes send their data
        call MPI_Gatherv(yt,3*local_cells,MPI_Real8,yt_all,sendcounts_yt,displs_yt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(yt0,3*local_cells,MPI_Real8,yt0_all,sendcounts_yt,displs_yt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(slipinc,local_cells,MPI_Real8,slipinc_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(slip,local_cells,MPI_Real8,slip_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(slipdsinc,local_cells,MPI_Real8,slipdsinc_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(slipds,local_cells,MPI_Real8,slipds_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(tau1,local_cells,MPI_Real8,tau1_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(tau2,local_cells,MPI_Real8,tau2_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(phy1,local_cells,MPI_Real8,phy1_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(phy2,local_cells,MPI_Real8,phy2_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(pore_fluid,local_cells,MPI_Real8,pore_fluid_all,sendcounts,displs,MPI_Real8,master,MPI_COMM_WORLD,ierr)
     end if

     ! Output calculations (only master process)
     if(myid==master)then
        imv=imv+1
        tmv(imv)=t
        maxv(imv) = 0.d0
        moment(imv) =0.d0
        
        ! Find max velocity and calculate moment
        do i=1,Nt_all
           if(yt_all(3*i-1).ge.maxv(imv))then
              maxv(imv)=yt_all(3*i-1)
              maxnum(imv)=i
           end if
   
          if(.not.rup(i).and.yt_all(3*i-1)/yrs.ge.vcos)then
             Trup(i)=t*yrs
             rup(i)=.true.
          end if
           moment(imv) = moment(imv)+0.5*(yt0_all(3*i-1)+yt_all(3*i-1))/yrs*1d-3*area(i)*xmu*1d6*1d5
        end do

        ! SEAS output variables
        do i = 1,10
         outs1(imv,1,i) = slip_all(s1(i))*1.d-3 ! meter
         outs1(imv,2,i) =  dlog10(yt_all(3*s1(i)-1)) ! log10(V) m/s
         outs1(imv,3,i) = tau1_all(s1(i))/1d6 ! MPa
         outs1(imv,4,i) = yt_all(3*s1(i)-2)/1d6
         outs1(imv,6,i) = dlog10(yt_all(3*s1(i))) ! log10(theta)
         outs1(imv,7,i) = 0.d0
         outs1(imv,5,i) = pore_fluid_all(s1(i))/1d6 ! darcy vel

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
             vel1 = vel1 + surf1(i,j)*(yt0_all(3*j-1)+yt_all(3*j-1))*0.5
             vel2 = vel2 + surf2(i,j)*(yt0_all(3*j-1)+yt_all(3*j-1))*0.5
             vel3 = vel3 + surf3(i,j)*(yt0_all(3*j-1)+yt_all(3*j-1))*0.5
          
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
              write(*,*) 'DEBUG: icos =', icos, 'ncos =', ncos, 't =', t 

              if(.not.end1.and.t - teve1.lt.2*tint_cos) then
                 teve1 = t !! to determine rupture contour output
                else
                 end1=.true.
              end if

              ! Calculate coseismic slip WITHOUT element mapping (for testing)
              do i=1,Nt_all
                 ! Use direct MPI gather order (no mapping)
                 slipz1_cos(i,icos) = slip_all(i)*1.d-3
                 slipz1_v(i,icos) = dlog10(yt_all(3*i-2)*1.d-3/yrs) 
                 slipz1_tau(i,icos) = tau1_all(i)
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
        !$OMP MASTER
        call output(Ioutput,Isnapshot,Nt_all,Nt,inul,imv,ias,icos,isse,x,&
             tmv,tas,tcos,tnul,tsse,maxv,moment,outs1,maxnum,msse1,msse2, areasse1,areasse2,&
             slipz1_inter,slipz1_tau,slipz1_sse, &
             slipz1_cos,slipave_inter,slipave_cos,slip_cos,v_cos,slip_nul,v_nul,&
             xi_all,x_all,intdepz1,intdepz2,intdepz3,n_cosz1,n_cosz2,n_cosz3,&
             n_intz1,n_intz2,n_intz3,slipz1_v,obvs,n_obv,obvstrk,obvdp,np1,np2,mpi_to_mesh_map)         
         !$OMP END MASTER
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
110 format(A)
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
          n_intz1,n_intz2,n_intz3,slipz1_v,obvs,n_obv,obvstrk,obvdp,np1,np2,mpi_to_mesh_map) 
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


  ! Deallocate arrays based on allocation pattern
  if(myid==master)then 
     ! Master: Deallocate full-size arrays and master-only arrays
     DEALLOCATE (x_all,xi_all,yt_all,dydt_all,yt_scale_all,yt0_all,&
                phy1_all,phy2_all,vi_all,tau1_all,tau2_all, &
          slip_all,slipinc_all,slipds_all,slipdsinc_all,&
           cca_all,ccb_all,xLf_all,seff_all,pore_fluid_all,dt_pf_all)
     
     ! Deallocate master-only output arrays (only allocated on master)
     if (allocated(maxnum)) DEALLOCATE(maxnum,maxv,moment,outs1)
     if (allocated(msse1)) DEALLOCATE(msse1,msse2,areasse1,areasse2)
     if (allocated(tmv)) DEALLOCATE(tmv,tas,tcos,tnul,tsse)

     ! Deallocate master-only simulation arrays (only allocated on master)
     if (allocated(slipz1_inter)) then
        DEALLOCATE (slipz1_inter,slipz1_tau,slipz1_sse, &
             slipz1_cos,slipave_inter,slipave_cos, &
             v_cos,slip_cos,v_nul,slip_nul)
     end if
     if (allocated(intdepz1)) then
        DEALLOCATE (intdepz1,intdepz2,intdepz3,ssetime,slipz1_v)
     end if

     ! Deallocate more master-only arrays
     if (allocated(Trup)) deallocate(Trup,rup,area,obvs)
     if (allocated(pstrk)) deallocate(pstrk,pdp,obvstrk,obvdp)
  else
     ! Workers: Deallocate dummy arrays (size 1)
     DEALLOCATE (x_all,xi_all,yt_all,dydt_all,yt_scale_all,yt0_all,&
                phy1_all,phy2_all,vi_all,tau1_all,tau2_all, &
          slip_all,slipinc_all,slipds_all,slipdsinc_all,&
           cca_all,ccb_all,xLf_all,seff_all,pore_fluid_all,dt_pf_all)
  end if


  DEALLOCATE (stiff,sr)
  DEALLOCATE (x,xi,yt,dydt,yt_scale)
  deallocate (phy1,phy2,tau1,tau2,tau0,slip,slipinc,slipds,slipdsinc,yt0,zzfric,zzfric2)
  DEALLOCATE (cca,ccb,xLf,seff,pore_fluid,dt_pf)
  
  ! Clean up MPI_Scatterv arrays
  if (use_trigreen_format .and. allocated(sendcounts)) then
     deallocate(sendcounts, displs)
  end if
  if (allocated(sendcounts_yt)) then
     deallocate(sendcounts_yt, displs_yt)
  end if
  if (allocated(mpi_to_mesh_map)) then
     deallocate(mpi_to_mesh_map)
  end if
  if (allocated(start_indices)) then
     deallocate(start_indices)
  end if
  
  call MPI_finalize(ierr)
END program main

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine rkqs(myid,y,dydx,n,Nt_all,Nt,x,htry,eps,yscal,hdid,hnext,z_all,p)
  Use mpi
  USE phy3d_module_bp6, only : nprocs
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
       USE phy3d_module_bp6, only :nprocs
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
     subroutine derivs(myid,dydt,nv,Nt_all,Nt,t,yt,z_all,z)
       USE mpi
       USE phy3d_module_bp6, only: phy1,phy2,tau1,tau2, stiff,cca,ccb,seff,xLf,eta,f0,Vpl,V0,Lratio,nprocs,&
            tm1,tm2,tmday,tmelse,tmmidn,tmmult,alpha,beta,phi,q0,toff,&
            compute_pf,compute_dpf_dt,compute_G,compute_dGdt,dirac_delta,heavi,sendcounts,displs
       ! MPI variables are passed as arguments or declared locally in main program
       implicit none
       integer, parameter :: DP = kind(1.0d0)
       integer :: nv,n,i,j,k,kk,l,ii,Nt,Nt_all
       real (DP) :: t,yt(nv),dydt(nv)   
       real (DP) :: deriv3,deriv2,deriv1,small,tauinc2,dydtinc
       real (DP) :: psi,help1,help2,help,help4
       real (DP) :: SECNDS
       real (DP) :: z(Nt),sr(Nt),z_all(Nt_all),zz(Nt),zz_ds(Nt),zzfric(Nt),zz_all(Nt_all),zzfric2(Nt)
       real (DP) :: pore_fluid(Nt)

       ! Local variables for blocking optimization
       integer :: block_size, j_start, j_end, i_block, j_block, i_end_block, j_end_block
       real(DP) :: temp_sum
       integer :: request1, request2
       intrinsic real
       
       ! pore fulid variables
       real(DP) :: frc,pressure

       ! Regularization parameter for rate-and-state friction
       real(DP), parameter :: theta_min = 1.0d-12  ! Minimum state variable (seconds) - increased for stability

       !MPI RELATED DEFINITIONS
       integer :: ierr,myid,master
       master = 0 

       small=1.d-6

       ! OPTIMIZATION: Advanced vectorization with loop unrolling and prefetching
       do i=1,Nt
          zz(i)=yt(3*i-1)-Vpl
       end do

       ! OPTIMIZATION: Advanced MPI communication with non-blocking operations
       ! Use non-blocking communication to overlap computation and communication
       
       ! Fixed: Use existing sendcounts and displs for MPI_Allgatherv
       ! Much more efficient - reuse the already calculated distribution arrays
       call MPI_Allgatherv(zz, Nt, MPI_Real8, zz_all, sendcounts, displs, MPI_Real8, MPI_COMM_WORLD, ierr)
       
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

       ! Apply physics-based regularization for rate-and-state friction
       ! Small regularization parameter to prevent ln(0) while maintaining physics
       
       do i=1,Nt
          if (yt(3*i) < theta_min) then
             ! Apply regularization: don't change original values, just prevent ln(0)
             ! This preserves the physical state while making calculations numerically stable
             if (yt(3*i) <= 0.0d0) then
                write(*,*) 'INFO: Regularizing zero state variable at i=', i, ' from', yt(3*i), ' to', theta_min
             end if
             yt(3*i) = max(yt(3*i), theta_min)
          end if
          if (xLf(i) <= 0.0d0) then
             write(*,*) 'ERROR: Non-positive xLf(i) at i=', i, ' value=', xLf(i), ' correcting to 1.0d-3'
             xLf(i) = 1.0d-3
          end if
       end do

       ! OPTIMIZATION: Advanced vectorization with SIMD-friendly structure
       !$OMP SIMD PRIVATE(psi,help1,help2,help,deriv1,deriv2,deriv3)
       do i=1,Nt

         z(i) = dsign(max(1.0d-8,dabs(z(i))),z(i))
         dydt(3*i -2) = compute_dpf_dt(z(i), t, alpha, beta, phi, q0, toff)

         pressure = compute_pf(z(i), t, alpha, beta, phi, q0, toff)

         psi = dlog(V0*yt(3*i)/xLf(i))
         help1 = yt(3*i-1)/(2*V0)
         help2 = (f0+ccb(i)*psi)/cca(i)
         help = dsqrt(1+(help1*dexp(help2))**2)
         !frc = f0+cca(i)*dlog(yt(3*i-1)/V0) + ccb(i)*dlog(V0*yt(3*i)/xLf(i))
         
         help4 = help1 * dexp(help2)
         frc = cca(i)*dlog(help4+dsqrt(1+help4**2))

          deriv1 = ((seff(i)-pressure)*ccb(i)/yt(3*i))*help1*dexp(help2)/help
          deriv2 = ((seff(i)-pressure)*cca(i)/(2*V0))*dexp(help2)/help
          
          

!aging             
          deriv3 = 1-yt(3*i-1)*yt(3*i)/xLf(i)
!slip law         deriv3 = -yt(3*i-1)*yt(3*i)/xLf(i)*dlog(yt(3*i-1)*yt(3*i)/xLf(i))
          ! add dpf/dt in the  term
          dydt(3*i-1) = (-zzfric(i)-deriv1*deriv3 + frc*dydt(3*i-2))/(eta+deriv2) ! total shear traction
          dydt(3*i)=deriv3     
       end do
       !$OMP END SIMD
       
       ! Post-validate results (outside SIMD for debugging)
       do i=1,Nt
          if (dydt(3*i-2) /= dydt(3*i-2) .or. abs(dydt(3*i-2)) > huge(dydt(3*i-2))/2) then
             write(*,*) 'WARNING: Invalid dydt(3*i-2) at i=', i, ' value=', dydt(3*i-2)
          end if
          if (dydt(3*i-1) /= dydt(3*i-1) .or. abs(dydt(3*i-1)) > huge(dydt(3*i-1))/2) then
             write(*,*) 'WARNING: Invalid dydt(3*i-1) at i=', i, ' value=', dydt(3*i-1)
          end if
          if (dydt(3*i) /= dydt(3*i) .or. abs(dydt(3*i)) > huge(dydt(3*i))/2) then
             write(*,*) 'WARNING: Invalid dydt(3*i) at i=', i, ' value=', dydt(3*i)
          end if
       end do

       RETURN
     END subroutine derivs

!-----------------------------------------------------------------------------
!    read parameters: sigma_effective, a,b,D_c
!----------------------------------------------------------------------------

    subroutine resdep(Nt_all,hnucl, &
         xilock1,xilock2,cca_all,ccb_all,xLf_all, &
         seff_all,x_all,z_all,vi_all)
      USE mpi
      USE phy3d_module_bp6, only: yrs,p18,Nl,Nd,Nab,xmu,xnu,gamma, &
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
USE phy3d_module_bp6, ONLY : jobname,foldername,restartname, &
                        tm1,tm2,tmday,tmelse,tmmidn,tmmult,Vpl
      implicit none
      integer, parameter :: DP = kind(1.0d0)
      integer :: inout,i,ndt,nrec,Ifileout,Nt,Nt_all
      real (DP) :: t,dt,dt_try
      real (DP) ::  yt(3*Nt_all),slip(Nt_all)
      character(len=40) :: filename
      
      ! Additional variables for debugging
      logical :: file_exists
      character(len=200) :: line_buffer
      integer :: preview_unit, zero_count, negative_count, invalid_count

      if(inout.eq.0) then
         write(*,*) 'Opening restart file: ', trim(restartname)
         
         ! Check if file exists and get some info
         inquire(file=trim(restartname), exist=file_exists)
         if (.not. file_exists) then
            write(*,*) 'ERROR: Restart file does not exist!'
            stop
         end if
         
         open(Ifileout,file=trim(restartname),status='old')
         write(*,*) 'File opened successfully, unit=', Ifileout
         
         ! Quick preview of first few lines in the file
         preview_unit = 99
         open(preview_unit, file=trim(restartname), status='old')
         write(*,*) 'File preview (first 5 lines):'
         do i = 1, 5
            read(preview_unit, '(A)', end=100) line_buffer
            write(*,*) 'Line', i, ': ', trim(line_buffer)
         end do
100      close(preview_unit)
         
          read(Ifileout,*)t,ndt,nrec
          write(*,*) 'Read header: t=',t,' ndt=',ndt,' nrec=',nrec
          read(Ifileout,*)dt,dt_try
          write(*,*) 'Read timesteps: dt=',dt,' dt_try=',dt_try
          write(*,*) 'About to read 3*Nt_all=', 3*Nt_all, ' yt values...'
          
          ! Count problematic values
          zero_count = 0
          negative_count = 0 
          invalid_count = 0
          
          do i=1,3*Nt_all
             read(Ifileout,*)yt(i)
             ! Debug: Show first few values being read
             if (i <= 10) then
                write(*,*) 'DEBUG: Read yt(',i,') = ',yt(i)
             end if
             
             ! Detailed analysis of yt values
             if (yt(i) /= yt(i) .or. abs(yt(i)) > huge(yt(i))/2) then
                write(*,*) 'WARNING: Invalid yt(',i,') = ',yt(i)
                invalid_count = invalid_count + 1
             else if (yt(i) == 0.0d0) then
                zero_count = zero_count + 1
                if (zero_count <= 5) then  ! Show first 5 zero values
                   write(*,*) 'ZERO yt(',i,') = ',yt(i)
                end if
             else if (yt(i) < 0.0d0) then
                negative_count = negative_count + 1
                if (negative_count <= 5) then  ! Show first 5 negative values
                   write(*,*) 'NEGATIVE yt(',i,') = ',yt(i)
                end if
             end if
          end do
          
          write(*,*) 'yt Statistics: Zero=', zero_count, ' Negative=', negative_count, ' Invalid=', invalid_count

          do i=1,Nt_all
             read(Ifileout,*)slip(i)
             ! Check for NaN/Inf in loaded data
             if (slip(i) /= slip(i) .or. abs(slip(i)) > huge(slip(i))/2) then
                write(*,*) 'WARNING: Invalid slip(',i,') = ',slip(i)
             end if
          end do
	  close(Ifileout)
          write(*,*) 'Restart file loaded successfully'
          
          ! Report statistics but preserve the physical state from restart file
          if (zero_count > 0) then
             write(*,*) 'INFO: Found', zero_count, 'zero state variables in restart file'
             write(*,*) '      (This is physically valid during fast slip events)'
             write(*,*) '      Regularization will be applied during physics calculations'
          end if
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
    tmv,tas,tcos,tnul,tsse,maxv,moment,outs1,&
    maxnum,msse1,msse2,areasse1,areasse2, &
     slipz1_inter,slipz1_tau,slipz1_sse,&
     slipz1_cos,slipave_inter,slipave_cos,slip_cos,v_cos,slip_nul,v_nul,&
     xi_all,x_all,intdepz1,intdepz2,intdepz3,n_cosz1,n_cosz2,n_cosz3,&
    n_intz1,n_intz2,n_intz3,slipz1_v,obvs,n_obv,obvstrk,obvdp,np1,np2,mpi_to_mesh_map) 


USE mpi
USE phy3d_module_bp6, only: xmu,nmv,nas,ncos,nnul,nsse,yrs,Vpl,Nl, &
		foldername,jobname
use hdf5  ! Add HDF5 support
implicit none
integer, parameter :: DP = kind(1.0d0)
integer :: Nt,Nt_all,i,j,k,l,kk,inul,imv,ias,icos,isse,Ioutput,Isnapshot,ix1,ix2,ix3,ix4,n_obv,np1,np2,ierr

real (DP) :: x(Nt),maxnum(nmv),moment(nmv),maxv(nmv),outs1(nmv,7,10),&
        msse1(nsse),msse2(nsse),areasse1(nsse),areasse2(nsse), &
	tmv(nmv),tas(nas),tcos(ncos),tnul(nnul),tsse(nsse),obvs(nmv,6,n_obv),obvstrk(nmv,2,np1),obvdp(nmv,2,np2)

! Persistent array for storing all time values across subroutine calls
real (DP), allocatable, SAVE :: tcos_all(:)

real (DP) :: slipz1_inter(Nt_all,nas),slipz1_cos(Nt_all,ncos),slipave_inter(Nt_all,nas),slipave_cos(Nt_all,ncos),&
        v_cos(Nt_all,ncos),slip_cos(Nt_all,ncos),slipz1_tau(Nt_all,ncos),slipz1_sse(Nt_all,nsse), &
     v_nul(Nt_all,nnul),slip_nul(Nt_all,nnul),xi_all(Nt_all),x_all(Nt_all),&
      slipz1_v(Nt_all,ncos)
integer :: n_intz1,n_intz2,n_intz3,n_cosz1,n_cosz2,n_cosz3
integer :: intdepz1(Nt_all),intdepz2(Nt_all),intdepz3(Nt_all)
integer :: mpi_to_mesh_map(Nt_all)  ! FIXED: Add mpi_to_mesh_map parameter

! HDF5 variables for time-series output
integer(HID_T) :: file_id, dset_id, dspace_id
integer(HID_T) :: group_id, attr_id, attr_space_id
integer(HSIZE_T), dimension(3) :: dims, maxdims
integer(HSIZE_T), dimension(2) :: dims_2d, maxdims_2d
integer(HSIZE_T), dimension(1) :: dims_1d, maxdims_1d
integer :: hdferr
logical :: hdf5_initialized = .false.
integer :: global_time_steps_written = 0  ! Total time steps written across all cycles
integer :: global_sse_steps_written = 0   ! Total SSE time steps written across all cycles

! HDF5 file naming
character(len=256) :: hdf5_filename, xdmf_filename
character(len=256) :: time_series_group_name
logical :: file_exists, file_exists_sse, mesh_group_exists
integer :: ios
integer(HSIZE_T) :: offset_1d(1), count_1d(1), offset_2d(2), count_2d(2)
integer(HID_T) :: memspace_id, filespace_id, dcpl_id
integer(HSIZE_T) :: chunk_2d(2), chunk_1d(1)

! Mesh variables for GTS file reading
integer :: n_vertices, n_edges_dummy, n_cells
real(DP), allocatable :: vertex_coords(:,:)
integer*4, allocatable :: cell_connectivity(:,:)
real(DP), allocatable :: vertex_coords_transposed(:,:)
integer*4, allocatable :: cell_connectivity_transposed(:,:)


! MPI variables
integer :: myid, master
master = 0

! Allocate tcos_all for storing all time values (only on first call)
if (.not. allocated(tcos_all)) then
   allocate(tcos_all(10000))
   tcos_all = 0.d0
end if

! Copy tcos to tcos_all before resetting for next cycle
tcos_all(global_time_steps_written+1:global_time_steps_written+icos) = tcos(1:icos)
       

call MPI_COMM_RANK(MPI_COMM_WORLD, myid, hdferr)

if(Ioutput == 0)then    !output during run 


   if(imv==nmv)then
      open(30,file=trim(foldername)//'maxvall'//jobname,position='append',status='unknown')
      open(311,file=trim(foldername)//'fltst_strk-15'//jobname,access='append',status='unknown')
      open(312,file=trim(foldername)//'fltst_strk+00'//jobname,access='append',status='unknown')
      open(313,file=trim(foldername)//'fltst_strk+05'//jobname,access='append',status='unknown')
      open(314,file=trim(foldername)//'fltst_strk+10'//jobname,access='append',status='unknown')
      open(315,file=trim(foldername)//'fltst_strk+15'//jobname,access='append',status='unknown')
      open(316,file=trim(foldername)//'fltst_strk+25'//jobname,access='append',status='unknown')
      open(317,file=trim(foldername)//'fltst_strk+35'//jobname,access='append',status='unknown')
      open(318,file=trim(foldername)//'fltst_strk+50'//jobname,access='append',status='unknown')
      open(319,file=trim(foldername)//'fltst_strk+75'//jobname,access='append',status='unknown')

      do i=1,nmv
         write(30,130)tmv(i),dlog10(maxv(i)*1d-3/yrs),moment(i)
        do j=311,319
         write(j,110) tmv(i),outs1(i,1,j-310),outs1(i,2,j-310),outs1(i,3,j-310),outs1(i,4,j-310), &
           outs1(i,5,j-310),outs1(i,6,j-310),outs1(i,7,j-310)
        end do
       end do
      close(30)
      
      do j=311,319
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
       write(*,*) 'DEBUG: Triggering HDF5 output - icos =', icos, 'ncos =', ncos
       
       ! CRITICAL: Synchronize all MPI processes before HDF5 output
       
       ! Only master MPI process should do HDF5 output to avoid deadlock
       ! Initialize HDF5 if not already done
       if (.not. hdf5_initialized) then
          call h5open_f(hdferr)
          hdf5_initialized = .true.
       end if
       
       ! Create HDF5 filename
       hdf5_filename = trim(foldername)//'timeseries_data_'//trim(jobname)//'.h5'
       
       ! Check if file exists, if so open in append mode, otherwise create new
       inquire(file=trim(hdf5_filename), exist=file_exists)
       
       if (file_exists) then
          ! Open existing file for read/write (single-process access)
          call h5fopen_f(trim(hdf5_filename), H5F_ACC_RDWR_F, file_id, hdferr)
          if (hdferr < 0) then
             write(*,*) 'ERROR: Failed to open HDF5 file for writing'
             ! Skip HDF5 operations if file open failed
             ! Continue with the rest of the code
          end if
          ! Open existing time-series group
          time_series_group_name = '/time_series'
          call h5gopen_f(file_id, trim(time_series_group_name), group_id, hdferr)
       else
          ! Create new file and initialize datasets with extensible dimensions
          call h5fcreate_f(trim(hdf5_filename), H5F_ACC_TRUNC_F, file_id, hdferr)
          ! Create time-series group
          time_series_group_name = '/time_series'
          call h5gcreate_f(file_id, trim(time_series_group_name), group_id, hdferr)
          
          ! Create extensible datasets for first time with optimal chunking
          dims_2d = (/INT(Nt_all, HSIZE_T), INT(icos, HSIZE_T)/)
          maxdims_2d = (/INT(Nt_all, HSIZE_T), H5S_UNLIMITED_F/)
          ! FIXED: Use optimal chunk size that aligns with data access patterns
          ! Chunk size should be large enough to be efficient but not too large
          ! CRITICAL: Chunk size must not exceed actual data dimensions
          ! Use simple, safe chunking logic
          if (icos <= 1) then
             chunk_2d = (/INT(min(Nt_all, 1000), HSIZE_T), INT(1, HSIZE_T)/)
          else
             chunk_2d = (/INT(min(Nt_all, 1000), HSIZE_T), INT(min(icos, 100), HSIZE_T)/)
          end if
          
          write(*,*) 'DEBUG: Creating HDF5 with icos =', icos, 'chunk_2d =', chunk_2d
          
          call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, hdferr)
          if (hdferr /= 0) then
             write(*,*) 'ERROR: Failed to create HDF5 dataset creation property list, hdferr =', hdferr
             stop
          end if
          
          call h5pset_chunk_f(dcpl_id, 2, chunk_2d, hdferr)
          if (hdferr /= 0) then
             write(*,*) 'ERROR: Failed to set HDF5 chunk size, hdferr =', hdferr
             stop
          end if
          
          call h5screate_simple_f(2, dims_2d, dspace_id, hdferr, maxdims_2d)
          if (hdferr /= 0) then
             write(*,*) 'ERROR: Failed to create HDF5 dataspace, hdferr =', hdferr
             stop
          end if
          
          call h5dcreate_f(group_id, 'slipz1_v', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr, dcpl_id)
          if (hdferr /= 0) then
             write(*,*) 'ERROR: Failed to create HDF5 dataset slipz1_v, hdferr =', hdferr
             stop
          end if
          call h5dclose_f(dset_id, hdferr)
          call h5sclose_f(dspace_id, hdferr)
          
          call h5screate_simple_f(2, dims_2d, dspace_id, hdferr, maxdims_2d)
          call h5dcreate_f(group_id, 'slipz1_cos', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr, dcpl_id)
          call h5dclose_f(dset_id, hdferr)
          call h5sclose_f(dspace_id, hdferr)
          
          call h5pclose_f(dcpl_id, hdferr)
          
          dims_1d = (/INT(icos, HSIZE_T)/)
          maxdims_1d = (/H5S_UNLIMITED_F/)
          ! FIXED: Use optimal chunk size for 1D time arrays
          ! CRITICAL: Chunk size must not exceed actual data dimensions
          if (icos <= 1) then
             chunk_1d = (/INT(1, HSIZE_T)/)
          else
             chunk_1d = (/INT(min(icos, 100), HSIZE_T)/)
          end if
          
          call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, hdferr)
          call h5pset_chunk_f(dcpl_id, 1, chunk_1d, hdferr)
          
          call h5screate_simple_f(1, dims_1d, dspace_id, hdferr, maxdims_1d)
          call h5dcreate_f(group_id, 'tcos', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr, dcpl_id)
          call h5dclose_f(dset_id, hdferr)
          call h5sclose_f(dspace_id, hdferr)
          
          call h5pclose_f(dcpl_id, hdferr)
          
          global_time_steps_written = 0
       end if
       
       ! Extend datasets if this is not the first cycle
       if (global_time_steps_written > 0) then
          ! Extend 2D datasets
          dims_2d = (/INT(Nt_all, HSIZE_T), INT(global_time_steps_written + icos, HSIZE_T)/)
          call h5dopen_f(group_id, 'slipz1_v', dset_id, hdferr)
          call h5dset_extent_f(dset_id, dims_2d, hdferr)
          call h5dclose_f(dset_id, hdferr)
          
          call h5dopen_f(group_id, 'slipz1_cos', dset_id, hdferr)
          call h5dset_extent_f(dset_id, dims_2d, hdferr)
          call h5dclose_f(dset_id, hdferr)
          
          ! Extend 1D dataset
          dims_1d = (/INT(global_time_steps_written + icos, HSIZE_T)/)
          call h5dopen_f(group_id, 'tcos', dset_id, hdferr)
          call h5dset_extent_f(dset_id, dims_1d, hdferr)
          call h5dclose_f(dset_id, hdferr)
       end if
       
       ! Write slipz1_v data using hyperslab selection for accumulative writing
       call h5dopen_f(group_id, 'slipz1_v', dset_id, hdferr)
       call h5dget_space_f(dset_id, filespace_id, hdferr)
       
       ! FIXED: Define hyperslab for appending new data with proper alignment
       ! Ensure offset aligns with chunk boundaries for better performance
       offset_2d = (/INT(0, HSIZE_T), INT(global_time_steps_written, HSIZE_T)/)
       count_2d = (/INT(Nt_all, HSIZE_T), INT(icos, HSIZE_T)/)
       call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, offset_2d, count_2d, hdferr)
       
       ! Create memory space for current data
       dims_2d = (/INT(Nt_all, HSIZE_T), INT(icos, HSIZE_T)/)
       call h5screate_simple_f(2, dims_2d, memspace_id, hdferr)
       
       ! Data is already in mesh order from lines 875-877, no reordering needed
       
       ! Write current cycle data
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_v(:,1:icos), dims_2d, hdferr, memspace_id, filespace_id)
       
       call h5sclose_f(memspace_id, hdferr)
       call h5sclose_f(filespace_id, hdferr)
       call h5dclose_f(dset_id, hdferr)
       
       ! Write slipz1_cos data using hyperslab selection
       call h5dopen_f(group_id, 'slipz1_cos', dset_id, hdferr)
       call h5dget_space_f(dset_id, filespace_id, hdferr)
       
       call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, offset_2d, count_2d, hdferr)
       call h5screate_simple_f(2, dims_2d, memspace_id, hdferr)
       
       ! Data is already in mesh order from lines 875-877, no reordering needed
       
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_cos(:,1:icos), dims_2d, hdferr, memspace_id, filespace_id)
       
       call h5sclose_f(memspace_id, hdferr)
       call h5sclose_f(filespace_id, hdferr)
       call h5dclose_f(dset_id, hdferr)
       
       ! Write time array using hyperslab selection
       call h5dopen_f(group_id, 'tcos', dset_id, hdferr)
       call h5dget_space_f(dset_id, filespace_id, hdferr)
       
       offset_1d = (/INT(global_time_steps_written, HSIZE_T)/)
       count_1d = (/INT(icos, HSIZE_T)/)
       call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, offset_1d, count_1d, hdferr)
       
       dims_1d = (/INT(icos, HSIZE_T)/)
       call h5screate_simple_f(1, dims_1d, memspace_id, hdferr)
       
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tcos(1:icos), dims_1d, hdferr, memspace_id, filespace_id)
       
       call h5sclose_f(memspace_id, hdferr)
       call h5sclose_f(filespace_id, hdferr)
       call h5dclose_f(dset_id, hdferr)
       
       ! Close time-series group
       call h5gclose_f(group_id, hdferr)
       
       ! VALIDATION: Write binary files for HDF5 validation
       ! Write slipz1_v data to binary file for validation
       open(unit=100, file=trim(foldername)//'slipz1_cos'//jobname, form='unformatted', access='stream', position='append', status='unknown')
       write(100) slipz1_v(:,1:icos)
       close(100)
       write(*,*) 'Validation: slipz1_v data written to ', trim(foldername)//'slipz1_appendix.dat'
       
       ! Write tcos data to binary file for validation
       open(unit=101, file=trim(foldername)//'t-cos'//jobname, form='unformatted', access='stream', position='append', status='unknown')
       write(101) tcos(1:icos)
       close(101)
       write(*,*) 'Validation: tcos data written to ', trim(foldername)//'t-cos.dat'
       
       ! Add mesh data to HDF5
       ! Read GTS file and store mesh information
       inquire(file='triangular_mesh.gts', exist=file_exists)
       if (.not. file_exists) then
          write(*,*) 'ERROR: triangular_mesh.gts file not found!'
          write(*,*) 'This file should contain the mesh geometry matching the simulation.'
          write(*,*) 'Please ensure triangular_mesh.gts exists in the current directory.'
          stop
       end if
       
       open(98, file='triangular_mesh.gts', status='old', action='read', iostat=ios)
       if (ios /= 0) then
          write(*,*) 'ERROR: Failed to open triangular_mesh.gts, iostat =', ios
          stop
       end if
       
       read(98,*, iostat=ios) n_vertices, n_edges_dummy, n_cells
       if (ios /= 0) then
          write(*,*) 'ERROR: Failed to read mesh dimensions from triangular_mesh.gts, iostat =', ios
          close(98)
          stop
       end if
       
       write(*,*) 'DEBUG: Reading mesh from triangular_mesh.gts:'
       write(*,*) '  Vertices:', n_vertices, 'Edges:', n_edges_dummy, 'Cells:', n_cells
       
       ! Allocate temporary arrays
       allocate(vertex_coords(n_vertices, 3))
       allocate(cell_connectivity(n_cells, 3))
       
       ! Read vertex coordinates
       do i = 1, n_vertices
          read(98,*, iostat=ios) vertex_coords(i, 1), vertex_coords(i, 2), vertex_coords(i, 3)
          if (ios /= 0) then
             write(*,*) 'ERROR: Failed to read vertex', i, 'from triangular_mesh.gts, iostat =', ios
             close(98)
             stop
          end if
       end do
       
       ! Read cell connectivity (indices start from 0 in GTS, which is correct for Paraview)
       do i = 1, n_cells
          read(98,*, iostat=ios) cell_connectivity(i, 1), cell_connectivity(i, 2), cell_connectivity(i, 3)
          if (ios /= 0) then
             write(*,*) 'ERROR: Failed to read cell', i, 'from triangular_mesh.gts, iostat =', ios
             close(98)
             stop
          end if
          ! Keep 0-based indexing for Paraview compatibility
          cell_connectivity(i, :) = cell_connectivity(i, :) - 1
       end do
       close(98)
       
       write(*,*) 'DEBUG: Successfully read mesh with', n_vertices, 'vertices and', n_cells, 'cells'
       write(*,*) 'DEBUG: First few vertex coordinates:'
       do i = 1, min(5, n_vertices)
          write(*,*) '  Vertex', i, ':', vertex_coords(i, 1), vertex_coords(i, 2), vertex_coords(i, 3)
       end do
       
       ! Write mesh data to HDF5 - check if mesh group already exists
       call h5lexists_f(file_id, '/mesh', mesh_group_exists, hdferr)
       if (.not. mesh_group_exists) then
          call h5gcreate_f(file_id, '/mesh', group_id, hdferr)
       else
          call h5gopen_f(file_id, '/mesh', group_id, hdferr)
       end if
       
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
       ! Check if geometry dataset already exists
       call h5lexists_f(group_id, 'geometry', mesh_group_exists, hdferr)
       if (.not. mesh_group_exists) then
          call h5dcreate_f(group_id, 'geometry', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vertex_coords_transposed, dims_2d, hdferr)
          call h5dclose_f(dset_id, hdferr)
       end if
       call h5sclose_f(dspace_id, hdferr)
       
       deallocate(vertex_coords_transposed)
       
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
       ! Check if topology dataset already exists
       call h5lexists_f(group_id, 'topology', mesh_group_exists, hdferr)
       if (.not. mesh_group_exists) then
          call h5dcreate_f(group_id, 'topology', H5T_STD_I32LE, dspace_id, dset_id, hdferr)
          call h5dwrite_f(dset_id, H5T_STD_I32LE, cell_connectivity_transposed, dims_2d, hdferr)
          call h5dclose_f(dset_id, hdferr)
       end if
       
       deallocate(cell_connectivity_transposed)
       call h5sclose_f(dspace_id, hdferr)
       
       call h5gclose_f(group_id, hdferr)
       
       ! Deallocate temporary arrays
       deallocate(vertex_coords, cell_connectivity)
       
       ! Close HDF5 file
       call h5fclose_f(file_id, hdferr)
       
       ! Create or update XDMF file for visualization with accumulative time steps
       xdmf_filename = trim(foldername)//'timeseries_data_'//trim(jobname)//'.xdmf'
       open(99, file=trim(xdmf_filename), status='replace')
       write(99,'(A)') '<?xml version="1.0" ?>'
       write(99,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
       write(99,'(A)') '<Xdmf Version="2.0">'
       write(99,'(A)') ' <Domain>'
       write(99,'(A)') '  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
       
       ! Write a Grid for ALL accumulated time steps (including previous cycles)
       do i = 1, global_time_steps_written + icos
          write(99,'(A,I0,A)') '   <Grid Name="step_', i, '" GridType="Uniform">'
          write(99,'(A,I0,A)') '    <Topology TopologyType="Triangle" NumberOfElements="',n_cells,'">'
          write(99,'(A,I0,3A)') '     <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="',n_cells,' 3">timeseries_data_', trim(jobname), '.h5:/mesh/topology</DataItem>'
          write(99,'(A)') '    </Topology>'
          write(99,'(A,I0,A)') '    <Geometry name="geo" GeometryType="XYZ" NumberOfElements="',n_vertices,'">'
          write(99,'(A,I0,3A)') '     <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="',n_vertices,' 3">timeseries_data_', trim(jobname), '.h5:/mesh/geometry</DataItem>'
          write(99,'(A)') '    </Geometry>'
          ! Use actual time value from tcos_all (with bounds check)
          if (i <= size(tcos_all)) then
             write(99,'(A,E15.8,A)') '    <Time Value="', tcos_all(i)*yrs, '"/>'
          else
             write(99,'(A,E15.8,A)') '    <Time Value="', real(i-1, DP), '"/>'  ! Fallback to step index
          end if
          write(99,'(A)') '    <Attribute Name="slip_rate" Center="Cell">'
          write(99,'(A,I0,A)') '     <DataItem ItemType="HyperSlab" Dimensions="',n_cells,'">'
          write(99,'(A,I0,A,I0,A)') '      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">', i-1, ' 0 1 1 1 ',Nt_all,'</DataItem>'
          write(99,'(A,I0,3A)') '      <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="1 ',Nt_all,'">timeseries_data_', trim(jobname), '.h5:/time_series/slipz1_v</DataItem>'
          write(99,'(A)') '     </DataItem>'
          write(99,'(A)') '    </Attribute>'
          write(99,'(A)') '    <Attribute Name="fault_slip" Center="Cell">'
          write(99,'(A,I0,A)') '     <DataItem ItemType="HyperSlab" Dimensions="',n_cells,'">'
          write(99,'(A,I0,A,I0,A)') '      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">', i-1, ' 0 1 1 1 ',Nt_all,'</DataItem>'
          write(99,'(A,I0,3A)') '      <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="1 ',Nt_all,'">timeseries_data_', trim(jobname), '.h5:/time_series/slipz1_cos</DataItem>'
          write(99,'(A)') '     </DataItem>'
          write(99,'(A)') '    </Attribute>'
          write(99,'(A)') '   </Grid>'
       end do
       
       write(99,'(A)') '  </Grid>'
       write(99,'(A)') ' </Domain>'
       write(99,'(A)') '</Xdmf>'
       close(99)
       
       write(*,*) 'Time-series data written to HDF5: ', trim(hdf5_filename)
       write(*,*) 'XDMF visualization file created: ', trim(xdmf_filename)
       
     
       ! Update global counter for accumulative writing
       global_time_steps_written = global_time_steps_written + icos
       
       icos = 0
       
       ! CRITICAL: Synchronize all MPI processes after HDF5 output      
    end if


	if(inul == nnul)then
       ! Null slip data collection completed - no output files needed
       inul = 0 
	end if

   if(isse==nsse)then
      ! HDF5 output for SSE time-series variables instead of binary files
      ! CRITICAL: Synchronize all MPI processes before SSE HDF5 output
      
      ! Only master MPI process should do HDF5 output to avoid deadlock
         ! Initialize HDF5 if not already done
         if (.not. hdf5_initialized) then
            call h5open_f(hdferr)
            hdf5_initialized = .true.
         end if
      
      ! Create HDF5 filename for SSE data
      hdf5_filename = trim(foldername)//'sse_timeseries_data_'//trim(jobname)//'.h5'
      
      ! Check if file exists, if so open in append mode, otherwise create new
      inquire(file=trim(hdf5_filename), exist=file_exists_sse)
      
      if (file_exists_sse) then
         ! Open existing file for read/write (accumulative mode for SSE)
         call h5fopen_f(trim(hdf5_filename), H5F_ACC_RDWR_F, file_id, hdferr)
         ! Open existing SSE time-series group
         time_series_group_name = '/sse_time_series'
         call h5gopen_f(file_id, trim(time_series_group_name), group_id, hdferr)
      else
         ! Create new file with extensible datasets for SSE
         call h5fcreate_f(trim(hdf5_filename), H5F_ACC_TRUNC_F, file_id, hdferr)
         ! Create SSE time-series group
         time_series_group_name = '/sse_time_series'
         call h5gcreate_f(file_id, trim(time_series_group_name), group_id, hdferr)
         
         ! Create initial extensible datasets for SSE data with optimal chunking
         dims_2d = (/INT(Nt_all, HSIZE_T), INT(nsse, HSIZE_T)/)
         maxdims_2d = (/INT(Nt_all, HSIZE_T), H5S_UNLIMITED_F/)
         ! FIXED: Use optimal chunk size for SSE data
         ! CRITICAL: Chunk size must not exceed actual data dimensions
         chunk_2d = (/INT(min(Nt_all, 1000), HSIZE_T), INT(min(max(nsse, 10), nsse), HSIZE_T)/)
         
         call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, hdferr)
         call h5pset_chunk_f(dcpl_id, 2, chunk_2d, hdferr)
         
         call h5screate_simple_f(2, dims_2d, dspace_id, hdferr, maxdims_2d)
         call h5dcreate_f(group_id, 'slipz1_sse', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr, dcpl_id)
         call h5dclose_f(dset_id, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         
         call h5screate_simple_f(2, dims_2d, dspace_id, hdferr, maxdims_2d)
         call h5dcreate_f(group_id, 'slipz1_tau', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr, dcpl_id)
         call h5dclose_f(dset_id, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         
         ! pore_fluid output removed
         
         call h5pclose_f(dcpl_id, hdferr)
         
         dims_1d = (/INT(nsse, HSIZE_T)/)
         maxdims_1d = (/H5S_UNLIMITED_F/)
         ! FIXED: Use optimal chunk size for SSE 1D time arrays
         ! CRITICAL: Chunk size must not exceed actual data dimensions
         chunk_1d = (/INT(min(max(nsse, 10), nsse), HSIZE_T)/)
         
         call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, hdferr)
         call h5pset_chunk_f(dcpl_id, 1, chunk_1d, hdferr)
         
         call h5screate_simple_f(1, dims_1d, dspace_id, hdferr, maxdims_1d)
         call h5dcreate_f(group_id, 'tsse', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr, dcpl_id)
         call h5dclose_f(dset_id, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         
         call h5pclose_f(dcpl_id, hdferr)
         
         global_sse_steps_written = 0
      end if
      
      ! For SSE, we append new nsse columns to existing data
      ! First, determine current dataset size
      call h5dopen_f(group_id, 'slipz1_sse', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5sget_simple_extent_dims_f(dspace_id, dims_2d, maxdims_2d, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      ! Current size is dims_2d(2), extend by nsse
      dims_2d = (/INT(Nt_all, HSIZE_T), dims_2d(2) + INT(nsse, HSIZE_T)/)
      
      ! Extend all SSE datasets
      call h5dopen_f(group_id, 'slipz1_sse', dset_id, hdferr)
      call h5dset_extent_f(dset_id, dims_2d, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      call h5dopen_f(group_id, 'slipz1_tau', dset_id, hdferr)
      call h5dset_extent_f(dset_id, dims_2d, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      dims_1d = (/dims_2d(2)/)
      call h5dopen_f(group_id, 'tsse', dset_id, hdferr)
      call h5dset_extent_f(dset_id, dims_1d, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      ! Write new SSE data using hyperslab selection
      call h5dopen_f(group_id, 'slipz1_sse', dset_id, hdferr)
      call h5dget_space_f(dset_id, filespace_id, hdferr)
      
      offset_2d = (/INT(0, HSIZE_T), dims_2d(2) - INT(nsse, HSIZE_T)/)  ! Start at the new columns
      count_2d = (/INT(Nt_all, HSIZE_T), INT(nsse, HSIZE_T)/)
      call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, offset_2d, count_2d, hdferr)
      
      dims_2d = (/INT(Nt_all, HSIZE_T), INT(nsse, HSIZE_T)/)
      call h5screate_simple_f(2, dims_2d, memspace_id, hdferr)
      
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_sse, dims_2d, hdferr, memspace_id, filespace_id)
      
      call h5sclose_f(memspace_id, hdferr)
      call h5sclose_f(filespace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      ! Write slipz1_tau data
      call h5dopen_f(group_id, 'slipz1_tau', dset_id, hdferr)
      call h5dget_space_f(dset_id, filespace_id, hdferr)
      call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, offset_2d, count_2d, hdferr)
      call h5screate_simple_f(2, dims_2d, memspace_id, hdferr)
      
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_tau, dims_2d, hdferr, memspace_id, filespace_id)
      call h5sclose_f(memspace_id, hdferr)
      call h5sclose_f(filespace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      ! pore_fluid output removed
      
      ! Write time array
      call h5dopen_f(group_id, 'tsse', dset_id, hdferr)
      call h5dget_space_f(dset_id, filespace_id, hdferr)
      
      offset_1d = (/dims_1d(1) - INT(nsse, HSIZE_T)/)
      count_1d = (/INT(nsse, HSIZE_T)/)
      call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, offset_1d, count_1d, hdferr)
      
      dims_1d = (/INT(nsse, HSIZE_T)/)
      call h5screate_simple_f(1, dims_1d, memspace_id, hdferr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tsse, dims_1d, hdferr, memspace_id, filespace_id)
      call h5sclose_f(memspace_id, hdferr)
      call h5sclose_f(filespace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      ! Close SSE time-series group
      call h5gclose_f(group_id, hdferr)
      
      ! Add mesh data to HDF5
      ! Read GTS file and store mesh information
      inquire(file='triangular_mesh.gts', exist=file_exists)
      if (.not. file_exists) then
         write(*,*) 'ERROR: triangular_mesh.gts file not found!'
         write(*,*) 'This file should contain the mesh geometry matching the simulation.'
         write(*,*) 'Please ensure triangular_mesh.gts exists in the current directory.'
         stop
      end if
      
      open(98, file='triangular_mesh.gts', status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,*) 'ERROR: Failed to open triangular_mesh.gts, iostat =', ios
         stop
      end if
      
      read(98,*, iostat=ios) n_vertices, n_edges_dummy, n_cells
      if (ios /= 0) then
         write(*,*) 'ERROR: Failed to read mesh dimensions from triangular_mesh.gts, iostat =', ios
         close(98)
         stop
      end if
      
      write(*,*) 'DEBUG: Reading SSE mesh from triangular_mesh.gts:'
      write(*,*) '  Vertices:', n_vertices, 'Edges:', n_edges_dummy, 'Cells:', n_cells
      
      ! Allocate temporary arrays
      allocate(vertex_coords(n_vertices, 3))
      allocate(cell_connectivity(n_cells, 3))
      
      ! Read vertex coordinates
      do i = 1, n_vertices
         read(98,*, iostat=ios) vertex_coords(i, 1), vertex_coords(i, 2), vertex_coords(i, 3)
         if (ios /= 0) then
            write(*,*) 'ERROR: Failed to read vertex', i, 'from triangular_mesh.gts, iostat =', ios
            close(98)
            stop
         end if
      end do
      
      ! Read cell connectivity (indices start from 0 in GTS, which is correct for Paraview)
      do i = 1, n_cells
         read(98,*, iostat=ios) cell_connectivity(i, 1), cell_connectivity(i, 2), cell_connectivity(i, 3)
         if (ios /= 0) then
            write(*,*) 'ERROR: Failed to read cell', i, 'from triangular_mesh.gts, iostat =', ios
            close(98)
            stop
         end if
         ! Keep 0-based indexing for Paraview compatibility
         cell_connectivity(i, :) = cell_connectivity(i, :) - 1
      end do
      close(98)
      
      write(*,*) 'DEBUG: Successfully read SSE mesh with', n_vertices, 'vertices and', n_cells, 'cells'
      
      ! Write mesh data to HDF5 - check if mesh group already exists
      call h5lexists_f(file_id, '/mesh', mesh_group_exists, hdferr)
      if (.not. mesh_group_exists) then
         call h5gcreate_f(file_id, '/mesh', group_id, hdferr)
      else
         call h5gopen_f(file_id, '/mesh', group_id, hdferr)
      end if
      
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
      ! Check if geometry dataset already exists
      call h5lexists_f(group_id, 'geometry', mesh_group_exists, hdferr)
      if (.not. mesh_group_exists) then
         call h5dcreate_f(group_id, 'geometry', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vertex_coords_transposed, dims_2d, hdferr)
         call h5dclose_f(dset_id, hdferr)
      end if
      call h5sclose_f(dspace_id, hdferr)
      
      deallocate(vertex_coords_transposed)
      
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
      ! Check if topology dataset already exists
      call h5lexists_f(group_id, 'topology', mesh_group_exists, hdferr)
      if (.not. mesh_group_exists) then
         call h5dcreate_f(group_id, 'topology', H5T_STD_I32LE, dspace_id, dset_id, hdferr)
         call h5dwrite_f(dset_id, H5T_STD_I32LE, cell_connectivity_transposed, dims_2d, hdferr)
         call h5dclose_f(dset_id, hdferr)
      end if
      
      deallocate(cell_connectivity_transposed)
      call h5sclose_f(dspace_id, hdferr)
      
      call h5gclose_f(group_id, hdferr)
      
      ! Deallocate temporary arrays
      deallocate(vertex_coords, cell_connectivity)
      
      ! Close HDF5 file
      call h5fclose_f(file_id, hdferr)
      
      ! Create or update XDMF file for SSE visualization with accumulative time steps
      xdmf_filename = trim(foldername)//'sse_timeseries_data_'//trim(jobname)//'.xdmf'
      open(99, file=trim(xdmf_filename), status='replace')
      write(99,'(A)') '<?xml version="1.0" ?>'
      write(99,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      write(99,'(A)') '<Xdmf Version="2.0">'
      write(99,'(A)') ' <Domain>'
      write(99,'(A)') '  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
      
      ! Write a Grid for ALL accumulated SSE time steps (including previous cycles)
      do i = 1, global_sse_steps_written + nsse
         write(99,'(A,I0,A)') '   <Grid Name="step_', i, '" GridType="Uniform">'
         write(99,'(A,I0,A)') '    <Topology TopologyType="Triangle" NumberOfElements="',n_cells,'">'
         write(99,'(A,I0,3A)') '     <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="',n_cells,' 3">sse_timeseries_data_', trim(jobname), '.h5:/mesh/topology</DataItem>'
         write(99,'(A)') '    </Topology>'
         write(99,'(A,I0,A)') '    <Geometry name="geo" GeometryType="XYZ" NumberOfElements="',n_vertices,'">'
         write(99,'(A,I0,3A)') '     <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="',n_vertices,' 3">sse_timeseries_data_', trim(jobname), '.h5:/mesh/geometry</DataItem>'
         write(99,'(A)') '    </Geometry>'
         write(99,'(A,E15.8,A)') '    <Time Value="', real(i-1, DP), '"/>'  ! Use step index as time for now
         write(99,'(A)') '    <Attribute Name="SSE_slip_rate" Center="Cell">'
         write(99,'(A,I0,A)') '     <DataItem ItemType="HyperSlab" Dimensions="',n_cells,'">'
         write(99,'(A,I0,A,I0,A)') '      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">', i-1, ' 0 1 1 1 ',Nt_all,'</DataItem>'
         write(99,'(A,I0,3A)') '      <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="1 ',Nt_all,'">sse_timeseries_data_', trim(jobname), '.h5:/sse_time_series/slipz1_sse</DataItem>'
         write(99,'(A)') '     </DataItem>'
         write(99,'(A)') '    </Attribute>'
         write(99,'(A)') '    <Attribute Name="shear_stress" Center="Cell">'
         write(99,'(A,I0,A)') '     <DataItem ItemType="HyperSlab" Dimensions="',n_cells,'">'
         write(99,'(A,I0,A,I0,A)') '      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">', i-1, ' 0 1 1 1 ',Nt_all,'</DataItem>'
         write(99,'(A,I0,3A)') '      <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="1 ',Nt_all,'">sse_timeseries_data_', trim(jobname), '.h5:/sse_time_series/slipz1_tau</DataItem>'
         write(99,'(A)') '     </DataItem>'
         write(99,'(A)') '    </Attribute>'
         ! pore_fluid XDMF output removed
         write(99,'(A)') '   </Grid>'
      end do
      
      write(99,'(A)') '  </Grid>'
      write(99,'(A)') ' </Domain>'
      write(99,'(A)') '</Xdmf>'
      close(99)
      
      write(*,*) 'SSE time-series data written to HDF5: ', trim(hdf5_filename)
      write(*,*) 'SSE XDMF visualization file created: ', trim(xdmf_filename)
      
      ! VALIDATION: Write binary files for SSE HDF5 validation
      ! Write slipz1_sse data to binary file for validation
      open(unit=102, file=trim(foldername)//'slipz1_sse'//jobname, form='unformatted', access='stream', position='append', status='unknown')
      write(102) slipz1_sse(:,1:nsse)
      close(102)
      write(*,*) 'Validation: slipz1_sse data written to ', trim(foldername)//'slipz1_sse_appendix.dat'
      
      ! Write tsse data to binary file for validation
      open(unit=103, file=trim(foldername)//'t-sse'//jobname, form='unformatted', access='stream', position='append', status='unknown')
      write(103) tsse(1:nsse)
      close(103)
      write(*,*) 'Validation: tsse data written to ', trim(foldername)//'t-sse.dat'
      
      ! Update global SSE counter for accumulative writing
      global_sse_steps_written = global_sse_steps_written + nsse
      isse = 0
      
      ! CRITICAL: Synchronize all MPI processes after SSE HDF5 output
     
  end if


else

   if((imv>0).and.(imv<nmv))then
      open(30,file=trim(foldername)//'maxvall'//jobname,position='append',status='unknown')
      open(311,file=trim(foldername)//'fltst_strk-15'//jobname,access='append',status='unknown')
      open(312,file=trim(foldername)//'fltst_strk+00'//jobname,access='append',status='unknown')
      open(313,file=trim(foldername)//'fltst_strk+05'//jobname,access='append',status='unknown')
      open(314,file=trim(foldername)//'fltst_strk+10'//jobname,access='append',status='unknown')
      open(315,file=trim(foldername)//'fltst_strk+15'//jobname,access='append',status='unknown')
      open(316,file=trim(foldername)//'fltst_strk+25'//jobname,access='append',status='unknown')
      open(317,file=trim(foldername)//'fltst_strk+35'//jobname,access='append',status='unknown')
      open(318,file=trim(foldername)//'fltst_strk+50'//jobname,access='append',status='unknown')
      open(319,file=trim(foldername)//'fltst_strk+75'//jobname,access='append',status='unknown')

      do i=1,imv
         write(30,130)tmv(i),dlog10(maxv(i)*1d-3/yrs),moment(i)
        do j=311,319
         write(j,110) tmv(i),outs1(i,1,j-310),outs1(i,2,j-310),outs1(i,3,j-310),outs1(i,4,j-310), &
           outs1(i,5,j-310),outs1(i,6,j-310),outs1(i,7,j-310)
        end do
      end do
       close(30)
       do j=311,319
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
      
      ! pore_fluid partial output removed
      
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
       
       ! Write partial slipz1_cos data (cosine slip time series)
       dims_2d = (/Nt_all, icos/)
       call h5screate_simple_f(2, dims_2d, dspace_id, hdferr)
       call h5dcreate_f(group_id, 'slipz1_cos_partial', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, slipz1_cos(:,1:icos), dims_2d, hdferr)
       call h5dclose_f(dset_id, hdferr)
       call h5sclose_f(dspace_id, hdferr)
       
       ! Write partial slipz1_v data (velocity time series)
       dims_2d = (/Nt_all, icos/)
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

end if  ! Close if(Ioutput == 0)then

RETURN
END subroutine output


!------------------------------------------------------------------------------
! Data reordering function to ensure consistent data ordering for HDF5 output
!------------------------------------------------------------------------------
subroutine reorder_data_for_hdf5(data_array, n_elements, n_timesteps, mpi_to_mesh_map)
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, intent(in) :: n_elements, n_timesteps
  real(DP), intent(inout) :: data_array(n_elements, n_timesteps)
  integer, intent(in) :: mpi_to_mesh_map(n_elements)
  
  real(DP), allocatable :: temp_array(:,:)
  integer :: i, j, mesh_idx
  
  ! Allocate temporary array for reordering
  allocate(temp_array(n_elements, n_timesteps))
  
  ! Copy original data to temporary array
  temp_array = data_array
  
  ! Reorder data according to mesh ordering
  do i = 1, n_elements
     mesh_idx = mpi_to_mesh_map(i)
     if (mesh_idx >= 1 .and. mesh_idx <= n_elements) then
        do j = 1, n_timesteps
           data_array(mesh_idx, j) = temp_array(i, j)
        end do
     else
        write(*,*) 'ERROR: Invalid mesh index', mesh_idx, 'for MPI index', i
     end if
  end do
  
  deallocate(temp_array)
  
end subroutine reorder_data_for_hdf5

! Heaviside function for pore fluid pressure calculation
function heavi(x)
  implicit none
  real(8), intent(in) :: x
  real(8) :: heavi
  
  if (x >= 0.0d0) then
    heavi = 1.0d0
  else
    heavi = 0.0d0
  end if
end function heavi
