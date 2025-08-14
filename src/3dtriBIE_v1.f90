!BSD 3-Clause License
!
!Copyright (c) 2019, SeisSol Group
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
!
!* Redistributions of source code must retain the above copyright notice, this
!  list of conditions and the following disclaimer.
!
!* Redistributions in binary form must reproduce the above copyright notice,
!  this list of conditions and the following disclaimer in the documentation
!  and/or other materials provided with the distribution.
!
!* Neither the name of the copyright holder nor the names of its
!  contributors may be used to endorse or promote products derived from
!  this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!------------------------------------------------
! code: TriBIE
! Purpose: 3D Subduction fault earthquake and aseismic sequence model
! Language: Fortran + OpenMPI 
! Authors:
! D. Li (liduoduo07@gmail.com)
! Last modified Jan. 2019

! Yajing Liu (yajing.liu@mcgill.ca)
! Last modified Jul. 4, 2012
!
! Notes:
! using varying plate convergence rate. pay attention to theta !!!!
!  yt(2*i) should not be negative!!!!
!  The fault could be set up with two buffer zones (e.g. 20 km wide each)
! cleaned up for higher efficiency
!------------------------------------------------

program main
  USE mpi
  USE phy3d_module_non
  implicit none
  integer, parameter :: DP=kind(1.d0)

  logical cyclecont
  integer :: Nt_all,Nt,ndt,ndtnext,kk,ii,n,l,ndt_v,ndt_inter,Itout,&  
       ndtmax,i,j,k,nm,nrec,ihrestart,Ifileout, &
       Iperb,record,Isnapshot,record_cor,s1,s2,s3,s4,s5,s6,s7,s8

  real (DP) :: accuracy,areatot, epsv,dt_try, dt,dtmin,dt_did,dt_next,dip, &
       hnucl, sigmadiff,sefffix,Lffix,xdepth,xlength,&
       t,tprint_inter, tint_out,tout,tmin_out,tint_cos,&   
       tslip_ave,tslip_aveint, tmax, &
       tslipsse,tslipcos,tstart1,tend1,tstart2,tend2,tstart3,tend3,tssestart,tsseend, &
       factor1,factor2,factor3,factor4,vtmp,fr_dummy,&
       xilock1,xilock2,xi_dummy,x_dummy,z_dummy,stiff_dummy,help

  real (DP) ::  tmbegin,tmrun

  real (DP), DIMENSION(:), ALLOCATABLE :: x,z,xi,yt,dydt,yt_scale,tau,tau0, &
       slip,slipinc

  !Arrays only defined at master cpu
  real (DP), DIMENSION(:), ALLOCATABLE :: x_all,xi_all,z_all,&
       yt_all,dydt_all,yt_scale_all,tau_all, &
       slip_all,slipinc_all,cca_all,ccb_all,xLf_all,seff_all

  !output related parameters
  integer :: imv,ias,icos,isse,Ioutput,inul,i_nul
  real (DP) :: vcos,vsse1,vsse2
  real (DP), DIMENSION(:), ALLOCATABLE :: maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,&
       maxnum,tmv,tas,tcos,tsse, maxv_s6,maxv_s7,maxv_s8

  real (DP), DIMENSION(:,:), ALLOCATABLE :: slipz1_inter, slipz1_v,slipz1_cos,slipz1_tau,slipz1_sse

  character(len=40) :: cTemp,filename,ct

  !MPI RELATED DEFINITIONS
  integer :: ierr,size,myid,master

  call MPI_Init(ierr)
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, size, ierr )

  master = 0 

  !read in intiliazation parameters

  open(12,file='./parameter.txt',form='formatted',status='old')

  read(12,'(a)')jobname
  read(12,'(a)')foldername
  read(12,'(a)')stiffname
  read(12,'(a)')restartname
  read(12,'(a)')profile
  read(12,*)Nab,Nt_all,Nt,Lratio,nprocs,hnucl
  read(12,*)Idin,Idout,Iprofile,Iperb,Isnapshot 
  read(12,*)factor1,factor2,factor3,factor4
  read(12,*)Vpl
  read(12,*)xilock1,xilock2
  read(12,*)sigmadiff,sefffix,Lffix
  read(12,*)tmax
  read(12,*)tslip_ave,tslip_aveint
  read(12,*)tint_out,tmin_out,tint_cos
  read(12,*)vcos,vsse1,vsse2
  read(12,*) tssestart,tsseend
  read(12,*) nmv,nas,ncos,nsse
  read(12,*) s1,s2,s3,s4,s5,s6,s7,s8

    110 format(A)
    close(12)


  if(mod(Nt_all,nprocs)/=0)then
     write(*,*)'Nd_all must be integer*nprocs. Change nprocs!'
     STOP
  else
     write(*,*)'Each cpu calculates',Nt_all/nprocs,'cells'
  end if
  
  ! hnucl = hstarfactor*hd            !h* (how about along strike?)
  if(myid == master)then
     ALLOCATE(x_all(Nt_all),xi_all(Nt_all),&
          cca_all(Nt_all),ccb_all(Nt_all),seff_all(Nt_all),xLf_all(Nt_all),&
          tau_all(Nt_all),slip_all(Nt_all),slipinc_all(Nt_all),yt_all(2*Nt_all), &
          dydt_all(2*Nt_all),yt_scale_all(2*Nt_all))

     ALLOCATE (maxv_s1(nmv),maxv_s2(nmv),maxv_s3(nmv),maxv_s4(nmv),maxv_s5(nmv),&
          maxv_s6(nmv),maxv_s7(nmv),maxv_s8(nmv),&
          maxv(nmv),maxnum(nmv), &
          tmv(nmv),tas(nas),tcos(ncos),tsse(nsse))

!!! modify output number
     ALLOCATE (slipz1_inter(Nt_all,nas),slipz1_cos(Nt_all,ncos), &
          slipz1_tau(Nt_all,nsse),slipz1_sse(Nt_all,nsse) )
     ALLOCATE( slipz1_v(Nt_all,ncos) )
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE (x(Nt),z(Nt),z_all(Nt_all),xi(Nt),cca(Nt),ccb(Nt),seff(Nt),&
       xLf(Nt),tau(Nt),tau0(Nt),slip(Nt),slipinc(Nt), &
       yt(2*Nt),dydt(2*Nt),yt_scale(2*Nt))

  ALLOCATE (stiff(Nt,Nt_all))   !!! stiffness of Stuart green calculation

  !Read in stiffness matrix, in nprocs segments

  write(cTemp,*) myid
  write(*,*) cTemp

  open(5, file=trim(stiffname)//'trigreen_'//trim(adjustl(cTemp))//'.bin',form='unformatted',access='stream')
  open(55,file=trim(stiffname)//'position.bin',form='unformatted',access='stream')

  record_cor=Nt*(myid)

  !-------------------------------------------------------------------------------------------
  ! read stiffness from Stuart green calculation.
  !-----------------------------------------------------------------------------------------
  if(myid==master)then
     do k=1,Nt_all
        read(55) xi_all(k),x_all(k),z_all(k) !xi is along the dip while x is along the strike
        xi_all(k)=xi_all(k)/1000
        x_all(k)=x_all(k)/1000
        z_all(k)=z_all(k)/1000
     end do
  end if

  do i=1,Nt !! observe
     do j=1,Nt_all !! source
        read(5) stiff(i,j)
     end do
  end do
  close(5)
  close(55)
  !!-----------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  if(myid==master)then
     CALL resdep(Nt_all,dip,hnucl,sigmadiff,sefffix,Lffix, &
           Iperb,factor1,factor2,factor3,factor4,&
          xilock1,xilock2,cca_all,ccb_all,xLf_all,seff_all,x_all,z_all)
  end if

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Scatter(cca_all,Nt,MPI_DOUBLE_PRECISION,cca,Nt,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  call MPI_Scatter(ccb_all,Nt,MPI_DOUBLE_PRECISION,ccb,Nt,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  call MPI_Scatter(xLf_all,Nt,MPI_DOUBLE_PRECISION,xLf,Nt,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  call MPI_Scatter(seff_all,Nt,MPI_DOUBLE_PRECISION,seff,Nt,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  call MPI_Scatter(x_all,Nt,MPI_DOUBLE_PRECISION,x,Nt,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(z_all,Nt_all,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)

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
!!! accuracy in Adaptive Rudge-Kutta
  accuracy = 1.d-3
  epsv = accuracy
  dtmin = 1.d-6
  dt_try=1.d-6


  if(myid==master)then

     open(1,file=trim(foldername)//'phypara'//jobname, status='unknown')
     write(1,*)'procs num = ', nprocs

  end if

  Ifileout = 60   !file index, after 47

  !----Initial values of velocity, state variable, shear stress and slip--
  !--SET INITIAL VPL FOR THE LOCKED PART TO BE 0 
  ! ! set plate convergence
  if(IDin.eq.0)then 
     do j=1,Nt
        yt(2*j-1)=0.8*Vpl 
        yt(2*j)=xLf(j)/(1.2*yt(2*j-1))	
        slip(j)=0.d0
     end do
  end if

  !------------------------------------------------------------------
  if(IDin.eq.1) then               !if this is a restart job
     if(myid==master)then
        call restart(0,'out',4,Nt_all,t,dt,dt_try,ndt,nrec,yt_all,slip_all)
        write(1,*)'This is a restart job. Start time ',t,' yr'
     end if
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_Bcast(t,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(dt,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(dt_try,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ndt,1,MPI_integer,master,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(nrec,1,MPI_integer,master,MPI_COMM_WORLD,ierr)
     call MPI_Scatter(yt_all,2*Nt,MPI_DOUBLE_PRECISION,yt,2*Nt,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
     call MPI_Scatter(slip_all,Nt,MPI_DOUBLE_PRECISION,slip,Nt,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)

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

     call derivs(myid,dydt,2*Nt,Nt_all,Nt,t,yt,z_all,x) 

     do j=1,2*Nt
        yt_scale(j)=dabs(yt(j))+dabs(dt_try*dydt(j))
     end do
     CALL rkqs(myid,yt,dydt,2*Nt,Nt_all,Nt,t,dt_try,accuracy,yt_scale, &
          dt_did,dt_next,z_all,x)

     dt = dt_did
     dt_try = dt_next

     do i=1,Nt
        if(yt(2*i-1).lt.epsv)then
           help=(yt(2*i-1)/(2*V0))*dexp((f0+ccb(i)*dlog(V0*yt(2*i)))/cca(i))
           tau(i)=seff(i)*cca(i)*dlog(help+dsqrt(1+help**2))
        else
           tau(i)=seff(i)*(f0+cca(i)*log(yt(2*i-1)/v0)+ccb(i)*log(V0*yt(2*i)/xLf(i)))
        end if
        slipinc(i) = yt(2*i-1)*dt
        slip(i) = slip(i) + slipinc(i)
     end do

     ndt = ndt + 1

     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_Gather(yt,2*Nt,MPI_DOUBLE_PRECISION,yt_all,2*Nt,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
     call MPI_Gather(slipinc,Nt,MPI_DOUBLE_PRECISION,slipinc_all,Nt,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
     call MPI_Gather(slip,Nt,MPI_DOUBLE_PRECISION,slip_all,Nt,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
     call MPI_Gather(tau,Nt,MPI_DOUBLE_PRECISION,tau_all,Nt,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)

     !-------------------
     !      Output:     (a single thread will do the writing while others 
     !      proceed 
     !-------------------

     if(myid==master)then
        imv=imv+1
        tmv(imv)=t
        maxv(imv) = 0.d0
        maxv_s1(imv)=0.d0
        maxv_s2(imv)=0.d0
        maxv_s3(imv)=0.d0
        maxv_s4(imv)=0.d0
        maxv_s5(imv)=0.d0
        maxv_s6(imv)=0.d0
        maxv_s7(imv)=0.d0
        maxv_s8(imv)=0.d0
        do i = 1,Nt_all      !!!!!!!!!!!!!!!!! find max velocity
           if(yt_all(2*i-1).ge.maxv(imv))then
              maxv(imv) = yt_all(2*i-1)
              maxnum(imv)=i
           end if
        end do
     write(*,*) ndt,t,maxv(imv)

        !! write fault observations
              maxv_s1(imv) = yt_all(2*s1-1)
              maxv_s2(imv) = yt_all(2*s2-1)
              maxv_s3(imv) = yt_all(2*s3-1)
              maxv_s4(imv) = yt_all(2*s4-1)
              maxv_s5(imv) = yt_all(2*s5-1)
              maxv_s6(imv) = yt_all(2*s6-1)
              maxv_s7(imv) = yt_all(2*s7-1)
              maxv_s8(imv) = yt_all(2*s8-1)

!!!!! checkout point  - not necessary now !!!!!!!!!!!!!!!

        !-----Interseismic slip every ? years----
        if (t.ge.tslip_ave)then
           ias = ias + 1 
           tas(ias)=t

           do i=1,Nt_all
              slipz1_inter(i,ias) = slip_all(i)*1.d-3           
           end do

           tslip_ave = tslip_ave + tslip_aveint
        end if

        !------velocity and slip of eq nucleation process  (deleted )-- 

        !----find SSE slip ------
        if(t.ge.tssestart.and.t.le.tsseend)then

           if((maxv(imv)/yrs).le.vcos.and.maxv(imv).ge.vsse1)then
              tslipsse = tslipsse + dt
              if(tslipsse.ge.0.002739726)then   !! output step = 1 day
            write(*,*)'SSE',t,maxv(imv)
                 isse = isse + 1
                 do i=1,Nt_all
                    slipz1_tau(i,isse)= tau_all(i)
                    slipz1_sse(i,isse)= dlog10(yt_all(2*i-1)*1.d-3/yrs)
                 end do
                 tsse(isse) = t
                 tslipsse = 0.d0
              end if
           end if
        end if

        !!!!!!!!!!!!!!!!------coseismic Slip   -------------------

        if((maxv(imv)/yrs).ge.vcos)then
           tslipcos = tslipcos+dt
           if(tslipcos.ge.tint_cos)then    !!!!tint_cos: every 5s output 
		write(*,*) 'coseis',t,maxv(imv)
              icos = icos +1
              tcos(icos) = t 

!!! find the total number of output elements and cell number of specific depth.

              do i=1,Nt_all
                 slipz1_cos(i,icos) = slip_all(i)*1.d-3
                 slipz1_v(i,icos) = dlog10(yt_all(2*i-1)*1.d-3/yrs) 
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
             tmv,tas,tcos,tsse,maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,maxv_s6,&
             maxv_s7,maxv_s8,maxnum, &
             slipz1_inter,slipz1_tau,slipz1_sse,slipz1_cos, xi_all,x_all,slipz1_v)         

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
          tmv,tas,tcos,tsse,maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,maxv_s6,maxv_s7,maxv_s8, &
          maxnum, &
          slipz1_inter,slipz1_tau,slipz1_sse, &
          slipz1_cos,&
          xi_all,x_all,slipz1_v) 

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
     write(10,'(T9,A,T45,F10.5)')'total running time(min)', tmrun/60.0
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
     DEALLOCATE (x_all,xi_all,yt_all,dydt_all,yt_scale_all,tau_all, &
          slip_all,slipinc_all,cca_all,ccb_all,xLf_all,seff_all, &
          maxnum,maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,maxv_s6,maxv_s7,maxv_s8,&
          tmv,tcos,tas,tsse)

     DEALLOCATE (slipz1_inter,slipz1_tau,slipz1_sse, &
          slipz1_cos,slipz1_v)
!     DEALLOCATE (intdepz1,intdepz2,intdepz3,ssetime,slipz1_v)
  end if

  DEALLOCATE (stiff)
  DEALLOCATE (x,z_all,xi,yt,dydt,yt_scale,tau,tau0,slip,slipinc)
  DEALLOCATE (cca,ccb,xLf,seff)
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
  real (DP) :: errmax,errmax1,h,htemp,xnew
  real (DP), dimension(:), allocatable :: yerr,ytemp
  real (DP), parameter :: SAFETY=0.9, PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4

  !MPI RELATED DEFINITIONS
  integer :: ierr,myid,master
  master = 0 

  nmax=2*Nt
  h=htry
  allocate (yerr(nmax),ytemp(nmax))
  !        write(*,*)'in rkqs, bf rkck, myid=',myid,y(n-1)
    1 call rkck(myid,dydx,h,n,Nt_all,Nt,y,yerr,ytemp,x,derivs,z_all,p)
  !        write(*,*)'in rkqs, af rkck, myid=',myid, ytemp(n-1)
  errmax=0.
  do  i=1,nmax
     j = int(ceiling(real(i)/2)) ! position within central part
     !   if(p(j).gt.-700.and.p(j).lt.180)then
     errmax = max(errmax,dabs(yerr(i)/yscal(i)))
     !   end if
  end do
  errmax=errmax/eps
!  write(*,*) errmax
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  call MPI_Allreduce(errmax, errmax1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

  if(errmax1.gt.1.)then
     htemp = SAFETY*h*(errmax1**PSHRNK)
     h = dsign(max(dabs(htemp),0.1*dabs(h)),h)
     xnew = x+h
     if(xnew.eq.x)pause 'stepsize underflow in rkqs'
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


       nmax = 2*Nt
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
       USE phy3d_module_non, only: stiff,cca,ccb,seff,xLf,eta,f0,Vpl,V0,Lratio,nprocs,&
            tm1,tm2,tmday,tmelse,tmmidn,tmmult
       implicit none
       integer, parameter :: DP = kind(1.0d0)
       integer :: nv,n,i,j,k,kk,l,ii,Nt,Nt_all
       real (DP) :: t,yt(nv),dydt(nv)   
       real (DP) :: deriv3,deriv2,deriv1,small
       real (DP) :: psi,help1,help2,help
       real (DP) :: SECNDS
       real (DP) :: z_all(Nt_all),x(Nt),zz(Nt),zzfric(Nt),zz_all(Nt_all)
       intrinsic imag,real

       !MPI RELATED DEFINITIONS
       integer :: ierr,myid,master
       master = 0 


       small=1.d-3

       do i=1,Nt
          zzfric(i) = 0.0
          zz(i)=yt(2*i-1)-Vpl
          if (z_all( myid*Nt+i ).lt.19.5)then
              zz(i) = 0
          end if
       end do


       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_Allgather(zz,Nt,MPI_DOUBLE_PRECISION,zz_all,Nt,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

       !----------------------------------------------------------------------
       !    summation of stiffness of all elements in slab
       !----------------------------------------------------------------------
       ! initilize zzfric
       call CPU_TIME(tm2)
       tmelse=tmelse+tm2-tm1
       tm1=tm2

    !!!! add the locked zone above 20 km for the seismogenic zone

    !       do j=1,Nt_all
    !          if(z_all(j).lt.19.5)then
    !             zz_all(j) = 0
    !          end if
    !       end do
    !!!!!!!!!!!!!!!!!!!!

       ! calculate stiffness*(Vkl-Vpl) of one proccessor.
       do i=1,Nt
          zzfric(i) = dot_product(stiff(i,1:Nt_all), zz_all(1:Nt_all))
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
          if(yt(2*i).lt.1.0d-5)then
             yt(2*i) = 1.0d-5
          end if

          if(yt(2*i-1).le.small)then 
             psi = dlog(V0*yt(2*i)/xLf(i))
             help1 = yt(2*i-1)/(2*V0)
             help2 = (f0+ccb(i)*psi)/cca(i)
             help = dsqrt(1+(help1*dexp(help2))**2)

             deriv1 = (seff(i)*ccb(i)/yt(2*i))*help1*dexp(help2)/help
             deriv2 = (seff(i)*cca(i)/(2*V0))*dexp(help2)/help
!aging             
	     deriv3 = 1 - yt(2*i-1)*yt(2*i)/xLf(i)
!slip law	     deriv3 = -yt(2*i-1)*yt(2*i)/xLf(i)*dlog(yt(2*i-1)*yt(2*i)/xLf(i))
             dydt(2*i-1)=-(zzfric(i)+deriv1*deriv3)/(eta+deriv2) 
             dydt(2*i)=deriv3     
          else
             deriv1 = seff(i)*ccb(i)/yt(2*i)
             deriv2 = seff(i)*cca(i)/yt(2*i-1)
!aging       
             deriv3 = 1 - yt(2*i-1)*yt(2*i)/xLf(i)
!slip law	     deriv3 = -yt(2*i-1)*yt(2*i)/xLf(i)*dlog(yt(2*i-1)*yt(2*i)/xLf(i))
             dydt(2*i-1) = -(zzfric(i)+deriv1*deriv3)/(eta+deriv2)
             dydt(2*i) = deriv3  
          end if
       end do

       RETURN
     END subroutine derivs

!-----------------------------------------------------------------------------
!    read parameters: sigma_effective, a,b,D_c
!----------------------------------------------------------------------------

    subroutine resdep(Nt_all,dip,hnucl,sigmadiff,sefffix,Lffix,Iperb,&
         factor1,factor2,factor3,factor4, &
         xilock1,xilock2,cca_all,ccb_all,xLf_all, &
         seff_all,x_all,z_all)
      USE mpi
      USE phy3d_module_non, only: p18,Nab,xmu,xnu,gamma, &
           Iprofile,foldername,jobname,profile
      implicit none
      integer, parameter :: DP = kind(1.0d0)
      integer, parameter :: DN=9
      integer :: k,i,j,kk,Iperb,record,l,m,nn,Nt,Nt_all

      real (DP) :: temp(DN),dep(DN),dist(DN),ptemp(Nt_all), &
           ccabmin(Nt_all),xLfmin(Nt_all),xilock1,xilock2, & 
           hnucl,sigmadiff,sefffix,Lffix,dip, &
           factor,factor1,factor2,factor3,factor4
      real (DP) :: cca_all(Nt_all),ccb_all(Nt_all),ccab_all(Nt_all), &
           xLf_all(Nt_all),seff_all(Nt_all),x_all(Nt_all),z_all(Nt_all)

      real (DP) ::a(Nab),tpr(Nab),zp(Nab),b(nab),ab(nab)

      factor = 1.0  !initial value of perturbation 

      do i = 1,Nt_all
         seff_all(i) = 28.0+180.d0*z_all(i)
         if (seff_all(i).ge.sigmadiff)  seff_all(i)=sigmadiff
      end do

      !open(444,file="sigma_seff_t5.txt",status='old')
      ! do i=1,Nt_all
      !  read(444,*) seff_all(i)
      ! end do
      !close(444)

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
         ab(4)=0.015
         ab(5)=0.025
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

      if(Iprofile.le.3)then   !use actual temp. 
         do i=1,Nab
            b(i) = a(i) - ab(i)
         end do
         open(9,file=profile,status='old')
         do k=1,DN
            read(9,*)temp(k),dep(k)
         end do
         close(9)

         do k=1,Nt_all
            if(z_all(k).ge.dep(DN))then
               i = DN-3
               j = DN
            else
               do l=1,(DN-1)
                  if(z_all(k).ge.dep(l).and.z_all(k).lt.dep(l+1)) then
                     i=l
                     j=l+1
                  end if
               end do
            end if
            ptemp(k)=temp(i)+(temp(j)-temp(i))*(z_all(k)-dep(i))/(dep(j)-dep(i))
            if(ptemp(k).ge.tpr(Nab))then
               m=Nab-1
               nn=Nab
            else
               do l=1,Nab-1
                  if (ptemp(k).ge.tpr(l).and.ptemp(k).lt.tpr(l+1)) then
                     m=l
                     nn=l+1
                  end if
               end do
            end if
            cca_all(k)=a(m)+(a(nn)-a(m))*(ptemp(k)-tpr(m))/(tpr(nn)-tpr(m))
            ccab_all(k)=ab(m)+(ab(nn)-ab(m))*(ptemp(k)-tpr(m))/(tpr(nn)-tpr(m))
            ccb_all(k)=b(m)+(b(nn)-b(m))*(ptemp(k)-tpr(m))/(tpr(nn)-tpr(m))
            xLf_all(k)=abs(ccab_all(k)**2)*seff_all(k)*hnucl*(1-xnu)/(gamma*xmu*abs(ccb_all(k)))
         end do
      end if

!!! check for minimum Dc

      do k=1,Nt_all
         xLfmin(k)=dabs(ab(4)**2)*seff_all(k)*hnucl*(1-xnu)/(gamma*xmu*abs(ccb_all(k)))
         if(xLf_all(k).lt.xLfmin(k))then        
            xLf_all(k) = xLfmin(k)
         end if
      end do

!!! set SSE depth effective normal stress and Dc
     do i=1,Nt_all
         if(z_all(i).ge.xilock1.and.z_all(i).le.xilock2)then
            seff_all(i) = sefffix
            xLf_all(i) = Lffix
         end if
      end do

!!!! add perturbation and buffer zone at both sides if necessary
!      if(Iperb.eq.1)then    !apply small perturbations in a,b
!         do k=1,Nt_all
!
!!            if(x_all(k).gt.170.0.and.x_all(k).lt.220.0)then
!!               factor = factor1
!!            else if(x_all(k).gt.-290.0.and.x_all(k).lt.-190.0)then
!!               factor = factor1
!!            else
!!               factor = 1.0
!!            end if
!
!!            if(factor.eq.factor1)then
!!               cca_all(k)=0.01
!!               ccab_all(k)=0.0035 ! bufferzone
!!               ccb_all(k)=cca_all(k)-ccab_all(k)
!!               xLf_all(k)=dabs(ccab_all(k)**2)*sigmadiff*hnucl*(1-xnu)/(gamma*xmu*abs(ccb_all(k)))
!!            else
!
!!               cca_all(k)=cca_all(k)*factor
!!               ccab_all(k)=ccab_all(k)*factor
!!               ccb_all(k)=cca_all(k)-ccab_all(k)
!!               xLf_all(k)=dabs(ccab_all(k)**2)*seff_all(k)*hnucl*(1-xnu)/(gamma*xmu*abs(ccb_all(k)))
!            end if
!         end do
!      end if


      !need to address when j=1 and j=Nd_all!! same in the old openmp f90 file!

!      do k=1,Nt_all
!         xLfmin(k)=dabs(ab(4)**2)*seff_all(k)*hnucl*(1-xnu)/(gamma*xmu*abs(ccb_all(k)))
!         if(xLf_all(k).lt.xLfmin(k))then	
!            xLf_all(k) = xLfmin(k)
!         end if
!      end do

      !open(444,file="sigma_seglf_t5.txt",status='old')
      ! do i=1,Nt_all
      !  read(444,*) xLf_all(i)
      ! end do
      !close(444)

      !open(444,file="newab2.txt",status='old')
      ! do i=1,Nt_all
      !  read(444,*) ccab_all(i)
      !  ccb_all(i) = cca_all(i) - ccab_all(i)
      ! end do
      !close(444)

      
      !      end do 

      !     To save info about some of the quantities
      open(2,file=trim(foldername)//'vardep'//jobname,status='unknown')
      !	write(2,300)'z','seff','Lf','ccab','cca'
      do i=1,Nt_all
         write(2,'(5(1x,e20.13))')z_all(i),seff_all(i),xLf_all(i), &
              ccab_all(i),cca_all(i)
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
    tmv,tas,tcos,tsse,maxv,maxv_s1,maxv_s2,maxv_s3,maxv_s4,maxv_s5,maxv_s6,&
    maxv_s7,maxv_s8,maxnum, &
     slipz1_inter,slipz1_tau,slipz1_sse,&
     slipz1_cos,&
     xi_all,x_all,slipz1_v) 


USE mpi
USE phy3d_module_non, only: xmu,nmv,nas,ncos,nsse,yrs,Vpl, &
		foldername,jobname
implicit none
integer, parameter :: DP = kind(1.0d0)
integer :: Nt,Nt_all,i,j,k,l,kk,inul,imv,ias,icos,isse,Ioutput,Isnapshot
real (DP) :: x(Nt),maxnum(nmv),maxv(nmv),maxv_s1(nmv),maxv_s2(nmv),maxv_s3(nmv),maxv_s4(nmv),maxv_s5(nmv),&
        maxv_s6(nmv),maxv_s7(nmv),maxv_s8(nmv), &
	tmv(nmv),tas(nas),tcos(ncos),tsse(nsse)

real (DP) :: slipz1_inter(Nt_all,nas),slipz1_cos(Nt_all,ncos),&
        slipz1_tau(Nt_all,nsse),slipz1_sse(Nt_all,nsse), &
     xi_all(Nt_all),x_all(Nt_all),&
      slipz1_v(Nt_all,ncos)
! integer :: n_intz1,n_intz2,n_intz3,n_cosz1,n_cosz2,n_cosz3
!integer :: intdepz1(Nt_all),intdepz2(Nt_all),intdepz3(Nt_all)

if(Ioutput == 0)then    !output during run 


   if(imv==nmv)then
      open(30,file=trim(foldername)//'maxv_all'//jobname,access='append',status='unknown')
      open(311,file=trim(foldername)//'maxv_s1'//jobname,access='append',status='unknown')
      open(312,file=trim(foldername)//'maxv_s2'//jobname,access='append',status='unknown')
      open(313,file=trim(foldername)//'maxv_s3'//jobname,access='append',status='unknown')
      open(314,file=trim(foldername)//'maxv_s4'//jobname,access='append',status='unknown')
      open(315,file=trim(foldername)//'maxv_s5'//jobname,access='append',status='unknown')
      open(316,file=trim(foldername)//'maxv_s6'//jobname,access='append',status='unknown')
      open(317,file=trim(foldername)//'maxv_s7'//jobname,access='append',status='unknown')
      open(318,file=trim(foldername)//'maxv_s8'//jobname,access='append',status='unknown')

      do i=1,nmv
         write(30,130)tmv(i),dlog10(maxv(i)*1d-3/yrs)
         write(311,130) dlog10(maxv_s1(i)*1d-3/yrs)
         write(312,130) dlog10(maxv_s2(i)*1d-3/yrs)
         write(313,130) dlog10(maxv_s3(i)*1d-3/yrs)
         write(314,130) dlog10(maxv_s4(i)*1d-3/yrs)
         write(315,130) dlog10(maxv_s5(i)*1d-3/yrs)
         write(316,130) dlog10(maxv_s6(i)*1d-3/yrs)
         write(317,130) dlog10(maxv_s7(i)*1d-3/yrs)
         write(318,130) dlog10(maxv_s8(i)*1d-3/yrs)
      end do
      close(30)
      do i=311,318
         close(i)
      end do 
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
      open(30,file=trim(foldername)//'maxv_all'//jobname,access='append',status='unknown')
      open(311,file=trim(foldername)//'maxv_s1'//jobname,access='append',status='unknown')
      open(312,file=trim(foldername)//'maxv_s2'//jobname,access='append',status='unknown')
      open(313,file=trim(foldername)//'maxv_s3'//jobname,access='append',status='unknown')
      open(314,file=trim(foldername)//'maxv_s4'//jobname,access='append',status='unknown')
      open(315,file=trim(foldername)//'maxv_s5'//jobname,access='append',status='unknown')
      open(316,file=trim(foldername)//'maxv_s6'//jobname,access='append',status='unknown')
      open(317,file=trim(foldername)//'maxv_s7'//jobname,access='append',status='unknown')
      open(318,file=trim(foldername)//'maxv_s8'//jobname,access='append',status='unknown')
     do i=1,imv
         write(30,130)tmv(i),dlog10(maxv(i)*1d-3/yrs)
         write(311,130) dlog10(maxv_s1(i)*1d-3/yrs)
         write(312,130) dlog10(maxv_s2(i)*1d-3/yrs)
         write(313,130) dlog10(maxv_s3(i)*1d-3/yrs)
         write(314,130) dlog10(maxv_s4(i)*1d-3/yrs)
         write(315,130) dlog10(maxv_s5(i)*1d-3/yrs)
         write(316,130) dlog10(maxv_s6(i)*1d-3/yrs)
         write(317,130) dlog10(maxv_s7(i)*1d-3/yrs)
         write(318,130) dlog10(maxv_s8(i)*1d-3/yrs)
      end do
      close(30)
      do i=311,315
         close(i)
      end do 
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


end if 

 120    format(E20.13,4X,E20.13,4X,I6)
 130    format(E20.13,2(1X,E13.6))
 140    format(E20.13)
 150    format(E15.8,3(1X,E15.8))
 160    format(E20.13,1x,E20.13)
 500    format(E15.8,1X,E20.13)
 600    format(E15.4,1X,E13.6,1X,E15.8)
 700    format(E13.6)
 900    format(E15.8)

RETURN
END subroutine output
