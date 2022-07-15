!------------------------------------------------
! 3D Subduction fault earthquake sequence model
! MPI 
! Yajing Liu
! Last modified Jul. 4, 2012
!
! LZ
! Adding normal stress change(const dc)&changing SV
! Oct.24,2018
!------------------------------------------------

program main
USE mpi
USE phy3d_module_non
implicit none
integer, parameter :: DP=kind(1.d0)

logical cyclecont
integer ::Nt_all,Nt, ndt,ndtnext,kk,ii,n,l,ndt_v,ndt_inter,Itout,&  
        ndtmax,i,j,k,nm,nrec,ihrestart,Ifileout,Iperb,record,Isnapshot,&
        ix1,ix2,ix3,ix4,iz1,iz2,iz3,&
        record_cor,n_cell_dummy,s1,s2,s3,s4,s5,s6,s7,s8

real (DP) :: accuracy,areatot, epsv,dt_try, dt,dtmin,dt_did,dt_next,dip, &
     hnucl,hl,hd, hstarfactor, &
     sigmadiff,sefffix,Lffix,xdepth,xlength,&
     t,tprint_inter, tint_out,tout,&
     tmin_out,tint_cos,&   
     tslip_ave,tslipend,tslip_aveint, tmax, &
     tslipcos,tstart1,tend1,tstart2,tend2,tstart3,tend3, &
     tssestart,tsseend, &
     factor1,factor2,factor3,factor4,vtmp,fr_dummy,&
     xilock1,xilock2,x1,x2,x3,x4,z1,z2,z3,&
     xi_dummy,x_dummy,z_dummy,stiff_dummy,help

real (DP) ::  tmbegin,tmrun
   


real (DP), DIMENSION(:), ALLOCATABLE :: x,z,xi,yt,dydt,yt_scale,tau,tau0, &
     slip,slipinc,zzseff

!Arrays only defined at master cpu
real (DP), DIMENSION(:), ALLOCATABLE :: x_all,xi_all,z_all,&
     yt_all,dydt_all,yt_scale_all,tau_all, &
     slip_all,slipinc_all,cca_all,ccb_all,xLf_all,seff_all


!output related parameters
integer :: imv,ias,icos,Ioutput
real (DP) :: vcos
real (DP), DIMENSION(:), ALLOCATABLE :: maxv,tmv,tas,tcos
real (DP), DIMENSION(:,:), ALLOCATABLE :: slipz1_inter,slipz1_cos,v_cos,slip_cos,&
      slipz1_tau,slipz1_v,slipz1_seff,seff_inter,xLf_inter,tau_inter,v_inter


character(len=40) :: cTemp,filename,ct

!MPI RELATED DEFINITIONS
integer :: ierr,size,myid,master

call MPI_Init(ierr)
CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
CALL MPI_COMM_SIZE( MPI_COMM_WORLD, size, ierr )

master = 0 

!read in intiliazation parameters

open(12,file='parameter_ns_const_dc.txt',status='old')

      read(12,110)jobname
      read(12,110)foldername
      read(12,110)stiffname
      read(12,110)restartname
      read(12,110)profile
!!!modified data read in
      read(12,*)Nab,Nt_all,Nt,Lratio,nprocs
      read(12,*)IDin,IDout,Iprofile,Iperb,Isnapshot 
      read(12,*)xlength,xdepth,dip,Vpl
      read(12,*)hnucl
      read(12,*)sigmadiff,sefffix,Lffix
      read(12,*)tmax
      read(12,*)tslip_ave,tslipend,tslip_aveint
      read(12,*)tint_out,tmin_out,tint_cos
      read(12,*)vcos
      read(12,*)nmv,nas,ncos

!!! modified data read in
 110  format(A)
      close(12)


if(mod(Nt_all,nprocs)/=0)then
   write(*,*)'Nd_all must be integer*nprocs. Change nprocs!'
   STOP

end if


if(myid == master)then
   ALLOCATE(x_all(Nt_all),xi_all(Nt_all),z_all(Nt_all),&
   cca_all(Nt_all),ccb_all(Nt_all),seff_all(Nt_all),xLf_all(Nt_all),&
   tau_all(Nt_all),slipinc_all(Nt_all),yt_all(2*Nt_all), &
   dydt_all(2*Nt_all),yt_scale_all(2*Nt_all))

   ALLOCATE (maxv(nmv),tmv(nmv),tas(nas),tcos(ncos))

!!! modify output number
   ALLOCATE (slipz1_inter(Nt_all,nas),slipz1_cos(Nt_all,ncos), &
	           v_cos(Nt_all,ncos),slip_cos(Nt_all,ncos), &
             slipz1_tau(Nt_all,ncos),slipz1_v(Nt_all,ncos),seff_inter(Nt_all,nas),&
            slipz1_seff(Nt_all,ncos),xLf_inter(Nt_all,nas),tau_inter(Nt_all,nas),&
             v_inter(Nt_all,nas))
   
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ALLOCATE (x(Nt),z(Nt),xi(Nt),cca(Nt),ccb(Nt),seff(Nt),&
        xLf(Nt),tau(Nt),tau0(Nt),slip(nt),slipinc(nt), &
         yt(2*Nt),dydt(2*Nt),yt_scale(2*Nt),zzseff(Nt))

ALLOCATE (slip_all(Nt_all))   !!!allocate slip_all in all proces

ALLOCATE (stiff(Nt,Nt_all))   !!! stiffness of Stuart green calculation

ALLOCATE (stiff_ns(Nt,Nt_all))

!Read in stiffness matrix, in nprocs segments

write(cTemp,*) myid
!write(*,*) cTemp
!write(*,*) stiffname

open(5, file=trim(stiffname)//'trigreen_'//trim(adjustl(cTemp))//'.bin',form='unformatted',access='stream')
open(51,file=trim(stiffname)//'trigreen_ns_'//trim(adjustl(cTemp))//'.bin',form='unformatted',access='stream')

open(55,file=trim(stiffname)//'position.bin',form='unformatted',access='stream')
! record = Nt_all*Nt*(myid)

record_cor=Nt*(myid)

!-------------------------------------------------------------------------------------------
! read stiffness from Stuart green calculation.
!-----------------------------------------------------------------------------------------

i=0
do while (i<record_cor)
 read(55) xi_dummy,x_dummy,z_dummy
 i=i+1
end do
write(*,*) 'i=',i
open(342,file='cell_position',access='append',status='unknown')
do k=1,Nt
   read(55) xi(k),x(k),z(k) !xi is along the dip(y) while x is along the strike
   xi(k)=xi(k)/1000  ! unit change from m to km
   x(k)=x(k)/1000
   z(k)=z(k)/1000
   write(342,*) xi(k), x(k), z(k)
end do
  close(342)  


do i=1,Nt !! observe
do j=1,Nt_all !! source
 read(5) stiff(i,j)
 read(51) stiff_ns(i,j)
end do
end do
 close(5)
 close(51)
 close(55)
!!-----------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------


call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_Gather(x,Nt,MPI_Real8,x_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
call MPI_Gather(xi,Nt,MPI_Real8,xi_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
call MPI_Gather(z,Nt,MPI_Real8,z_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)


if(myid==master)then
call resdep(Nt_all,dip,hnucl,sigmadiff,sefffix,Lffix, &  
           cca_all,ccb_all,xLf_all,seff_all,xi_all,z_all,x_all)

end if


call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_Scatter(cca_all,Nt,MPI_Real8,cca,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
call MPI_Scatter(ccb_all,Nt,MPI_Real8,ccb,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
call MPI_Scatter(xLf_all,Nt,MPI_Real8,xLf,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
call MPI_Scatter(seff_all,Nt,MPI_Real8,seff,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)



call CPU_TIME(tmbegin)

tm1=tmbegin
tmday=86400.0
tmelse=0.0
tmmult=0.0
tmmidn=0.0
  
t=0.d0
tprint_inter = 0.d0
tslipcos = 0.d0

tout =0.d0

ndtnext=0
ndt_v=0
ndt_inter=0
ndt =0
nrec=0  

Ioutput = 0    !initial always 0 (output)

imv=0   !counter for maxv, sliptot, slipthresh1, slipthresh2,slipthresh3 output 
ias=0 !counter for slip at iz3 and s.z. average slip output 
icos = 0 

accuracy = 1.d-6
epsv = 1.d-4 
dtmin = 1.d-10
dt_try=dtmin


if(myid==master)then                
  open(1,file=trim(foldername)//'phypara'//jobname, status='unknown')
end if

Ifileout = 60   !file index, after 47

!----Initial values of velocity, state variable, shear stress and slip--
!--SET INITIAL VPL FOR THE LOCKED PART TO BE 0 
!open(343,file='yt_value',access='append',status='unknown')
   if(IDin.eq.0)then 
      do j=1,Nt
!        yt(2*j-1)=Vpl
!------- planar model 0-270km
!         yt(2*j-1)=(1.536d-7)*((x(j)-93)**3)-(2.777d-5)*((x(j)-93)**2)-0.00643*(x(j)-93)+1.967 
!--------lms fault model -93-180km
         yt(2*j-1)=(1.536d-7)*(x(j)**3)-(2.777d-5)*(x(j)**2)-0.00643*x(j)+1.967
	     yt(2*j)=xLf(j)/(1.1*yt(2*j-1))	
         slip(j)=0.d0
!         write(343,*) yt(2*j-1), yt(2*j)
      end do
    end if 
!close(343)


!------------------------------------------------------------------
  if(IDin.eq.1) then 
!if this is a restart job
   if(myid==master)then
   call restart(0,'out',4,Nt_all,t,dt,dt_try,ndt,nrec,yt_all,slip_all,seff_all,xLf_all)
   write(1,*)'This is a restart job. Start time ',t,' yr'
   end if
   
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   call MPI_Bcast(t,1,MPI_Real8,master,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(dt,1,MPI_Real8,master,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(dt_try,1,MPI_Real8,master,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ndt,1,MPI_integer,master,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nrec,1,MPI_integer,master,MPI_COMM_WORLD,ierr)
   call MPI_Scatter(yt_all,2*Nt,MPI_Real8,yt,2*Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
   call MPI_Scatter(slip_all,Nt,MPI_Real8,slip,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
   call MPI_Scatter(seff_all,Nt,MPI_Real8,seff,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
   call MPI_Scatter(xLf_all,Nt,MPI_Real8,xLf,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
   
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
!
   call derivs(myid,dydt,2*Nt,Nt_all,Nt,t,yt,x,zzseff) 

        do j=1,2*Nt
            yt_scale(j)=dabs(yt(j))+dabs(dt_try*dydt(j))
        end do
        do i=1,Nt
            seff(i)=seff(i)+dt_try*zzseff(i)
        end do

   call rkqs(myid,yt,dydt,2*Nt,Nt_all,Nt,t,dt_try,accuracy,yt_scale, &
            dt_did,dt_next,x,zzseff)

        dt = dt_did
        dt_try = dt_next

          do i=1,Nt
              if(yt(2*i-1).lt.epsv)then
                 help=(yt(2*i-1)/(2*V0))*dexp((f0+ccb(i)*dlog(V0*yt(2*i)/xLf(i)))/cca(i))
                 tau(i)=seff(i)*cca(i)*dlog(help+dsqrt(1+help**2))
              else
                 tau(i)=seff(i)*(f0+cca(i)*dlog(yt(2*i-1)/v0)+ccb(i)*dlog(V0*yt(2*i)/xLf(i)))
              end if
             slipinc(i) = yt(2*i-1)*dt
             slip(i) = slip(i) + slipinc(i)
          end do
          ndt = ndt + 1

call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_Gather(yt,2*Nt,MPI_Real8,yt_all,2*Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
call MPI_Gather(slipinc,Nt,MPI_Real8,slipinc_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
call MPI_Gather(slip,Nt,MPI_Real8,slip_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
call MPI_Gather(tau,Nt,MPI_Real8,tau_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)
call MPI_Gather(seff,Nt,MPI_Real8,seff_all,Nt,MPI_Real8,master,MPI_COMM_WORLD,ierr)

!-------------------
!      Output:     (a single thread will do the writing while others 
!      proceed 
!-------------------

if(myid==master)then
   imv=imv+1
   tmv(imv)=t
   maxv(imv) = 0.d0


   do i=1,Nt_all      !!!!!!!!!!!!!!!!! find max velocity
      if(yt_all(2*i-1).ge.maxv(imv))then
	       maxv(imv)=yt_all(2*i-1)
      end if
   end do

!-----Interseismic slip every ? years----

   if (t.ge.tslip_ave)then
         ias = ias + 1 
         tas(ias)=t

      do i=1,Nt_all
        slipz1_inter(i,ias)=slip_all(i)*1.d-3
        seff_inter(i,ias) =seff_all(i)
        tau_inter(i,ias)=tau_all(i)
        v_inter(i,ias)=dlog10(yt_all(2*i-1)*1.d-3/yrs)
        xLf_inter(i,ias) =xLf_all(i)
      end do
       tslip_ave = tslip_ave + tslip_aveint
   end if

!------coseismic Slip   -------------------

  if((maxv(imv)/yrs).ge.vcos)then
      tslipcos = tslipcos+dt
    if(tslipcos.ge.tint_cos)then    !!!!tint_cos: every 5s output 
      icos = icos +1
      tcos(icos) = t 
  
!--------------------------------------------------------------
         do i=1,Nt_all
            slipz1_tau(i,icos)=tau_all(i)
            slipz1_v(i,icos)=dlog10(yt_all(2*i-1)*1.d-3/yrs)
            slipz1_cos(i,icos)=slip_all(i)*1.d-3
            slipz1_seff(i,icos)=seff_all(i)
         end do
!--------------------------------------------------------------
         if(Isnapshot == 1)then
            do k=1,Nt_all
              v_cos(k,icos)=yt_all(2*k-1)
              slip_cos(k,icos)=slip_all(k)*(1.d-3)
            end do 
         end if 
        tslipcos = 0.d0
     end if
   end if 
 
end if

!----Output restart files -------------
  if(myid==master)then
          if(mod(ndt,500).eq.0)ihrestart=1

          if(IDout.eq.1.and.ihrestart.eq.1)then
	          filename='out0'
           call restart(1,filename,4,Nt_all,t,dt,dt_try,ndt,nrec,yt_all,slip_all,seff_all,xLf_all)  
           ihrestart=0
          end if

          if(abs(t-tout).le.tmin_out)then
	           ihrestart = 1
	           Itout=int(tout)
	           write(ct,*)Itout
	           ct=adjustl(ct)
	           filename='out'//trim(ct)
             call restart(1,filename,Ifileout,Nt_all,t,dt,dt_try,ndt,nrec,yt_all,slip_all,seff_all,xLf_all)
             ihrestart = 0
             tout = tout+tint_out 
          end if

       if (dlog10(maxv(imv)/Vpl) .gt. 10) then
          open(34,file=trim(foldername)//'maxv_infinite'//jobname,access='append',status='unknown')
          do i=1,Nt_all
            write(34,*) yt_all(i*2-1)
          end do
         close(34)
       end if
  end if

!----Output velocity and slip records --- 
!----velocity in mm/yr ; slip in meter, moment in 10^{14} Nm -- 
  if(myid==master)then 
    Ioutput = 0 
    call output(Ioutput,Isnapshot,Nt_all,Nt,hd,imv,ias,icos,x,&
     tmv,tas,tcos,maxv,slipz1_inter,slipz1_cos,slip_cos,v_cos,&
     xi_all,x_all,slipz1_tau,slipz1_v,seff_inter,slipz1_seff,xLf_inter,tau_inter,v_inter)         

  end if

!      end of output 
!----------------------------------------------------
!   End of output commands.
!----------------------------------------------------  

!----------------------------------------------------
!   Go to next cycle if t < tmax or ndt > ndtmax 
!---------------------------------------------------

   if (t>tmax)cyclecont = .false.

end do

!--- Final output ------- 

if(myid==master)then 
   filename='outlast'
   call restart(1,filename,Ifileout,Nt_all,t,dt,dt_try,ndt,nrec,yt_all,slip_all,seff_all,xLf_all)
    Ioutput = 1
   call output(Ioutput,Isnapshot,Nt_all,Nt,hd,imv,ias,icos,x,&
    tmv,tas,tcos,maxv,slipz1_inter,slipz1_cos,slip_cos,v_cos,&
     xi_all,x_all,slipz1_tau,slipz1_v,seff_inter,slipz1_seff,xLf_inter,tau_inter,v_inter) 

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
DEALLOCATE (x_all,xi_all,z_all,yt_all,dydt_all,yt_scale_all,tau_all, &
     slip_all,slipinc_all,cca_all,ccb_all,xLf_all,seff_all, &
     maxv,tmv,tcos,tas)

DEALLOCATE (slipz1_inter,slipz1_cos,v_cos,slip_cos,slipz1_tau,slipz1_v,seff_inter,&
slipz1_seff,xLf_inter,tau_inter,v_inter)
end if 


DEALLOCATE (stiff,stiff_ns,zzseff,xLf)
DEALLOCATE (x,z,xi,yt,dydt,yt_scale,tau,tau0,slip,slipinc)

call MPI_finalize(ierr)

END

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine rkqs(myid,y,dydx,n,Nt_all,Nt,x,htry,eps,yscal,hdid,hnext,xx,zzseff)
Use mpi
USE phy3d_module_non, only : nprocs
implicit none
      integer, parameter :: DP = kind(1.0d0)   
       integer :: n,i,j,k,NMAX,Nt,Nt_all
       real (DP) :: eps,hdid,hnext,htry,x
	real (DP) :: dydx(n),y(n),yscal(n),xx(Nt),zzseff(Nt)
       external derivs
       real (DP) :: errmax,errmax1,h,htemp,xnew,errmax_all(nprocs)
	real (DP), dimension(:), allocatable :: yerr,ytemp
      real (DP), parameter :: SAFETY=0.90, PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4

!MPI RELATED DEFINITIONS
integer :: ierr,myid,master
master = 0 

      nmax=2*Nt
      h=htry
      allocate (yerr(nmax),ytemp(nmax))
!     write(*,*)'in rkqs, bf rkck, myid=',myid,y(n-1)
  1    call rkck(myid,dydx,h,n,Nt_all,Nt,y,yerr,ytemp,x,derivs,xx,zzseff)
!        write(*,*)'in rkqs, af rkck, myid=',myid, ytemp(n-1)
       errmax=0.
       do  i=1,nmax
             errmax=max(errmax,dabs(yerr(i)/yscal(i)))
       end do 
       errmax=errmax/eps
call MPI_Barrier(MPI_COMM_WORLD,ierr)

call MPI_Gather(errmax,1,MPI_Real8,errmax_all,1,MPI_Real8,master,MPI_COMM_WORLD,ierr)

   if(myid==master)then
      errmax1=maxval(errmax_all)
   end if
CALL MPI_BCAST(errmax1,1,MPI_REAL8,master,MPI_COMM_WORLD, ierr)

       if(errmax1.gt.1.)then
            htemp=SAFETY*h*(errmax1**PSHRNK)
            h=dsign(max(dabs(htemp),0.1*dabs(h)),h)
            xnew=x+h
           ! if(xnew.eq.x)pause 'stepsize underflow in rkqs'
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
       end
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine rkck(myid,dydx,h,n,Nt_all,Nt,y,yerr,yout,x,derivs,xx,zzseff)
USE phy3d_module_non, only :nprocs
implicit none
integer, parameter :: DP = kind(1.0d0)   
integer :: n,i,NMAX,myid,Nt_all,Nt
external derivs
real (DP) :: h,x,dydx(n),y(n),yerr(n),yout(n),xx(Nt),zzseff(Nt)
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
       call derivs(myid,ak2,n,Nt_all,Nt,x+A2*h,ytemp,xx,zzseff)

       do  i=1,n
            ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
       end do 
      call derivs(myid,ak3,n,Nt_all,Nt,x+A3*h,ytemp,xx,zzseff)

       do  i=1,n
            ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
       end do 
       call derivs(myid,ak4,n,Nt_all,Nt,x+A4*h,ytemp,xx,zzseff)

       do  i=1,n
            ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+  &
                    B54*ak4(i))
       end do 
       call derivs(myid,ak5,n,Nt_all,Nt,x+A5*h,ytemp,xx,zzseff)

       do  i=1,n
            ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+  &
                  B64*ak4(i)+B65*ak5(i))
       end do 
       call derivs(myid,ak6,n,Nt_all,Nt,x+A6*h,ytemp,xx,zzseff)

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
  end
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine derivs(myid,dydt,nv,Nt_all,Nt,t,yt,x,zz_seff)
USE mpi
USE phy3d_module_non, only: stiff,cca,ccb,seff,xLf,eta,f0,Vpl,V0,Lratio,nprocs,&
                        tm1,tm2,tmday,tmelse,tmmidn,tmmult,stiff_ns
       implicit none
      integer, parameter :: DP = kind(1.0d0)
       integer :: nv,n,i,j,k,kk,l,ii,Nt,Nt_all
       real (DP) :: t,yt(nv),dydt(nv)   
       real (DP) :: deriv3,deriv2,deriv1,small
       real (DP) :: psi,help1,help2,help,dd
       real (DP) :: SECNDS
       real (DP) :: zz(Nt),zzfric(Nt),zzfric_all(Nt_all), zz_all(Nt_all),&
                    x(Nt),zz_seff(Nt)
       intrinsic imag,real

!MPI RELATED DEFINITIONS
integer :: ierr,myid,master
master = 0 

!      write(*,*) myid
       small=1.d-4
!      open(341,file='zz(i)',access='append',status='unknown')
        do i=1,Nt
!          zz(i)=yt(2*i-1)-Vpl 
!----------planar model 0-270km
!          dd=(1.536d-7)*((x(i)-93)**3)-(2.777d-5)*((x(i)-93)**2)-0.00643*(x(i)-93)+1.967-0.001
!----------lms fault model -93-180km
           dd=(1.536d-7)*(x(i)**3)-(2.777d-5)*(x(i)**2)-0.00643*x(i)+1.967-0.001
          zz(i)=yt(2*i-1)-dd
!          write(341,'(E20.13)') zz(i)
        end do
!     close(341)

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

 do j=1,Nt
  zzfric(j)=0
  zz_seff(j)=0
 end do

! calculate stiffness*(Vkl-Vpl) of one proccessor.
 do i=1,Nt
   do j=1, Nt_all
   ! number in all element, total Nt_all
    zzfric(i)=zzfric(i) + stiff(i,j)*zz_all(j)
    zz_seff(i)=zz_seff(i) + stiff_ns(i,j)*zz_all(j)
   end do
 end do

 call CPU_TIME(tm2)

  tmmult=tmmult+tm2-tm1

  tm1=tm2


      do i=1,Nt
          if (yt(2*i-1).le.1e-5) then
              yt(2*i-1)=yt(2*i-1)+1e-5
          end if
               yt(2*i)=abs(yt(2*i))

           if(yt(2*i-1).le.small)then 
              psi = dlog(V0*yt(2*i)/xLf(i))
              help1 = yt(2*i-1)/(2*V0)
              help2 = (f0+ccb(i)*psi)/cca(i)
              help = dsqrt(1+(help1*dexp(help2))**2)
              
              deriv1 = (seff(i)*ccb(i)/yt(2*i))*help1*dexp(help2)/help
              deriv2 = (seff(i)*cca(i)/(2*V0))*dexp(help2)/help
              deriv3 = 1 - yt(2*i-1)*yt(2*i)/xLf(i)
              dydt(2*i-1)=-(zzfric(i)+deriv1*deriv3)/(eta+deriv2) 
              dydt(2*i)=deriv3     
            else
              deriv1 = seff(i)*ccb(i)/yt(2*i)
              deriv2 = seff(i)*cca(i)/yt(2*i-1)
              deriv3 = 1 - yt(2*i-1)*yt(2*i)/xLf(i)
              dydt(2*i-1) = -(zzfric(i)+deriv1*deriv3)/(eta+deriv2)
              dydt(2*i) = deriv3  
            end if
      end do

!DEALLOCATE  (zz,zzfric,data,data1,data2,fft1,fft2,zr,zi,cr,ci)

      RETURN
 END

!-----------------------------------------------------------------------------
!    FOUR1
!----------------------------------------------------------------------------
! --- End of subroutine FOUR1 -------------------------------------------

! ---- Subroutine TWOFFT to do FFT of two real data streams ------------
! ----from Press et al., Numerical Recipes ---- calls FOUR1 above -------
! --- End of subroutine TWOFFT --------------------------------------

subroutine resdep(Nt_all,dip,hnucl,sigmadiff,sefffix,Lffix, &
                  cca_all,ccb_all,xLf_all, &
                  seff_all,xi_all,z_all,x_all)
USE mpi
USE phy3d_module_non, only: p18,Nl,Nd_all,Nd,Nab,xmu,xnu,gamma, &
     Iprofile,foldername,jobname,profile
implicit none
integer, parameter :: DP = kind(1.0d0)
integer, parameter :: DN=35 !4020 !34 
integer :: k,i,j,kk,Iperb,record,l,m,nn,Nt,Nt_all

real (DP) :: temp(DN),dep(DN),dist(DN),ptemp(Nt_all), &
             ccabmin(Nl),xLfmin(Nl),hnucl,hl,hd,sigmadiff,sefffix,Lffix,dip
             
real (DP) :: cca_all(Nt_all),ccb_all(Nt_all),ccab_all(Nt_all), &
       xLf_all(Nt_all),seff_all(Nt_all),xi_all(Nt_all),z_all(Nt_all),&
       x_all(Nt_all)

real (DP) ::a(Nab),tpr(Nab),zp(Nab),b(nab),ab(nab)

do i=1,Nt_all
   seff_all(i)=28.0+180.d0*z_all(i)
   if (seff_all(i).ge.sigmadiff)seff_all(i)=sigmadiff
   if (x_all(i).ge.-44.0 .and. x_all(i).le.-41.0 .and. z_all(i).lt.3.0) then  !set locked zone of XYD area
      seff_all(i) = 5000     ! normal stress 500MPa
   end if

end do 


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
      if(Iprofile.eq.1)then   !wet granite  used in Liu&Rice(2009)
         tpr(1)=0
         tpr(2)=100  ! 8km in depth
         tpr(3)=380  !350->380;
         tpr(4)=450
         tpr(5)=540  !500->540;
!	a(1)=0.01365   
!	a(2)=0.01865   
!	a(3)=0.03115   
!	a(4)=0.03615  
!	a(5)=0.03865
         a(1) = 0.01 !0.015->0.010
         a(2) = 0.01
         a(3) = 0.01
         a(4) = 0.01
         a(5) = 0.02
	ab(1)=0.008 !0.004->0.008
	ab(2)=-0.006 !-0.004 to -0.006 
	ab(3)=-0.006
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
               do l=1,Nab-1
                  if (ptemp(k).ge.tpr(l).and.ptemp(k).lt.tpr(l+1)) then
                     m=l
                     nn=l+1
                  end if
               end do
                   if (ptemp(k).ge.tpr(Nab)) then
                     m=4
                     nn=5
                  end if
                  cca_all(k)=a(m)+(a(nn)-a(m))*(ptemp(k)-tpr(m))/(tpr(nn)-tpr(m))
                  ccab_all(k)=ab(m)+(ab(nn)-ab(m))*(ptemp(k)-tpr(m))/(tpr(nn)-tpr(m))
                  ccb_all(k)=b(m)+(b(nn)-b(m))*(ptemp(k)-tpr(m))/(tpr(nn)-tpr(m))
               
               xLf_all(k)=abs(ab(3))*seff_all(k)*hnucl*(1-xnu)/(gamma*xmu)
            end do
       end if
	
!--------buffer zone--------------------
!     do i=1,Nt_all
!         if (x_all(i).gt.-100.and.x_all(i).lt.-90.0) then  !set the vs zone of longmenshan fault
!            cca_all(i)=a(1);
!            ccab_all(i)=ab(1);
!            ccb_all(i)=b(1);
!         end if
!     end do
   !---------------------------------------
!     To save info about some of the quantities
      open(2,file=trim(foldername)//'vardep'//jobname,status='unknown')
      write(2,300)'z','seff','Lf','ccab','cca'
      do i=1,Nt_all
         write(2,'(5(1x,e20.13))')z_all(i),seff_all(i),xLf_all(i), &
                    ccab_all(i),cca_all(i)
      end do 
      close(2)
 300  format(5(1x,A20))
      RETURN
      END
!       
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

 subroutine restart(inout,filename,Ifileout,Nt_all,t,dt,dt_try,ndt,nrec,yt,slip,seff,xLf)
      USE phy3d_module_non, ONLY : jobname,foldername,restartname, &
                        tm1,tm2,tmday,tmelse,tmmidn,tmmult
      implicit none
      integer, parameter :: DP = kind(1.0d0)
      integer :: inout,i,ndt,nrec,Ifileout,Nt,Nt_all
      real (DP) :: t,dt,dt_try
      real (DP) ::  yt(2*Nt_all),slip(Nt_all),seff(Nt_all),xLf(Nt_all)
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
          do i=1,Nt_all
             read(Ifileout,*)seff(i)
          end do
          do i=1,Nt_all
             read(Ifileout,*)xLf(i)
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
         do i=1,Nt_all
              write(Ifileout,*)seff(i)
         end do
         do i=1,Nt_all
             write(Ifileout,*)xLf(i)
         end do
         close(Ifileout)
     end if

      RETURN
 END

!------Output -------------------------------------------
!--------------------------------------------------------
subroutine output(Ioutput,Isnapshot,Nt_all,Nt,hd,imv,ias,icos,x,&
    tmv,tas,tcos,maxv,slipz1_inter,slipz1_cos,slip_cos,v_cos,&
     xi_all,x_all,slipz1_tau,slipz1_v,seff_inter,slipz1_seff,xLf_inter,&
     tau_inter,v_inter) 

USE mpi
USE phy3d_module_non, only: xmu,nmv,nas,ncos,nnul,nsse,Vpl,Nd_all,Nl, &
		foldername,jobname
implicit none
integer, parameter :: DP = kind(1.0d0)
integer :: Nt,Nt_all,i,j,k,l,kk,imv,ias,icos,Ioutput,Isnapshot
real (DP) :: hl,hd
real (DP) :: x(Nt),maxv(nmv),tmv(nmv),tas(nas),tcos(ncos)

real (DP) :: slipz1_inter(Nt_all,nas),slipz1_cos(Nt_all,ncos),&
        v_cos(Nt_all,ncos),slip_cos(Nt_all,ncos),slipz1_tau(Nt_all,ncos),slipz1_v(Nt_all,ncos), &
        xi_all(Nt_all),x_all(Nt_all),slipz1_seff(Nt_all,ncos),&
        seff_inter(Nt_all,nas),xLf_inter(Nt_all,nas),tau_inter(Nt_all,nas),v_inter(Nt_all,nas)


if(Ioutput == 0)then    !output during run 


   if(imv==nmv)then
      open(30,file=trim(foldername)//'maxv'//jobname,access='append',status='unknown')
      !open(311,file=trim(foldername)//'errmax'//jobname,access='append',status='unknown')
     
      do i=1,nmv
         write(30,130)tmv(i),dlog10(maxv(i)/Vpl)
      end do
      close(30)
      !close(311)
       imv = 0
   end if

    if(ias==nas)then
       open(31,form='unformatted',file=trim(foldername)//'slipz1-inter'//jobname,access='append',status='unknown')
       open(32,file=trim(foldername)//'t-inter'//jobname,access='append',status='unknown')
       open(39,form='unformatted',file=trim(foldername)//'seff-inter'//jobname,access='append',status='unknown')
      ! open(50,form='unformatted',file=trim(foldername)//'xLf'//jobname,access='append',status='unknown')
       open(51,form='unformatted',file=trim(foldername)//'tau-inter'//jobname,access='append',status='unknown')
       open(52,form='unformatted',file=trim(foldername)//'v-inter'//jobname,access='append',status='unknown')
      do j=1,nas
        do i=1,Nt_all
        write(31) slipz1_inter(i,j)
        write(39) seff_inter(i,j)
        !write(50) xLf_inter(i,j)
        write(51) tau_inter(i,j)
        write(52) v_inter(i,j)
        end do
      end do
     
      do i=1,nas
      write(32,*) tas(i)
      end do

       do i=31,32
         close(i)
       end do 
       close(39)
     !  close(50)
       close(51)
       close(52)
       ias = 0 
    end if

    if(icos==ncos)then
        if(Isnapshot == 1)then
           open(40,file=trim(foldername)//'vs-cos'//jobname,access='append',status='unknown')
           do j=1,ncos
                 do kk=1,Nt_all
                    write(40,160)v_cos(kk,j),slip_cos(kk,j)
                 end do 
           end do 
          close(40)
        end if

        open(45,form='unformatted',file=trim(foldername)//'slipz1-cos'//jobname,access='append',status='unknown')
        open(46,file=trim(foldername)//'t-cos'//jobname,access='append',status='unknown')
        open(47,form='unformatted',file=trim(foldername)//'slipz1-v'//jobname,access='append',status='unknown')
        open(48,form='unformatted',file=trim(foldername)//'slipz1-tau'//jobname,access='append',status='unknown')
        open(49,form='unformatted',file=trim(foldername)//'slipz1-seff'//jobname,access='append',status='unknown')


       do j=1, ncos
         do i=1, Nt_all
           write(45) slipz1_cos(i,j)
           write(47) slipz1_v(i,j)
           write(48) slipz1_tau(i,j)
           write(49) slipz1_seff(i,j)
         end do
       end do

       do j=1,ncos
         write(46,*) tcos(j)
       end do
       do i=45,49
         close(i)
       end do 
       icos = 0 

   end if

else 

   if((imv>0).and.(imv<nmv))then
      open(30,file=trim(foldername)//'maxv'//jobname,access='append',status='unknown')
      do i=1,imv
         write(30,130)tmv(i),dlog10(maxv(i)/Vpl)
      end do
      close(30)
      imv = 0
   end if

   if((ias>0).and.(ias<nas))then
        open(31,form='unformatted',file=trim(foldername)//'slipz1-inter'//jobname,access='append',status='unknown')
        open(32,file=trim(foldername)//'t-inter'//jobname,access='append',status='unknown')
        open(39,form='unformatted',file=trim(foldername)//'seff-inter'//jobname,access='append',status='unknown')
       ! open(50,form='unformatted',file=trim(foldername)//'xLf'//jobname,access='append',status='unknown')
        open(51,form='unformatted',file=trim(foldername)//'tau-inter'//jobname,access='append',status='unknown')
        open(52,form='unformatted',file=trim(foldername)//'v-inter'//jobname,access='append',status='unknown')
     do j=1,ias
       do i=1, Nt_all
         write(31) slipz1_inter(i,j)
         write(39) seff_inter(i,j)
        ! write(50) xLf_inter(i,j)
         write(51) tau_inter(i,j)
         write(52) v_inter(i,j)
       end do
     end do

     do i=1,ias
        write(32,*) tas(i)
     end do
     close(31)
     close(32) 
     close(39)
    ! close(50)
     close(51)
     close(52)
       ias = 0 
   end if


   if((icos>0).and.(icos<ncos))then
        open(45,form='unformatted',file=trim(foldername)//'slipz1-cos'//jobname,access='append',status='unknown')
        open(46,file=trim(foldername)//'t-cos'//jobname,access='append',status='unknown')
        open(47,form='unformatted',file=trim(foldername)//'slipz1-v'//jobname,access='append',status='unknown')
        open(48,form='unformatted',file=trim(foldername)//'slipz1-tau'//jobname,access='append',status='unknown')
        open(49,form='unformatted',file=trim(foldername)//'slipz1-seff'//jobname,access='append',status='unknown')

       do j=1,icos
          do i=1,Nt_all
            write(45) slipz1_cos(i,j)
            write(47) slipz1_v(i,j)
            write(48) slipz1_tau(i,j)
            write(49) slipz1_seff(i,j)
          end do
       end do

       do i=1,icos
         write(46,*) tcos(i)
       end do

       do i=45,49
         close(i)
       end do 

        if(Isnapshot == 1)then
           open(40,file=trim(foldername)//'vs-cos'//jobname,access='append',status='unknown')
           do j=1,icos
              do kk=1,Nt_all
                 write(40,160)v_cos(kk,j),slip_cos(kk,j)
              end do 
           end do 
          close(40)
        end if
         icos = 0 
    end if   
            
end if 

 120    format(E20.13,4X,E20.13,4X,I6)
 130    format(E20.13,1X,E13.6)
 140    format(E20.13)
 150    format(E15.8,3(1X,E15.8))
 160    format(E20.13,1x,E20.13)
 500    format(E15.8,1X,E20.13)
 600    format(E15.4,1X,E13.6,1X,E15.8)
 700    format(E13.6)
 900    format(E15.8)

RETURN
END 
