!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module mod_dtrigreen
  implicit none
  
  ! Module variables (replacing common blocks)
  integer :: nwarn = 0
  
  contains
  
subroutine dstuart (nu, xo, tridlc, ss,ds,op, u,t)

!     Calls Stuart's angle subroutine to get triangle dislocation results.

!     Comninou and Dundurs use a 1,2,3 = NED coordinate system.
!       Their simple angular dislocation is in the 1,3 plane
!       with its tip under the origin and its arm pointing toward +x1.
  implicit none
  integer, parameter :: DP=kind(1.d0)

  integer :: i,k,iflag
  real(DP) :: nu,ss,ds,op
  real(DP) :: xo(3), tridlc(9), burg(3), u(3), t(9)
  real(DP) :: utmp(3), ttmp(9)
  real(DP) :: tri(3,4)
  ! nwarn is now a module variable, no need for common block

!      nu         = poissons ratio
!      xo(3)      = observation point in global NED=123 coordinate system
!      tridlc(9)  = x2,x1,x3 coordinates for the three triangle corners.
!      ss,ds,op   = burg(3)    = burgers vector
!                     (ss = + for LL; ds = + for normal; op = + for expansion)
!      u(3)       = displacements = u1,u2,u3
!      t(9)       = tilts = t11,t12,t13,t21,t22,t23,t31,t32,t33
!                     (where tij = d(ui)/dxj)

  real(DP), parameter :: pi=3.14159265358979323846d0
  real(DP), parameter :: deg2rad=pi/180.d0
  real(DP), parameter :: rad2deg=180.d0/pi

!.....Initialize displacement and tilt arrays.
  u = 0.d0
  t = 0.d0

!.....Fill array with triangle corner coordinates.
!.......(Corner 4 is corner 1 repeated.)
  do i=1,4
    do  k=1,3
      tri(k,i) = tridlc(k+3*mod(i-1,3))
    end do
  end do


!.....Express burger's vector in global 123=NED coordinates.
      call burger2global (ss,ds,op, tri, burg)
!        print '(a,6f8.3)', 'ss,op,ds,burg =(stu) ', ss,op,ds,burg

!.....For each of 3 sides, calculate fields.
  do k = 1,3
    call twoangles (nu,xo,tri(1,k),tri(1,k+1),burg,utmp,ttmp)
    u = u + utmp
    t = t + ttmp
  end do

!.....Apply additional displacement if pt is under the triangle.
  call undertriangle (xo, tri, iflag)

  if (iflag .eq. 1) then
!        Point lies under interior of triangle.
    u = u - burg
  else if (iflag .ne. 0) then
!        Point lies on triangle or on edge of prism below triangle.
    nwarn = nwarn + 1
    if (nwarn.lt.10) then
      print *, ' *** Warning, possible singular location'
      print *, '      for point = ', xo
    else if (nwarn.eq.11) then
      print *, ' *** Warning, possible singular location etc...'
    endif
    u = u - 0.5*burg
  endif

  return
end subroutine dstuart

subroutine burger2global (ss,ds,op, tri, burg)

!.....Express burger's vector in global 123=NED coordinates.
!     Note special definition of ss and ds needed for horizontal dislocations
!         and of ds for vertical dislocations.
  implicit none
  integer, parameter :: DP=kind(1.d0)
  real(DP) :: ss,ds,op
  real(DP) :: tri(3,3), burg(3)
  real(DP) :: side12(3), side13(3)
  real(DP) :: perp(3), hor(3), dwndip(3), b(3)
  real(DP), parameter :: pi=3.14159265358979323846d0
  real(DP), parameter :: deg2rad=pi/180.d0
  real(DP), parameter :: rad2deg=180.d0/pi
  real(DP) :: strikerad,strike,dip

!     tri(3,i) = coordinates of ith corner in 123 coordinate system.
!                   (i=4 returns to i=1)

!.....For a horizontal triangle define strike- and dip-directions by fiat.
  if (tri(3,1).eq.tri(3,2) .and. tri(3,1).eq.tri(3,3)) then
!        Let 1 = ss-axis, 2 = ds-axis, 3 = op-axis
    burg(1) = ss
    burg(2) = ds
    burg(3) = -op
    return
  endif

!.....Otherwise, construct a coordinate system attached to the triangle.

!.....Calculate a normal vector.
  side12(:) = tri(:,2) - tri(:,1)
  side13(:) = tri(:,3) - tri(:,1)
  call cross (side12,side13,perp)
  call unit (perp)

!.....Vertical triangle.
!       By definition, for a vertical triangle, side from which corners
!          appear clockwise is considered top side.

!.....For dipping triangle...
  if (perp(3).lt.0.d0) then
!       Dipping triangle has CCW order to corners viewed from above.
!        Switch corners 2 and 3 and reverse perp.
    call switch3(tri(1,2),tri(1,3))
    perp(:) = -1.d0*perp(:)
!       Perp should now be a unit normal vector pointing down,
!        (away from top), unless triangle is vertical.
  endif

!.....Otherwise define a horizontal vector along strike.
!      (Strike is to left hand as one stands looking down dip
!         -- perp points opposite down-dip.)
  strikerad = atan2(perp(2),perp(1)) + pi*0.5d0
  hor(1) = cos(strikerad)
  hor(2) = sin(strikerad)
  hor(3) = 0.d0

!.....Next a vector in the plane of the triangle pointing in the
!       down-dip direction.
  call cross (perp,hor,dwndip)

!.....Express Burgers vector in this coordinate system.
  burg(:) = ss*hor(:) + ds*dwndip(:) - op*perp(:)

!.....Do a check.
  strike = strikerad * rad2deg
  dip = acos(-perp(3)) * rad2deg
  call global2triburger (burg, strike, dip, b)

  return
end subroutine burger2global

subroutine global2triburger (burg, strike, dip, b)

!.....Express burger's vector as b="ss,op,ds" (i.e., in local coordinates
!      attached to dlc.)  Order is here set to match Comninou & Dunders
!      conventions.
  implicit none
  integer, parameter :: DP=kind(1.d0)
  real(DP) :: burg(3), b(3)
  real(DP) :: perp(3), hor(3), dwndip(3)
  real(DP), parameter :: pi=3.14159265358979323846d0
  real(DP), parameter :: deg2rad=pi/180.d0
  real(DP), parameter :: rad2deg=180.d0/pi
  real(DP) :: strike,dip,strikerad,diprad,dipdirrad

!.....Construct a coordinates system attached to the triangle.

!.....For a horizontal triangle define strike- and dip-directions by fiat.
  if (dip .eq. 0.d0) then
!        Let 1 = ss-axis, 2 = op-axis, 3 = ds-axis
    strike = 0.d0
  endif

!.....Define a horizontal vector along strike.
  strikerad = strike * deg2rad
  hor(1) = cos(strikerad)
  hor(2) = sin(strikerad)
  hor(3) = 0.d0

!.....Next a vector in the plane of the triangle pointing in the
!       down-dip dipdirection.
  diprad = dip * deg2rad
  dipdirrad = strikerad + pi*0.5
  dwndip(1) = cos (diprad) * cos(dipdirrad)
  dwndip(2) = cos (diprad) * sin(dipdirrad)
  dwndip(3) = sin (diprad)

!.....Finally a perpendicular vector
  call cross (hor,dwndip,perp)

!.....Get local burgers vector.  
  b(1) =   dot_product(burg,hor)
  b(2) = - dot_product(burg,perp)
  b(3) =   dot_product(burg,dwndip)

  return
end subroutine global2triburger

subroutine twoangles (nu,xo,xcorna,xcornb, burg, u,t)

!     Calculates the effect of a vertical side formed by two angles
!       located at points xcorna and xcornb.
  implicit none
 integer, parameter :: DP=kind(1.d0)
  real(DP) :: nu, xo(3), xcorna(3), xcornb(3), burg(3), u(3), t(9)
  real(DP) :: xop(3), xca(3), xcb(3), b(3)
  real(DP) :: v(3,3), dv(3,3,3)
  real(DP) :: va(3,3), dva(3,3,3)
  real(DP) :: vb(3,3), dvb(3,3,3)

!     nu     = Poissons ratio
!     xo(3)  = observation point in x1,x2,x3 = NED coordinates.
!     xca(3) = x1,x2,x3 coordinates of corner A.
!     xcb(3) = x1,x2,x3 coordinates of corner B.
  real(DP), parameter :: pi=3.14159265358979323846d0
  real(DP), parameter :: deg2rad=pi/180.d0
  real(DP), parameter :: rad2deg=180.d0/pi
  integer :: iswitch,j,k,iret
  real(DP) :: theta,strike,dip,rot,deptha,depthb,xlen,betarad,x1,x2,x3

!.....Copy corners.
  xca(:) = xcorna(:)
  xcb(:) = xcornb(:)

!.....See which is shallower, A or B, and switch if necessary.
  if (xca(3).gt.xcb(3)) then
!       B is shallower than A
    call switch3 (xca, xcb)
    iswitch = 1
  else
    iswitch = 0
  endif

!.....Calculate angle from A to B.
  theta = atan2(xcb(2)-xca(2), xcb(1)-xca(1)) * rad2deg

!.....Set parameters.
  strike = theta
  deptha = xca(3)
  depthb = xcb(3)
  xlen = sqrt((xcb(1)-xca(1))**2 + (xcb(2)-xca(2))**2)
  betarad =  atan2 (xlen, xcb(3)-xca(3))

!.....Locate the observation point relative to the shallower of A or B.
  xop(1) =  xo(1) - xca(1)
  xop(2) =  xo(2) - xca(2)
  xop(3) =  xo(3)
!     If strike .ne. 0, need to rotate obs points and output fields.
  rot = -strike
  call vector_rot3(xop,-rot)
  x1 = xop(1)
  x2 = xop(2)
  x3 = xop(3)

!.....Compute the displacement and tilts.
  call comdun2 ( nu, x1,x2,x3, deptha, betarad, va, dva, iret )
  call comdun2 ( nu, x1-xlen,x2,x3, depthb, betarad, vb, dvb, iret )

!.....Combine the results from both angles.
  v(:,:) = va(:,:) - vb(:,:)
  dv(:,:,:) = dva(:,:,:) - dvb(:,:,:)

!.....Calculate displacements(u) and tilts(t)
  dip = 90.d0
  call global2triburger (burg, strike, dip, b)

  if (iswitch.eq.1)  b(:) = -1d0*b(:)

  u = matmul(b,v)
  do j=1,3
    do k=1,3
      t(k + 3*(j-1)) = dot_product(b(:),dv(:,j,k))
    end do
  end do

!.....Rotate fields back to global NED=123 coordinate system.
  call vector_rot3 (u,rot)
  call tensor_rot3 (t,rot)

  return
end subroutine twoangles

subroutine comdun2 ( nu, x1, x2, x3, a, beta, v, dv, iret )
!     Performs some checks before calling comdun
use comdun_module
  implicit none
 integer, parameter :: DP=kind(1.d0)
  real(DP) :: nu,x1,x2,x3,a,beta
  integer :: iret,nwarn
  real(DP) :: v(3,3), dv(3,3,3)
  real(DP) :: vp(3,3), dvp(3,3,3)
  real(DP) :: vm(3,3), dvm(3,3,3)
  real(DP) :: betatol = 1.d-4
  real(DP) :: tol = 1.d-4
  real(DP) :: x2p,x2m,tan2angle

  ! nwarn is now a module variable

!.....Check for angle beta - if close to zero, return zero result.
  if (beta .lt. betatol) then
    v(:,:) = 0.d0
    dv(:,:,:) = 0.d0
    return
  endif

!.....Check for distance from x1-x3 plane.  If too close, average.
  if (abs(x2).lt. tol) then
    x2p = abs(x2) + tol
    x2m = -x2p
    call comdun ( nu, x1, x2p, x3, a, beta, vp, dvp  )
    call comdun ( nu, x1, x2m, x3, a, beta, vm, dvm )
    v(:,:) = 0.5d0*(vp(:,:)+vm(:,:))
    dv(:,:,:) = 0.5d0*(dvp(:,:,:)+dvm(:,:,:))
  else
!        Normal call to comdun
    call comdun ( nu, x1, x2, x3, a, beta, v, dv )
  endif

!c.....Check to see if point is inside angle.
  tan2angle = atan2(x1,x3-a)
  if (abs(x2).lt.tol .and. (tan2angle.le.beta.and.tan2angle.ge.0.d0)) then
!        Close to the displacement singularity.
    nwarn = nwarn + 1
    if (nwarn.lt.10) then
      print *, '  *** Displacement singularity...'
    else if (nwarn.eq.11) then
      print *, '  *** Displacement singularity...etc,etc'
    endif
    iret = 1
  endif

  return
end subroutine comdun2

function inside (x0,y0, px,py,n)

!     From J.Berger, et al. BSSA v. 74, p. 1849-1862, October 1984.

!     x0,y0 = point to test
!     px,py,n = corners of polygon of n sides.
!     Return value = 0 if point is outside
!                  = +/-1 if point is inside
!                  = 2 if point is on an edge or vertex.
  implicit none
  integer, parameter :: DP=kind(1.d0)
  integer :: inside,i, isicr
  integer, intent(in) :: n
  real(DP), intent(in) :: x0, y0
  real(DP), intent(in) :: px(n), py(n)

  inside = 0
  do i=1,n-1
    isicr = ksicr(px(i)-x0,py(i)-y0,px(i+1)-x0,py(i+1)-y0)
    if (isicr.eq.4) then
      inside = 2
      return
    endif
    inside = inside + isicr
  end do

  isicr = ksicr(px(n)-x0,py(n)-y0,px(1)-x0,py(1)-y0)

  if(isicr.eq.4) then
    inside = 2
    return
  endif

  inside = (inside + isicr)/2
  return
end function inside

subroutine undertriangle (xo, tri, iflag)

!     Returns a flag = 1 if the point xo lies under the triangle interior.
!     Returns a flag = 2 if the point xo lies under a triangle edge or corner.
!     Returns a flag = 3 if the point xo lies on a triangle edge or corner.
!     Returns a flag = 0 otherwise.
  implicit none
 integer, parameter :: DP=kind(1.d0)
  integer :: iflag
  real(DP) :: xo(3), tri(3,4)
  real(DP) :: xc(3), yc(3)
  real(DP) :: side12(3), side13(3), perp(3)
  real(DP) :: tol = 1.d-6
  real(DP) :: xpt,ypt,zpt,zontri,ht
 !  integer :: inside

  iflag=0
!     Make arrays of x,y coords of triangle corners.
  xc(1:3) = tri(1,1:3)
  yc(1:3) = tri(2,1:3)

  xpt = xo(1)
  ypt = xo(2)
  zpt = xo(3)
  iflag = inside (xpt,ypt, xc,yc,3)

!     Point is outside vertical prism....
  if (iflag.eq.0)  return

!     Point is inside vertical prism.... is it below triangle?

!     Form perpendicular to the triangle.
  side12(:) = tri(:,2) - tri(:,1)
  side13(:) = tri(:,3) - tri(:,1)
  call cross (side12,side13,perp)
  call unit (perp)

  if (perp(3).eq.0.0) then
!        Triangle is vertical, return edge flag.
    iflag = 2
    return
  endif

!     Calculate height (+) below or (-) above triangle surface.
  zontri = (dot_product(tri(:,1),perp(:)) - xpt*perp(1) - ypt*perp(2))/perp(3)
  ht = zpt - zontri
  if (abs(ht).lt.tol) then
    iflag = 3
  else if (ht.lt.-tol) then
    iflag = 0
  else
    iflag = abs(iflag)
  endif

  return
end subroutine undertriangle



function ksicr (x1,y1,x2,y2)
  implicit none
 integer, parameter :: DP=kind(1.d0)
  integer :: ksicr
  real(DP), intent(in) :: x1,x2,y1,y2
  
  ! Convert goto-based logic to proper Fortran 90
  if (y1*y2 .gt. 0.0d0) then
    ksicr = 0
    return
  endif

  if (x1*y2 .ne. x2*y1 .or. x1*x2 .gt. 0.0d0) then
    if (y1*y2 .lt. 0.0d0) then
      if (y1.gt.0.0d0) then
        if (x1*y2 .ge. y1*x2) then
          ksicr = 0
        else
          ksicr = 2
        endif
      else
        if (x1*y2 .ge. y1*x2) then
          ksicr = 0
        else
          ksicr = 2
        endif
      endif
    else
      if (y2.eq.0.0d0) then
        if (y1.eq.0.0d0 .or. x2.gt.0.0d0) then
          ksicr = 0
        else
          if (y1.gt.0.0d0) then
            ksicr = -1
          else
            ksicr = 1
          endif
        endif
      else
        if (x1.gt.0.0d0) then
          ksicr = 0
        else
          if (y2.gt.0.0d0) then
            ksicr = 1
          else
            ksicr = -1
          endif
        endif
      endif
    endif
  else
    ksicr = 4
  endif
  
  return
end function ksicr

subroutine vector_rot3 (u,phi)

!     Modified from Laurie Erickson's version of the DIS3D program.

!     Rotates vector u around the x3 axis in a CW sense
!        by angle phi in the x1-x2 plane.
  implicit none
 integer, parameter :: DP=kind(1.d0)
  real(DP) :: phi
  real(DP) :: u(3), utmp(3)
  real(DP), parameter :: pi=3.14159265358979323846d0
  real(DP), parameter :: deg2rad=pi/180.d0
  real(DP), parameter :: rad2deg=180.d0/pi

  real(DP) :: phirad,cp,sp

  if (phi.eq.0.0) return

  phirad = phi*deg2rad
  cp = cos(phirad)
  sp = sin(phirad)

  utmp(1:2) = u(1:2)

  u(1) =  cp*utmp(1) + sp*utmp(2)
  u(2) = -sp*utmp(1) + cp*utmp(2)

  return
end subroutine vector_rot3

subroutine tensor_rot3 (t,phi)

!     Modified from Laurie Erickson's version of the DIS3D program.

!     Rotates a general tensor t by the angle phi around the
!       x3 axis (within the x1-x2 plane).
  implicit none
integer, parameter :: DP=kind(1.d0)
  
  real(DP), intent(inout) :: phi
  real(DP), intent(inout) :: t(9)
  real(DP) :: ttmp(9)
  real(DP), parameter :: pi=3.14159265358979323846d0
  real(DP), parameter :: deg2rad=pi/180.d0
  real(DP), parameter :: rad2deg=180.d0/pi
  real(DP) :: phirad,cp,sp

  if (phi.eq.0.d0) return

  phirad = phi*deg2rad
  cp = cos(phirad)
  sp = sin(phirad)

  ttmp(1:8) = t(1:8)

  t(1) =  (cp**2)*ttmp(1)+cp*sp*ttmp(2)+cp*sp*ttmp(4)+(sp**2)*ttmp(5)
  t(2) = -cp*sp*ttmp(1)+(cp**2)*ttmp(2)-(sp**2)*ttmp(4)+sp*cp*ttmp(5)
  t(3) =  cp*ttmp(3) + sp*ttmp(6)
  t(4) = -sp*cp*ttmp(1)-(sp**2)*ttmp(2)+(cp**2)*ttmp(4)+cp*sp*ttmp(5)
  t(5) =  (sp**2)*ttmp(1)-sp*cp*ttmp(2)-sp*cp*ttmp(4)+(cp**2)*ttmp(5)
  t(6) = -sp*ttmp(3) + cp*ttmp(6)
  t(7) =  cp*ttmp(7) + sp*ttmp(DP)
  t(DP) = -sp*ttmp(7) + cp*ttmp(DP)

  return
end subroutine tensor_rot3

subroutine switch (a, b)
!     Switches real numbers.
  implicit none
integer, parameter :: DP=kind(1.d0)
  real(DP) :: a,b,temp
  temp = a
  a = b
  b = temp

  return
end subroutine switch

subroutine switch3 (xca, xcb)
!     Interchanges two vectors
  implicit none
integer, parameter :: DP=kind(1.d0)
  real(DP) :: xca(3), xcb(3), xtmp(3)

  xtmp(:) = xca(:)
  xca(:) = xcb(:)
  xcb(:) = xtmp(:)

  return
end subroutine switch3


subroutine cross (a,b,acrsb)

!     Calculates the cross product, acrsb, of vectors a and b.
  implicit none
 integer, parameter :: DP=kind(1.d0)
  real(DP) :: a(3), b(3), acrsb(3)

  acrsb(1)= a(2)*b(3) - a(3)*b(2)
  acrsb(2)= a(3)*b(1) - a(1)*b(3)
  acrsb(3)= a(1)*b(2) - a(2)*b(1)

  return
end subroutine cross

subroutine unit(x)
!     Makes a unit vector out of x.
  implicit none
 integer, parameter :: DP=kind(1.d0)
  real(DP) :: x(3)
  real(DP) :: d
  d = sqrt( dot_product(x,x))

 if(d.gt.0.d0) then
   x(:) = x(:)/d
  else
   print *, ' ** Warning... zero unit vector in sub unit.'
   stop
  endif

  return
end subroutine unit

end module mod_dtrigreen
