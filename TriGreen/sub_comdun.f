
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Calls Stuart's angle subroutine to get triangle dislocation results.

c     Comninou and Dundurs use a 1,2,3 = NED coordinate system.
c       Their simple angular dislocation is in the 1,3 plane
c       with its tip under the origin and its arm pointing toward +x1.

c---+-|--1----+----2----+----3----+----4----+----5----+----6----+----7--

      subroutine comdun ( nu, x1, x2, x3, a, beta, v, dv )

c     Programmed by Bill Stuart (1997) from
c       Comninou and Dundurs, 1975, "The angular dislocation in a half space",
c       Journal of Elasticity, v.5, n.3,4, p. 203-216.

c       nu             - Poisson ratio
c       x1,x2,x3       - obs pt: x1 +rt,  x2 out,  x3 dwn
c       a              - depth, >0, to ang disloc vertex
c       beta           - acute angle, rad, of disloc
c       b1,b2,b3       - burg. vecs. in x(i)
c       v(b,comp)      - obs. displ
c       dv(b,comp,dx)  - displ deriv

      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 nu, x1, x2, x3, a, beta, v, dv

      dimension v(3,3),vc(3,3), dv(3,3,3),dvc(3,3,3)

      cb  = cos(beta)
      sb  = sin(beta)
      ctb = cb/sb

      y1  = x1
      y2  = x2
      y3  = x3 - a
      yb3 = y3 + 2.*a

      z1  = y1*cb - y3*sb
      z2  = y2
      z3  = y1*sb + y3*cb

      zb1 =  y1*cb + yb3*sb
      zb3 = -y1*sb + yb3*cb

      r   = sqrt(y1**2 + y2**2 + y3**2)
      rb  = sqrt(y1*y1 + y2*y2 + yb3*yb3 )

      f   = -atan2 (y2,y1) + atan2 (y2,z1) +
     .      atan2 (y2*r*sb , (y1*z1 + y2**2*cb) )

      fb  = -atan2 (y2,y1) + atan2 (y2,zb1) +
     .       atan2 (y2*rb*sb , (y1*zb1 + y2**2*cb) )


c     derivs f, fb

      dfdy1 =
     .  y2*(1/(y1**2 + y2**2) - cb/(y2**2 + z1**2) +
     -    (sb*(cb*y1*(-r**2 + y2**2) + (-r**2 + y1**2)*z1))/
     -     (r*(r**2*sb**2*y2**2 + (cb*y2**2 + y1*z1)**2)))


      dfdy2 =
     .  -(y1/(y1**2 + y2**2)) + z1/(y2**2 + z1**2) +
     -  (sb*(cb*(-(r**2*y2**2) + y2**4) + y1*(r**2 + y2**2)*z1))/
     -   (r*(r**2*sb**2*y2**2 + (cb*y2**2 + y1*z1)**2))


      dfdy3 =
     .  sb*y2*(1/(y2**2 + z1**2) +
     -    (r**2*sb*y1 + cb*y2**2*y3 + y1*y3*z1)/
     -     (r*(r**2*sb**2*y2**2 + (cb*y2**2 + y1*z1)**2)))


      dfbdy1 =
     .  y2*(1/(y1**2 + y2**2) - cb/(y2**2 + zb1**2) +
     -    (sb*(cb*y1*(-rb**2 + y2**2) + (-rb**2 + y1**2)*zb1))/
     -     (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2)))


      dfbdy2 =
     .  -(y1/(y1**2 + y2**2)) + zb1/(y2**2 + zb1**2) +
     -  (sb*(cb*(-(rb**2*y2**2) + y2**4) + y1*(rb**2 + y2**2)*zb1))/
     -   (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2))


      dfbdy3 =
     .  sb*y2*(-(1/(y2**2 + zb1**2)) +
     -    (-(rb**2*sb*y1) + cb*y2**2*yb3 + y1*yb3*zb1)/
     -     (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2)))


c     v  displacements

      v(1,1) =
     .  -(y1*y2*(1/(r*(r - y3)) + 1/(rb*(rb + yb3)))) -
     -  cb*y2*((r*sb - y1)/(r*(r - z3)) +
     -     (rb*sb - y1)/(rb*(rb + zb3))) +
     -  2*(1 - nu)*(-2*atan2(y2,y1) + atan2(y2,z1) +
     -     atan2((r*sb*y2),(cb*y2**2 + y1*z1)) + atan2(y2,zb1) +
     -     atan2((rb*sb*y2),(cb*y2**2 + y1*zb1)))


      v(1,2) =
     .  -(y2**2*(1/(r*(r - y3)) + 1/(rb*(rb + yb3)) -
     -       cb*(1/(r*(r - z3)) + 1/(rb*(rb + zb3))))) +
     -  (1 - 2*nu)*(log(r - y3) + log(rb + yb3) -
     -     cb*(log(r - z3) + log(rb + zb3)))


      v(1,3) =
     .  y2*(1/r - 1/rb - cb*((cb*r - y3)/(r*(r - z3)) -
     -       (cb*rb + yb3)/(rb*(rb + zb3))))


      v(2,1) =
     .  y1**2*(1/(r*(r - y3)) + 1/(rb*(rb + yb3))) +
     -  ((r*sb - y1)*z1)/(r*(r - z3)) +
     -  ((rb*sb - y1)*zb1)/(rb*(rb + zb3)) +
     -  (-1 + 2*nu)*(log(r - y3) + log(rb + yb3) -
     -     cb*(log(r - z3) + log(rb + zb3)))


      v(2,2) =
     .  y1*y2*(1/(r*(r - y3)) + 1/(rb*(rb + yb3))) -
     -  y2*(z1/(r*(r - z3)) + zb1/(rb*(rb + zb3))) +
     -  2*(1 - nu)*(-2*atan2(y2,y1) + atan2(y2,z1) +
     -     atan2((r*sb*y2),(cb*y2**2 + y1*z1)) + atan2(y2,zb1) +
     -     atan2((rb*sb*y2),(cb*y2**2 + y1*zb1)))




      v(2,3) =
     .  -((1/r - 1/rb)*y1) + ((cb*r - y3)*z1)/(r*(r - z3)) -
     -  ((cb*rb + yb3)*zb1)/(rb*(rb + zb3)) +
     -  (-1 + 2*nu)*sb*(log(r - z3) - log(rb + zb3))


      v(3,1) =
     .sb*y2*((r*sb - y1)/(r*(r - z3)) +
     . (rb*sb - y1)/(rb*(rb + zb3)))


      v(3,2) =
     .  -(sb*y2**2*(1/(r*(r - z3)) + 1/(rb*(rb + zb3)))) +
     -  (1 - 2*nu)*sb*(log(r - z3) + log(rb + zb3))


      v(3,3) =
     .  sb*y2*((cb*r - y3)/(r*(r - z3)) -
     -     (cb*rb + yb3)/(rb*(rb + zb3))) +
     -  2*(1 - nu)*(atan2(y2,z1) +
     -     atan2((r*sb*y2),(cb*y2**2 + y1*z1)) - atan2(y2,zb1) -
     -     atan2((rb*sb*y2),(cb*y2**2 + y1*zb1)))


c     vc  displacements

      vc(1,1) =
     .  (a*ctb*y2*(-a + yb3))/rb**3 +
     -  ((1 - 2*nu)*y2*(ctb*(1 - 2*nu - a/rb) -
     -       ((nu + a/rb)*y1)/(rb + yb3)))/(rb + yb3) +
     -  (y2*(-a + yb3)*(-(ctb*(1 - 2*nu)) + (a*y1)/rb**2 +
     -       ((2*nu + a/rb)*y1)/(rb + yb3)))/(rb*(rb + yb3)) +
     -  (cb*ctb*(1 - 2*nu)*(cb + a/rb)*y2)/(rb + zb3) +
     -  (y2*(-a + yb3)*(-((a*cb*ctb*yb3)/rb**2) +
     -       (cb*(2*cb*(1 - nu)*(rb*sb - y1) +
     -            ctb*(cb*(1 - 2*nu) - a/rb)*(cb*rb + yb3)))/(rb + zb3))
     -     )/(rb*(rb + zb3)) -
     -  2*ctb**2*(1 - 2*nu)*(1 - nu)*
     -   (-atan2(y2,y1) + atan2(y2,zb1) +
     -     atan2((rb*sb*y2),(cb*y2**2 + y1*zb1)))


      vc(1,2) =
     .  -((a*ctb*y1*(-a + yb3))/rb**3) -
     -  ((1 - 2*nu)*(-a + ctb*(1 - 2*nu - a/rb)*y1 + nu*yb3 +
     -       ((nu + a/rb)*y2**2)/(rb + yb3)))/(rb + yb3) +
     -  ((-a + yb3)*(-2*nu + (-a + ctb*(1 - 2*nu)*y1)/rb +
     -       (a*y2**2)/rb**3 + ((2*nu + a/rb)*y2**2)/(rb*(rb + yb3))))/
     -   (rb + yb3) - (ctb*(1 - 2*nu)*(cb + a/rb)*zb1)/(rb + zb3) +
     -  ((-a + yb3)*(cb**2 + (a*ctb*yb3*zb1)/rb**3 -
     -       (a*cb + ctb*(1 - 2*nu)*zb1)/rb -
     -       (cb**2*y2**2 - (a*ctb*(cb*rb + yb3)*zb1)/rb)/
     -        (rb*(rb + zb3))))/(rb + zb3) +
     -  (1 - 2*nu)*((2*ctb**2*(1 - nu) - nu)*log(rb + yb3) -
     -     cb*(1 + 2*ctb**2*(1 - nu) - 2*nu)*log(rb + zb3))


      vc(1,3) =
     .  (y2*(-a + yb3)*(a/rb**2 + (2*nu)/(rb + yb3)))/rb +
     -  (cb*y2*(-a + yb3)*(1 - 2*nu - (a*yb3)/rb**2 -
     -       ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/(rb*(rb + zb3)) +
     -  2*(1 - nu)*(((2*nu + a/rb)*y2)/(rb + yb3) -
     -     (cb*(cb + a/rb)*y2)/(rb + zb3) +
     -     ctb*(1 - 2*nu)*(-atan2(y2,y1) + atan2(y2,zb1) +
     -        atan2((rb*sb*y2),(cb*y2**2 + y1*zb1))))


      vc(2,1) =
     .  -((a*ctb*y1*(-a + yb3))/rb**3) +
     -  ((1 - 2*nu)*(-a + ctb*(-1 + 2*nu)*y1 + (a*ctb*y1)/rb + nu*yb3 +
     -       ((nu + a/rb)*y1**2)/(rb + yb3)))/(rb + yb3) +
     -  ((-a + yb3)*(2*nu - (a*y1**2)/rb**3 +
     -       (a + ctb*(1 - 2*nu)*y1)/rb -
     -       ((2*nu + a/rb)*y1**2)/(rb*(rb + yb3))))/(rb + yb3) -
     -  (ctb*(1 - 2*nu)*(-((a*(rb*sb - y1))/(cb*rb)) + cb*zb1))/
     -   (rb + zb3) + (ctb*(-a + yb3)*
     -     (-(cb*sb) + (a*y1*yb3)/(cb*rb**3) +
     -       ((rb*sb - y1)*(2*cb*(1 - nu) -
     -            ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/rb))/
     -   (rb + zb3) + (1 - 2*nu)*
     -   ((2*ctb**2*(1 - nu) + nu)*log(rb + yb3) -
     -     cb*(1 + 2*ctb**2*(1 - nu))*log(rb + zb3))


      vc(2,2) =
     .  -((a*ctb*y2*(-a + yb3))/rb**3) +
     -  ((1 - 2*nu)*y2*(ctb*(-1 + 2*nu + a/rb) +
     -       ((nu + a/rb)*y1)/(rb + yb3)))/(rb + yb3) +
     -  (y2*(-a + yb3)*(ctb*(1 - 2*nu) - (2*nu*y1)/(rb + yb3) -
     -       (a*y1*(1/rb + 1/(rb + yb3)))/rb))/(rb*(rb + yb3)) -
     -  (ctb*(1 - 2*nu)*(1 + a/(cb*rb))*y2)/(rb + zb3) +
     -  (ctb*y2*(-a + yb3)*(-2*cb*(1 - nu) + (a*yb3)/(cb*rb**2) +
     -       ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb*(rb + zb3)) + 2*ctb**2*(1 - 2*nu)*(1 - nu)*
     -   (-atan2(y2,y1) + atan2(y2,zb1) +
     -     atan2((rb*sb*y2),(cb*y2**2 + y1*zb1)))


      vc(2,3) =
     .  (-2*(1 - nu)*(2*nu + a/rb)*y1)/(rb + yb3) +
     -  ((-a + yb3)*(ctb*(1 - 2*nu) - (a*y1)/rb**2 -
     -       (2*nu*y1)/(rb + yb3)))/rb +
     -  (2*(1 - nu)*(cb + a/rb)*zb1)/(rb + zb3) -
     -  ((-a + yb3)*(cb*sb +
     -       (ctb*(cb*rb + yb3)*
     -          (2*cb*(1 - nu) - (cb*rb + yb3)/(rb + zb3)))/rb +
     -   (a*(sb - (yb3*zb1)/rb**2 -
     -   ((cb*rb + yb3)*zb1)/(rb*(rb + zb3))))/rb))/(rb + zb3)
     -   - 2*ctb*(1 - 2*nu)*(1 - nu)*(log(rb + yb3)
     .   - cb*log(rb + zb3))


      vc(3,1) =
     .  -((y2*(-a + yb3)*(a/rb**2 + 1/(rb + yb3)))/rb) +
     -  (1 - 2*nu)*(((1 + a/rb)*y2)/(rb + yb3) -
     -     (cb*(cb + a/rb)*y2)/(rb + zb3)) +
     -  (cb*y2*(-a + yb3)*((a*yb3)/rb**2 +
     -       ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/(rb*(rb + zb3))


      vc(3,2) =
     .  (y1*(-a + yb3)*(a/rb**2 + 1/(rb + yb3)))/rb -
     -  ((-a + yb3)*((cb - a/rb)*sb + ((1 + (a*yb3)/rb**2)*zb1)/rb -
     -       (cb*sb*y2**2 - (a*(cb*rb + yb3)*zb1)/rb)/(rb*(rb + zb3))))/
     -   (rb + zb3) + (1 - 2*nu)*
     -   (-(((1 + a/rb)*y1)/(rb + yb3)) +
     -     ((cb + a/rb)*zb1)/(rb + zb3) - sb*log(rb + zb3))


      vc(3,3) =
     .  (sb*y2*(-a + yb3)*(1 + (a*yb3)/rb**2 +
     -       ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/(rb*(rb + zb3)) +
     -  2*(1 - nu)*(((cb + a/rb)*sb*y2)/(rb + zb3) - atan2(y2,y1) +
     -     atan2(y2,zb1) + atan2((rb*sb*y2),(cb*y2**2 + y1*zb1)))


c     v   derivatives

      dv(1,1,1) =
     .  -(y2*(1/(r*(r - y3)) + 1/(rb*(rb + yb3)))) -
     -  y1*y2*(-(y1/(r**2*(r - y3)**2)) - y1/(r**3*(r - y3)) -
     -     y1/(rb**2*(rb + yb3)**2) - y1/(rb**3*(rb + yb3))) +
     -  2*(1 - nu)*(y2*(1/(y1**2 + y2**2) - cb/(y2**2 + z1**2) +
     -        (sb*(cb*y1*(-r**2 + y2**2) + (-r**2 + y1**2)*z1))/
     -         (r*(r**2*sb**2*y2**2 + (cb*y2**2 + y1*z1)**2))) +
     -     y2*(1/(y1**2 + y2**2) - cb/(y2**2 + zb1**2) +
     -        (sb*(cb*y1*(-rb**2 + y2**2) + (-rb**2 + y1**2)*zb1))/
     -         (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2)))) -
     -  cb*y2*(-(((r*sb - y1)*(-sb + y1/r))/(r*(r - z3)**2)) -
     -     ((r*sb - y1)*y1)/(r**3*(r - z3)) +
     -     (-1 + (sb*y1)/r)/(r*(r - z3)) -
     -     ((rb*sb - y1)*(-sb + y1/rb))/(rb*(rb + zb3)**2) -
     -     ((rb*sb - y1)*y1)/(rb**3*(rb + zb3)) +
     -     (-1 + (sb*y1)/rb)/(rb*(rb + zb3)))


      dv(1,1,2) =
     .  -(y1*(1/(r*(r - y3)) + 1/(rb*(rb + yb3)))) -
     -  y1*y2*(-(y2/(r**2*(r - y3)**2)) - y2/(r**3*(r - y3)) -
     -     y2/(rb**2*(rb + yb3)**2) - y2/(rb**3*(rb + yb3))) +
     -  2*(1 - nu)*((-2*y1)/(y1**2 + y2**2) + z1/(y2**2 + z1**2) +
     -     (sb*(cb*(-(r**2*y2**2) + y2**4) + y1*(r**2 + y2**2)*z1))/
     -      (r*(r**2*sb**2*y2**2 + (cb*y2**2 + y1*z1)**2)) +
     -     zb1/(y2**2 + zb1**2) +
     -     (sb*(cb*(-(rb**2*y2**2) + y2**4) + y1*(rb**2 + y2**2)*zb1))/
     -      (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2))) -
     -  cb*((r*sb - y1)/(r*(r - z3)) + (rb*sb - y1)/(rb*(rb + zb3))) -
     -  cb*y2*(-(((r*sb - y1)*y2)/(r**2*(r - z3)**2)) +
     -     (sb*y2)/(r**2*(r - z3)) - ((r*sb - y1)*y2)/(r**3*(r - z3)) -
     -     ((rb*sb - y1)*y2)/(rb**2*(rb + zb3)**2) +
     -     (sb*y2)/(rb**2*(rb + zb3)) -
     -     ((rb*sb - y1)*y2)/(rb**3*(rb + zb3)))


      dv(1,1,3) =
     . -(y1*y2*(-(y3/(r**3*(r - y3))) - (-1 + y3/r)/(r*(r - y3)**2) -
     - yb3/(rb**3*(rb + yb3)) - (1 + yb3/rb)/(rb*(rb + yb3)**2)))
     -   + 2*(1 - nu)*(sb*y2*
     -      (1/(y2**2 + z1**2) +
     -        (r**2*sb*y1 + cb*y2**2*y3 + y1*y3*z1)/
     -         (r*(r**2*sb**2*y2**2 + (cb*y2**2 + y1*z1)**2))) +
     -     sb*y2*(-(1/(y2**2 + zb1**2)) +
     -        (-(rb**2*sb*y1) + cb*y2**2*yb3 + y1*yb3*zb1)/
     -         (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2)))) -
     -  cb*y2*(-(((r*sb - y1)*(-cb + y3/r))/(r*(r - z3)**2)) +
     -  (sb*y3)/(r**2*(r - z3)) - ((r*sb - y1)*y3)/(r**3*(r - z3)) -
     -     ((rb*sb - y1)*(cb + yb3/rb))/(rb*(rb + zb3)**2) +
     -     (sb*yb3)/(rb**2*(rb + zb3)) -
     -     ((rb*sb - y1)*yb3)/(rb**3*(rb + zb3)))


      dv(1,2,1) =
     .  -(y2**2*(-(y1/(r**2*(r - y3)**2)) - y1/(r**3*(r - y3)) -
     -       y1/(rb**2*(rb + yb3)**2) - y1/(rb**3*(rb + yb3)) -
     -       cb*(-((-sb + y1/r)/(r*(r - z3)**2)) - y1/(r**3*(r - z3)) -
     -          (-sb + y1/rb)/(rb*(rb + zb3)**2) - y1/(rb**3*(rb + zb3))
     -          ))) + (1 - 2*nu)*
     -   (y1/(r*(r - y3)) + y1/(rb*(rb + yb3)) -
     -     cb*((-sb + y1/r)/(r - z3) + (-sb + y1/rb)/(rb + zb3)))


      dv(1,2,2) =
     .  -2*y2*(1/(r*(r - y3)) + 1/(rb*(rb + yb3)) -
     -     cb*(1/(r*(r - z3)) + 1/(rb*(rb + zb3)))) -
     -  y2**2*(-(y2/(r**2*(r - y3)**2)) - y2/(r**3*(r - y3)) -
     -     y2/(rb**2*(rb + yb3)**2) - y2/(rb**3*(rb + yb3)) -
     -     cb*(-(y2/(r**2*(r - z3)**2)) - y2/(r**3*(r - z3)) -
     -        y2/(rb**2*(rb + zb3)**2) - y2/(rb**3*(rb + zb3)))) +
     -  (1 - 2*nu)*(y2/(r*(r - y3)) + y2/(rb*(rb + yb3)) -
     -     cb*(y2/(r*(r - z3)) + y2/(rb*(rb + zb3))))


      dv(1,2,3) =
     .  -(y2**2*(-(y3/(r**3*(r - y3))) - (-1 + y3/r)/(r*(r - y3)**2) -
     -       yb3/(rb**3*(rb + yb3)) - (1 + yb3/rb)/(rb*(rb + yb3)**2) -
     -       cb*(-((-cb + y3/r)/(r*(r - z3)**2)) - y3/(r**3*(r - z3)) -
     -          (cb + yb3/rb)/(rb*(rb + zb3)**2) -
     -          yb3/(rb**3*(rb + zb3))))) +
     -  (1 - 2*nu)*((-1 + y3/r)/(r - y3) + (1 + yb3/rb)/(rb + yb3) -
     -     cb*((-cb + y3/r)/(r - z3) + (cb + yb3/rb)/(rb + zb3)))


      dv(1,3,1) =
     .  y2*(-(y1/r**3) + y1/rb**3 -
     -    cb*(-(((-sb + y1/r)*(cb*r - y3))/(r*(r - z3)**2)) +
     -       (cb*y1)/(r**2*(r - z3)) -
     -       (y1*(cb*r - y3))/(r**3*(r - z3)) +
     -       ((-sb + y1/rb)*(cb*rb + yb3))/(rb*(rb + zb3)**2) -
     -       (cb*y1)/(rb**2*(rb + zb3)) +
     -       (y1*(cb*rb + yb3))/(rb**3*(rb + zb3))))


      dv(1,3,2) =
     .  1/r - 1/rb - cb*((cb*r - y3)/(r*(r - z3)) -
     -     (cb*rb + yb3)/(rb*(rb + zb3))) +
     -  y2*(-(y2/r**3) + y2/rb**3 -
     -     cb*(-((y2*(cb*r - y3))/(r**2*(r - z3)**2)) +
     -        (cb*y2)/(r**2*(r - z3)) -
     -        (y2*(cb*r - y3))/(r**3*(r - z3)) +
     -        (y2*(cb*rb + yb3))/(rb**2*(rb + zb3)**2) -
     -        (cb*y2)/(rb**2*(rb + zb3)) +
     -        (y2*(cb*rb + yb3))/(rb**3*(rb + zb3))))


      dv(1,3,3) =
     .  y2*(-(y3/r**3) + yb3/rb**3 -
     -    cb*(-(((cb*r - y3)*(-cb + y3/r))/(r*(r - z3)**2)) -
     -       ((cb*r - y3)*y3)/(r**3*(r - z3)) +
     -       (-1 + (cb*y3)/r)/(r*(r - z3)) +
     -       ((cb*rb + yb3)*(cb + yb3/rb))/(rb*(rb + zb3)**2) +
     -       (yb3*(cb*rb + yb3))/(rb**3*(rb + zb3)) -
     -       (1 + (cb*yb3)/rb)/(rb*(rb + zb3))))


      dv(2,1,1) =
     .  2*y1*(1/(r*(r - y3)) + 1/(rb*(rb + yb3))) +
     -  y1**2*(-(y1/(r**2*(r - y3)**2)) - y1/(r**3*(r - y3)) -
     -     y1/(rb**2*(rb + yb3)**2) - y1/(rb**3*(rb + yb3))) -
     -  ((r*sb - y1)*(-sb + y1/r)*z1)/(r*(r - z3)**2) +
     -  (cb*(r*sb - y1))/(r*(r - z3)) -
     -  ((r*sb - y1)*y1*z1)/(r**3*(r - z3)) +
     -  ((-1 + (sb*y1)/r)*z1)/(r*(r - z3)) -
     -  ((rb*sb - y1)*(-sb + y1/rb)*zb1)/(rb*(rb + zb3)**2) +
     -  (cb*(rb*sb - y1))/(rb*(rb + zb3)) -
     -  ((rb*sb - y1)*y1*zb1)/(rb**3*(rb + zb3)) +
     -  ((-1 + (sb*y1)/rb)*zb1)/(rb*(rb + zb3)) +
     -  (-1 + 2*nu)*(y1/(r*(r - y3)) + y1/(rb*(rb + yb3)) -
     -     cb*((-sb + y1/r)/(r - z3) + (-sb + y1/rb)/(rb + zb3)))


      dv(2,1,2) =
     .  y1**2*(-(y2/(r**2*(r - y3)**2)) - y2/(r**3*(r - y3)) -
     -     y2/(rb**2*(rb + yb3)**2) - y2/(rb**3*(rb + yb3))) -
     -  ((r*sb - y1)*y2*z1)/(r**2*(r - z3)**2) +
     -  (sb*y2*z1)/(r**2*(r - z3)) -
     -  ((r*sb - y1)*y2*z1)/(r**3*(r - z3)) -
     -  ((rb*sb - y1)*y2*zb1)/(rb**2*(rb + zb3)**2) +
     -  (sb*y2*zb1)/(rb**2*(rb + zb3)) -
     -  ((rb*sb - y1)*y2*zb1)/(rb**3*(rb + zb3)) +
     -  (-1 + 2*nu)*(y2/(r*(r - y3)) + y2/(rb*(rb + yb3)) -
     -     cb*(y2/(r*(r - z3)) + y2/(rb*(rb + zb3))))


      dv(2,1,3) =
     .  y1**2*(-(y3/(r**3*(r - y3))) - (-1 + y3/r)/(r*(r - y3)**2) -
     -     yb3/(rb**3*(rb + yb3)) - (1 + yb3/rb)/(rb*(rb + yb3)**2)) -
     -  ((r*sb - y1)*(-cb + y3/r)*z1)/(r*(r - z3)**2) -
     -  (sb*(r*sb - y1))/(r*(r - z3)) + (sb*y3*z1)/(r**2*(r - z3)) -
     -  ((r*sb - y1)*y3*z1)/(r**3*(r - z3)) -
     -  ((rb*sb - y1)*(cb + yb3/rb)*zb1)/(rb*(rb + zb3)**2) +
     -  (sb*(rb*sb - y1))/(rb*(rb + zb3)) +
     -  (sb*yb3*zb1)/(rb**2*(rb + zb3)) -
     -  ((rb*sb - y1)*yb3*zb1)/(rb**3*(rb + zb3)) +
     -  (-1 + 2*nu)*((-1 + y3/r)/(r - y3) + (1 + yb3/rb)/(rb + yb3) -
     -     cb*((-cb + y3/r)/(r - z3) + (cb + yb3/rb)/(rb + zb3)))


      dv(2,2,1) =
     .  y2*(1/(r*(r - y3)) + 1/(rb*(rb + yb3))) +
     -  y1*y2*(-(y1/(r**2*(r - y3)**2)) - y1/(r**3*(r - y3)) -
     -     y1/(rb**2*(rb + yb3)**2) - y1/(rb**3*(rb + yb3))) +
     -  2*(1 - nu)*(y2*(1/(y1**2 + y2**2) - cb/(y2**2 + z1**2) +
     -        (sb*(cb*y1*(-r**2 + y2**2) + (-r**2 + y1**2)*z1))/
     -         (r*(r**2*sb**2*y2**2 + (cb*y2**2 + y1*z1)**2))) +
     -     y2*(1/(y1**2 + y2**2) - cb/(y2**2 + zb1**2) +
     -        (sb*(cb*y1*(-rb**2 + y2**2) + (-rb**2 + y1**2)*zb1))/
     -         (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2)))) -
     -  y2*(-(((-sb + y1/r)*z1)/(r*(r - z3)**2)) + cb/(r*(r - z3)) -
     -     (y1*z1)/(r**3*(r - z3)) -
     -     ((-sb + y1/rb)*zb1)/(rb*(rb + zb3)**2) +
     -     cb/(rb*(rb + zb3)) - (y1*zb1)/(rb**3*(rb + zb3)))


      dv(2,2,2) =
     .  y1*(1/(r*(r - y3)) + 1/(rb*(rb + yb3))) +
     -  y1*y2*(-(y2/(r**2*(r - y3)**2)) - y2/(r**3*(r - y3)) -
     -     y2/(rb**2*(rb + yb3)**2) - y2/(rb**3*(rb + yb3))) -
     -  z1/(r*(r - z3)) + 2*(1 - nu)*
     -   ((-2*y1)/(y1**2 + y2**2) + z1/(y2**2 + z1**2) +
     -     (sb*(cb*(-(r**2*y2**2) + y2**4) + y1*(r**2 + y2**2)*z1))/
     -      (r*(r**2*sb**2*y2**2 + (cb*y2**2 + y1*z1)**2)) +
     -     zb1/(y2**2 + zb1**2) +
     -     (sb*(cb*(-(rb**2*y2**2) + y2**4) + y1*(rb**2 + y2**2)*zb1))/
     -      (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2))) -
     -  zb1/(rb*(rb + zb3)) -
     -  y2*(-((y2*z1)/(r**2*(r - z3)**2)) - (y2*z1)/(r**3*(r - z3)) -
     -     (y2*zb1)/(rb**2*(rb + zb3)**2) - (y2*zb1)/(rb**3*(rb + zb3)))


      dv(2,2,3) =
     .  y1*y2*(-(y3/(r**3*(r - y3))) - (-1 + y3/r)/(r*(r - y3)**2) -
     -     yb3/(rb**3*(rb + yb3)) - (1 + yb3/rb)/(rb*(rb + yb3)**2)) +
     -  2*(1 - nu)*(sb*y2*(1/(y2**2 + z1**2) +
     -        (r**2*sb*y1 + cb*y2**2*y3 + y1*y3*z1)/
     -         (r*(r**2*sb**2*y2**2 + (cb*y2**2 + y1*z1)**2))) +
     -     sb*y2*(-(1/(y2**2 + zb1**2)) +
     -        (-(rb**2*sb*y1) + cb*y2**2*yb3 + y1*yb3*zb1)/
     -         (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2)))) -
     -  y2*(-(((-cb + y3/r)*z1)/(r*(r - z3)**2)) - sb/(r*(r - z3)) -
     -     (y3*z1)/(r**3*(r - z3)) -
     -     ((cb + yb3/rb)*zb1)/(rb*(rb + zb3)**2) +
     -     sb/(rb*(rb + zb3)) - (yb3*zb1)/(rb**3*(rb + zb3)))


      dv(2,3,1) =
     .  -(1/r) + 1/rb - y1*(-(y1/r**3) + y1/rb**3) -
     -  ((-sb + y1/r)*(cb*r - y3)*z1)/(r*(r - z3)**2) +
     -  (cb*(cb*r - y3))/(r*(r - z3)) + (cb*y1*z1)/(r**2*(r - z3)) -
     -  (y1*(cb*r - y3)*z1)/(r**3*(r - z3)) +
     -  ((-sb + y1/rb)*(cb*rb + yb3)*zb1)/(rb*(rb + zb3)**2) -
     -  (cb*(cb*rb + yb3))/(rb*(rb + zb3)) -
     -  (cb*y1*zb1)/(rb**2*(rb + zb3)) +
     -  (y1*(cb*rb + yb3)*zb1)/(rb**3*(rb + zb3)) +
     -  (-1 + 2*nu)*sb*((-sb + y1/r)/(r - z3) -
     -     (-sb + y1/rb)/(rb + zb3))


      dv(2,3,2) =
     .  -(y1*(-(y2/r**3) + y2/rb**3)) -
     -  (y2*(cb*r - y3)*z1)/(r**2*(r - z3)**2) +
     -  (cb*y2*z1)/(r**2*(r - z3)) -
     -  (y2*(cb*r - y3)*z1)/(r**3*(r - z3)) +
     -  (y2*(cb*rb + yb3)*zb1)/(rb**2*(rb + zb3)**2) -
     -  (cb*y2*zb1)/(rb**2*(rb + zb3)) +
     -  (y2*(cb*rb + yb3)*zb1)/(rb**3*(rb + zb3)) +
     -  (-1 + 2*nu)*sb*(y2/(r*(r - z3)) - y2/(rb*(rb + zb3)))


      dv(2,3,3) =
     .  -(y1*(-(y3/r**3) + yb3/rb**3)) -
     -  ((cb*r - y3)*(-cb + y3/r)*z1)/(r*(r - z3)**2) -
     -  (sb*(cb*r - y3))/(r*(r - z3)) -
     -  ((cb*r - y3)*y3*z1)/(r**3*(r - z3)) +
     -  ((-1 + (cb*y3)/r)*z1)/(r*(r - z3)) +
     -  ((cb*rb + yb3)*(cb + yb3/rb)*zb1)/(rb*(rb + zb3)**2) -
     -  (sb*(cb*rb + yb3))/(rb*(rb + zb3)) +
     -  (yb3*(cb*rb + yb3)*zb1)/(rb**3*(rb + zb3)) -
     -  ((1 + (cb*yb3)/rb)*zb1)/(rb*(rb + zb3)) +
     -  (-1 + 2*nu)*sb*((-cb + y3/r)/(r - z3) -
     -     (cb + yb3/rb)/(rb + zb3))


      dv(3,1,1) =
     .  sb*y2*(-(((r*sb - y1)*(-sb + y1/r))/(r*(r - z3)**2)) -
     -    ((r*sb - y1)*y1)/(r**3*(r - z3)) +
     -    (-1 + (sb*y1)/r)/(r*(r - z3)) -
     -    ((rb*sb - y1)*(-sb + y1/rb))/(rb*(rb + zb3)**2) -
     -    ((rb*sb - y1)*y1)/(rb**3*(rb + zb3)) +
     -    (-1 + (sb*y1)/rb)/(rb*(rb + zb3)))


      dv(3,1,2) =
     .  sb*((r*sb - y1)/(r*(r - z3)) + (rb*sb - y1)/(rb*(rb + zb3))) +
     -  sb*y2*(-(((r*sb - y1)*y2)/(r**2*(r - z3)**2)) +
     -     (sb*y2)/(r**2*(r - z3)) - ((r*sb - y1)*y2)/(r**3*(r - z3)) -
     -     ((rb*sb - y1)*y2)/(rb**2*(rb + zb3)**2) +
     -     (sb*y2)/(rb**2*(rb + zb3)) -
     -     ((rb*sb - y1)*y2)/(rb**3*(rb + zb3)))


      dv(3,1,3) =
     .  sb*y2*(-(((r*sb - y1)*(-cb + y3/r))/(r*(r - z3)**2)) +
     -    (sb*y3)/(r**2*(r - z3)) - ((r*sb - y1)*y3)/(r**3*(r - z3)) -
     -    ((rb*sb - y1)*(cb + yb3/rb))/(rb*(rb + zb3)**2) +
     -    (sb*yb3)/(rb**2*(rb + zb3)) -
     -    ((rb*sb - y1)*yb3)/(rb**3*(rb + zb3)))


      dv(3,2,1) =
     .  -(sb*y2**2*(-((-sb + y1/r)/(r*(r - z3)**2)) -
     -       y1/(r**3*(r - z3)) - (-sb + y1/rb)/(rb*(rb + zb3)**2) -
     -       y1/(rb**3*(rb + zb3)))) +
     -  (1 - 2*nu)*sb*((-sb + y1/r)/(r - z3) + (-sb + y1/rb)/(rb + zb3))


      dv(3,2,2) =
     .  -2*sb*y2*(1/(r*(r - z3)) + 1/(rb*(rb + zb3))) -
     -  sb*y2**2*(-(y2/(r**2*(r - z3)**2)) - y2/(r**3*(r - z3)) -
     -     y2/(rb**2*(rb + zb3)**2) - y2/(rb**3*(rb + zb3))) +
     -  (1 - 2*nu)*sb*(y2/(r*(r - z3)) + y2/(rb*(rb + zb3)))


      dv(3,2,3) =
     .  -(sb*y2**2*(-((-cb + y3/r)/(r*(r - z3)**2)) -
     -       y3/(r**3*(r - z3)) - (cb + yb3/rb)/(rb*(rb + zb3)**2) -
     -       yb3/(rb**3*(rb + zb3)))) +
     -  (1 - 2*nu)*sb*((-cb + y3/r)/(r - z3) + (cb + yb3/rb)/(rb + zb3))


      dv(3,3,1) =
     .  2*(1 - nu)*(y2*(1/(y1**2 + y2**2) - cb/(y2**2 + z1**2) +
     -        (sb*(cb*y1*(-r**2 + y2**2) + (-r**2 + y1**2)*z1))/
     -         (r*(r**2*sb**2*y2**2 + (cb*y2**2 + y1*z1)**2))) -
     -     y2*(1/(y1**2 + y2**2) - cb/(y2**2 + zb1**2) +
     -        (sb*(cb*y1*(-rb**2 + y2**2) + (-rb**2 + y1**2)*zb1))/
     -         (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2)))) +
     -  sb*y2*(-(((-sb + y1/r)*(cb*r - y3))/(r*(r - z3)**2)) +
     -     (cb*y1)/(r**2*(r - z3)) - (y1*(cb*r - y3))/(r**3*(r - z3)) +
     -     ((-sb + y1/rb)*(cb*rb + yb3))/(rb*(rb + zb3)**2) -
     -     (cb*y1)/(rb**2*(rb + zb3)) +
     -     (y1*(cb*rb + yb3))/(rb**3*(rb + zb3)))


      dv(3,3,2) =
     .  2*(1 - nu)*(z1/(y2**2 + z1**2) +
     -     (sb*(cb*(-(r**2*y2**2) + y2**4) + y1*(r**2 + y2**2)*z1))/
     -      (r*(r**2*sb**2*y2**2 + (cb*y2**2 + y1*z1)**2)) -
     -     zb1/(y2**2 + zb1**2) -
     -     (sb*(cb*(-(rb**2*y2**2) + y2**4) + y1*(rb**2 + y2**2)*zb1))/
     -      (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2))) +
     -  sb*((cb*r - y3)/(r*(r - z3)) - (cb*rb + yb3)/(rb*(rb + zb3))) +
     -  sb*y2*(-((y2*(cb*r - y3))/(r**2*(r - z3)**2)) +
     -     (cb*y2)/(r**2*(r - z3)) - (y2*(cb*r - y3))/(r**3*(r - z3)) +
     -     (y2*(cb*rb + yb3))/(rb**2*(rb + zb3)**2) -
     -     (cb*y2)/(rb**2*(rb + zb3)) +
     -     (y2*(cb*rb + yb3))/(rb**3*(rb + zb3)))


      dv(3,3,3) =
     .  2*(1 - nu)*(sb*y2*(1/(y2**2 + z1**2) +
     -        (r**2*sb*y1 + cb*y2**2*y3 + y1*y3*z1)/
     -         (r*(r**2*sb**2*y2**2 + (cb*y2**2 + y1*z1)**2))) -
     -     sb*y2*(-(1/(y2**2 + zb1**2)) +
     -        (-(rb**2*sb*y1) + cb*y2**2*yb3 + y1*yb3*zb1)/
     -         (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2)))) +
     -  sb*y2*(-(((cb*r - y3)*(-cb + y3/r))/(r*(r - z3)**2)) -
     -     ((cb*r - y3)*y3)/(r**3*(r - z3)) +
     -     (-1 + (cb*y3)/r)/(r*(r - z3)) +
     -     ((cb*rb + yb3)*(cb + yb3/rb))/(rb*(rb + zb3)**2) +
     -     (yb3*(cb*rb + yb3))/(rb**3*(rb + zb3)) -
     -     (1 + (cb*yb3)/rb)/(rb*(rb + zb3)))


c     vc   derivatives

      dvc(1,1,1) =
     .  (-3*a*ctb*y1*y2*(-a + yb3))/rb**5 -
     -  ((1 - 2*nu)*y1*y2*(ctb*(1 - 2*nu - a/rb) -
     -       ((nu + a/rb)*y1)/(rb + yb3)))/(rb*(rb + yb3)**2) -
     -  (y1*y2*(-a + yb3)*(-(ctb*(1 - 2*nu)) + (a*y1)/rb**2 +
     -       ((2*nu + a/rb)*y1)/(rb + yb3)))/(rb**2*(rb + yb3)**2) -
     -  (y1*y2*(-a + yb3)*(-(ctb*(1 - 2*nu)) + (a*y1)/rb**2 +
     -       ((2*nu + a/rb)*y1)/(rb + yb3)))/(rb**3*(rb + yb3)) +
     -  (y2*(-a + yb3)*(a/rb**2 - (2*a*y1**2)/rb**4 -
     -       ((2*nu + a/rb)*y1**2)/(rb*(rb + yb3)**2) +
     -       (2*nu + a/rb)/(rb + yb3) - (a*y1**2)/(rb**3*(rb + yb3))))/
     -   (rb*(rb + yb3)) + ((1 - 2*nu)*y2*
     -     ((a*ctb*y1)/rb**3 + ((nu + a/rb)*y1**2)/(rb*(rb + yb3)**2) -
     -       (nu + a/rb)/(rb + yb3) + (a*y1**2)/(rb**3*(rb + yb3))))/
     -   (rb + yb3) - 2*ctb**2*(1 - 2*nu)*(1 - nu)*y2*
     -   (1/(y1**2 + y2**2) - cb/(y2**2 + zb1**2) +
     -     (sb*(cb*y1*(-rb**2 + y2**2) + (-rb**2 + y1**2)*zb1))/
     -      (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2))) -
     -  (cb*ctb*(1 - 2*nu)*(cb + a/rb)*(-sb + y1/rb)*y2)/
     -   (rb + zb3)**2 - (a*cb*ctb*(1 - 2*nu)*y1*y2)/
     -   (rb**3*(rb + zb3)) -
     -  ((-sb + y1/rb)*y2*(-a + yb3)*
     -     (-((a*cb*ctb*yb3)/rb**2) +
     -       (cb*(2*cb*(1 - nu)*(rb*sb - y1) +
     -            ctb*(cb*(1 - 2*nu) - a/rb)*(cb*rb + yb3)))/(rb + zb3))
     -     )/(rb*(rb + zb3)**2) -
     -  (y1*y2*(-a + yb3)*(-((a*cb*ctb*yb3)/rb**2) +
     -       (cb*(2*cb*(1 - nu)*(rb*sb - y1) +
     -            ctb*(cb*(1 - 2*nu) - a/rb)*(cb*rb + yb3)))/(rb + zb3))
     -     )/(rb**3*(rb + zb3)) +
     -  (y2*(-a + yb3)*((2*a*cb*ctb*y1*yb3)/rb**4 -
     -       (cb*(-sb + y1/rb)*
     -          (2*cb*(1 - nu)*(rb*sb - y1) +
     -            ctb*(cb*(1 - 2*nu) - a/rb)*(cb*rb + yb3)))/
     -        (rb + zb3)**2 +
     -       (cb*((cb*ctb*(cb*(1 - 2*nu) - a/rb)*y1)/rb +
     -            2*cb*(1 - nu)*(-1 + (sb*y1)/rb) +
     -            (a*ctb*y1*(cb*rb + yb3))/rb**3))/(rb + zb3)))/
     -   (rb*(rb + zb3))

      dvc(1,1,2) =
     .  (a*ctb*(-a + yb3))/rb**3 - (3*a*ctb*y2**2*(-a + yb3))/rb**5 -
     -  ((1 - 2*nu)*y2**2*(ctb*(1 - 2*nu - a/rb) -
     -       ((nu + a/rb)*y1)/(rb + yb3)))/(rb*(rb + yb3)**2) +
     -  ((1 - 2*nu)*(ctb*(1 - 2*nu - a/rb) -
     -       ((nu + a/rb)*y1)/(rb + yb3)))/(rb + yb3) -
     -  (y2**2*(-a + yb3)*(-(ctb*(1 - 2*nu)) + (a*y1)/rb**2 +
     -       ((2*nu + a/rb)*y1)/(rb + yb3)))/(rb**2*(rb + yb3)**2) +
     -  ((-a + yb3)*(-(ctb*(1 - 2*nu)) + (a*y1)/rb**2 +
     -       ((2*nu + a/rb)*y1)/(rb + yb3)))/(rb*(rb + yb3)) -
     -  (y2**2*(-a + yb3)*(-(ctb*(1 - 2*nu)) + (a*y1)/rb**2 +
     -       ((2*nu + a/rb)*y1)/(rb + yb3)))/(rb**3*(rb + yb3)) +
     -  (y2*(-a + yb3)*((-2*a*y1*y2)/rb**4 -
     -       ((2*nu + a/rb)*y1*y2)/(rb*(rb + yb3)**2) -
     -       (a*y1*y2)/(rb**3*(rb + yb3))))/(rb*(rb + yb3)) +
     -  ((1 - 2*nu)*y2*((a*ctb*y2)/rb**3 +
     -       ((nu + a/rb)*y1*y2)/(rb*(rb + yb3)**2) +
     -       (a*y1*y2)/(rb**3*(rb + yb3))))/(rb + yb3) -
     -  2*ctb**2*(1 - 2*nu)*(1 - nu)*
     -   (-(y1/(y1**2 + y2**2)) + zb1/(y2**2 + zb1**2) +
     -     (sb*(cb*(-(rb**2*y2**2) + y2**4) + y1*(rb**2 + y2**2)*zb1))/
     -      (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2))) -
     -  (cb*ctb*(1 - 2*nu)*(cb + a/rb)*y2**2)/(rb*(rb + zb3)**2) +
     -  (cb*ctb*(1 - 2*nu)*(cb + a/rb))/(rb + zb3) -
     -  (a*cb*ctb*(1 - 2*nu)*y2**2)/(rb**3*(rb + zb3)) -
     -  (y2**2*(-a + yb3)*(-((a*cb*ctb*yb3)/rb**2) +
     -       (cb*(2*cb*(1 - nu)*(rb*sb - y1) +
     -            ctb*(cb*(1 - 2*nu) - a/rb)*(cb*rb + yb3)))/(rb + zb3))
     -     )/(rb**2*(rb + zb3)**2) +
     -  ((-a + yb3)*(-((a*cb*ctb*yb3)/rb**2) +
     -       (cb*(2*cb*(1 - nu)*(rb*sb - y1) +
     -            ctb*(cb*(1 - 2*nu) - a/rb)*(cb*rb + yb3)))/(rb + zb3))
     -     )/(rb*(rb + zb3)) -
     -  (y2**2*(-a + yb3)*(-((a*cb*ctb*yb3)/rb**2) +
     -       (cb*(2*cb*(1 - nu)*(rb*sb - y1) +
     -            ctb*(cb*(1 - 2*nu) - a/rb)*(cb*rb + yb3)))/(rb + zb3))
     -     )/(rb**3*(rb + zb3)) +
     -  (y2*(-a + yb3)*((2*a*cb*ctb*y2*yb3)/rb**4 -
     -       (cb*y2*(2*cb*(1 - nu)*(rb*sb - y1) +
     -            ctb*(cb*(1 - 2*nu) - a/rb)*(cb*rb + yb3)))/
     -        (rb*(rb + zb3)**2) +
     -       (cb*((cb*ctb*(cb*(1 - 2*nu) - a/rb)*y2)/rb +
     -            (2*cb*(1 - nu)*sb*y2)/rb +
     -            (a*ctb*y2*(cb*rb + yb3))/rb**3))/(rb + zb3)))/
     -   (rb*(rb + zb3))

      dvc(1,1,3) =
     .  (a*ctb*y2)/rb**3 - (3*a*ctb*y2*yb3*(-a + yb3))/rb**5 -
     -  ((1 - 2*nu)*y2*(1 + yb3/rb)*
     -     (ctb*(1 - 2*nu - a/rb) - ((nu + a/rb)*y1)/(rb + yb3)))/
     -   (rb + yb3)**2 + (y2*
     -     (-(ctb*(1 - 2*nu)) + (a*y1)/rb**2 +
     -       ((2*nu + a/rb)*y1)/(rb + yb3)))/(rb*(rb + yb3)) -
     -  (y2*yb3*(-a + yb3)*(-(ctb*(1 - 2*nu)) + (a*y1)/rb**2 +
     -       ((2*nu + a/rb)*y1)/(rb + yb3)))/(rb**3*(rb + yb3)) -
     -  (y2*(-a + yb3)*(1 + yb3/rb)*
     -     (-(ctb*(1 - 2*nu)) + (a*y1)/rb**2 +
     -       ((2*nu + a/rb)*y1)/(rb + yb3)))/(rb*(rb + yb3)**2) +
     -  ((1 - 2*nu)*y2*((a*ctb*yb3)/rb**3 +
     -       (a*y1*yb3)/(rb**3*(rb + yb3)) +
     -       ((nu + a/rb)*y1*(1 + yb3/rb))/(rb + yb3)**2))/(rb + yb3) +
     -  (y2*(-a + yb3)*((-2*a*y1*yb3)/rb**4 -
     -       (a*y1*yb3)/(rb**3*(rb + yb3)) -
     -       ((2*nu + a/rb)*y1*(1 + yb3/rb))/(rb + yb3)**2))/
     -   (rb*(rb + yb3)) - 2*ctb**2*(1 - 2*nu)*(1 - nu)*sb*y2*
     -   (-(1/(y2**2 + zb1**2)) +
     -     (-(rb**2*sb*y1) + cb*y2**2*yb3 + y1*yb3*zb1)/
     -      (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2))) -
     -  (cb*ctb*(1 - 2*nu)*(cb + a/rb)*y2*(cb + yb3/rb))/
     -   (rb + zb3)**2 - (a*cb*ctb*(1 - 2*nu)*y2*yb3)/
     -   (rb**3*(rb + zb3)) -
     -  (y2*(-a + yb3)*(cb + yb3/rb)*
     -     (-((a*cb*ctb*yb3)/rb**2) +
     -       (cb*(2*cb*(1 - nu)*(rb*sb - y1) +
     -            ctb*(cb*(1 - 2*nu) - a/rb)*(cb*rb + yb3)))/(rb + zb3))
     -     )/(rb*(rb + zb3)**2) +
     -  (y2*(-((a*cb*ctb*yb3)/rb**2) +
     -       (cb*(2*cb*(1 - nu)*(rb*sb - y1) +
     -            ctb*(cb*(1 - 2*nu) - a/rb)*(cb*rb + yb3)))/(rb + zb3))
     -     )/(rb*(rb + zb3)) -
     -  (y2*yb3*(-a + yb3)*(-((a*cb*ctb*yb3)/rb**2) +
     -       (cb*(2*cb*(1 - nu)*(rb*sb - y1) +
     -            ctb*(cb*(1 - 2*nu) - a/rb)*(cb*rb + yb3)))/(rb + zb3))
     -     )/(rb**3*(rb + zb3)) +
     -  (y2*(-a + yb3)*(-((a*cb*ctb)/rb**2) +
     -       (2*a*cb*ctb*yb3**2)/rb**4 -
     -       (cb*(cb + yb3/rb)*
     -          (2*cb*(1 - nu)*(rb*sb - y1) +
     -            ctb*(cb*(1 - 2*nu) - a/rb)*(cb*rb + yb3)))/
     -        (rb + zb3)**2 +
     -       (cb*((2*cb*(1 - nu)*sb*yb3)/rb +
     -            (a*ctb*yb3*(cb*rb + yb3))/rb**3 +
     -            ctb*(cb*(1 - 2*nu) - a/rb)*(1 + (cb*yb3)/rb)))/
     -        (rb + zb3)))/(rb*(rb + zb3))



      dvc(1,2,1) =
     .  -((a*ctb*(-a + yb3))/rb**3) +
     -  (3*a*ctb*y1**2*(-a + yb3))/rb**5 +
     -  ((1 - 2*nu)*y1*(-a + ctb*(1 - 2*nu - a/rb)*y1 + nu*yb3 +
     -       ((nu + a/rb)*y2**2)/(rb + yb3)))/(rb*(rb + yb3)**2) -
     -  (y1*(-a + yb3)*(-2*nu + (-a + ctb*(1 - 2*nu)*y1)/rb +
     -       (a*y2**2)/rb**3 + ((2*nu + a/rb)*y2**2)/(rb*(rb + yb3))))/
     -   (rb*(rb + yb3)**2) -
     -  ((1 - 2*nu)*(ctb*(1 - 2*nu - a/rb) + (a*ctb*y1**2)/rb**3 -
     -       ((nu + a/rb)*y1*y2**2)/(rb*(rb + yb3)**2) -
     -       (a*y1*y2**2)/(rb**3*(rb + yb3))))/(rb + yb3) +
     -  ((-a + yb3)*((ctb*(1 - 2*nu))/rb -
     -       (y1*(-a + ctb*(1 - 2*nu)*y1))/rb**3 -
     -       (3*a*y1*y2**2)/rb**5 -
     -       ((2*nu + a/rb)*y1*y2**2)/(rb**2*(rb + yb3)**2) -
     -       (a*y1*y2**2)/(rb**4*(rb + yb3)) -
     -       ((2*nu + a/rb)*y1*y2**2)/(rb**3*(rb + yb3))))/(rb + yb3) +
     -  (ctb*(1 - 2*nu)*(cb + a/rb)*(-sb + y1/rb)*zb1)/(rb + zb3)**2 -
     -  (cb*ctb*(1 - 2*nu)*(cb + a/rb))/(rb + zb3) +
     -  (a*ctb*(1 - 2*nu)*y1*zb1)/(rb**3*(rb + zb3)) +
     -  (1 - 2*nu)*(((2*ctb**2*(1 - nu) - nu)*y1)/(rb*(rb + yb3)) -
     -     (cb*(1 + 2*ctb**2*(1 - nu) - 2*nu)*(-sb + y1/rb))/(rb + zb3))
     -    - ((-sb + y1/rb)*(-a + yb3)*
     -     (cb**2 + (a*ctb*yb3*zb1)/rb**3 -
     -       (a*cb + ctb*(1 - 2*nu)*zb1)/rb -
     -       (cb**2*y2**2 - (a*ctb*(cb*rb + yb3)*zb1)/rb)/
     -        (rb*(rb + zb3))))/(rb + zb3)**2 +
     -  ((-a + yb3)*(-((cb*ctb*(1 - 2*nu))/rb) + (a*cb*ctb*yb3)/rb**3 -
     -       (3*a*ctb*y1*yb3*zb1)/rb**5 +
     -       (y1*(a*cb + ctb*(1 - 2*nu)*zb1))/rb**3 +
     -       ((-sb + y1/rb)*(cb**2*y2**2 -
     -            (a*ctb*(cb*rb + yb3)*zb1)/rb))/(rb*(rb + zb3)**2) +
     -       (y1*(cb**2*y2**2 - (a*ctb*(cb*rb + yb3)*zb1)/rb))/
     -        (rb**3*(rb + zb3)) -
     -       (-((a*cb*ctb*(cb*rb + yb3))/rb) -
     -          (a*cb*ctb*y1*zb1)/rb**2 +
     -          (a*ctb*y1*(cb*rb + yb3)*zb1)/rb**3)/(rb*(rb + zb3))))/
     -   (rb + zb3)



      dvc(1,2,2) =
     .  (3*a*ctb*y1*y2*(-a + yb3))/rb**5 +
     -  ((1 - 2*nu)*y2*(-a + ctb*(1 - 2*nu - a/rb)*y1 + nu*yb3 +
     -       ((nu + a/rb)*y2**2)/(rb + yb3)))/(rb*(rb + yb3)**2) -
     -  (y2*(-a + yb3)*(-2*nu + (-a + ctb*(1 - 2*nu)*y1)/rb +
     -       (a*y2**2)/rb**3 + ((2*nu + a/rb)*y2**2)/(rb*(rb + yb3))))/
     -   (rb*(rb + yb3)**2) -
     -  ((1 - 2*nu)*((a*ctb*y1*y2)/rb**3 -
     -       ((nu + a/rb)*y2**3)/(rb*(rb + yb3)**2) +
     -       (2*(nu + a/rb)*y2)/(rb + yb3) -
     -       (a*y2**3)/(rb**3*(rb + yb3))))/(rb + yb3) +
     -  ((-a + yb3)*((2*a*y2)/rb**3 -
     -       ((-a + ctb*(1 - 2*nu)*y1)*y2)/rb**3 - (3*a*y2**3)/rb**5 -
     -       ((2*nu + a/rb)*y2**3)/(rb**2*(rb + yb3)**2) +
     -       (2*(2*nu + a/rb)*y2)/(rb*(rb + yb3)) -
     -       (a*y2**3)/(rb**4*(rb + yb3)) -
     -       ((2*nu + a/rb)*y2**3)/(rb**3*(rb + yb3))))/(rb + yb3) +
     -  (ctb*(1 - 2*nu)*(cb + a/rb)*y2*zb1)/(rb*(rb + zb3)**2) +
     -  (a*ctb*(1 - 2*nu)*y2*zb1)/(rb**3*(rb + zb3)) +
     -  (1 - 2*nu)*(((2*ctb**2*(1 - nu) - nu)*y2)/(rb*(rb + yb3)) -
     -     (cb*(1 + 2*ctb**2*(1 - nu) - 2*nu)*y2)/(rb*(rb + zb3))) -
     -  (y2*(-a + yb3)*(cb**2 + (a*ctb*yb3*zb1)/rb**3 -
     -       (a*cb + ctb*(1 - 2*nu)*zb1)/rb -
     -       (cb**2*y2**2 - (a*ctb*(cb*rb + yb3)*zb1)/rb)/
     -        (rb*(rb + zb3))))/(rb*(rb + zb3)**2) +
     -  ((-a + yb3)*((-3*a*ctb*y2*yb3*zb1)/rb**5 +
     -       (y2*(a*cb + ctb*(1 - 2*nu)*zb1))/rb**3 +
     -       (y2*(cb**2*y2**2 - (a*ctb*(cb*rb + yb3)*zb1)/rb))/
     -        (rb**2*(rb + zb3)**2) +
     -       (y2*(cb**2*y2**2 - (a*ctb*(cb*rb + yb3)*zb1)/rb))/
     -        (rb**3*(rb + zb3)) -
     -       (2*cb**2*y2 - (a*cb*ctb*y2*zb1)/rb**2 +
     -          (a*ctb*y2*(cb*rb + yb3)*zb1)/rb**3)/(rb*(rb + zb3))))/
     -   (rb + zb3)


      dvc(1,2,3) =
     .  -((a*ctb*y1)/rb**3) + (3*a*ctb*y1*yb3*(-a + yb3))/rb**5 +
     -  ((1 - 2*nu)*(1 + yb3/rb)*
     -     (-a + ctb*(1 - 2*nu - a/rb)*y1 + nu*yb3 +
     -       ((nu + a/rb)*y2**2)/(rb + yb3)))/(rb + yb3)**2 +
     -  (-2*nu + (-a + ctb*(1 - 2*nu)*y1)/rb + (a*y2**2)/rb**3 +
     -     ((2*nu + a/rb)*y2**2)/(rb*(rb + yb3)))/(rb + yb3) -
     -  ((-a + yb3)*(1 + yb3/rb)*
     -     (-2*nu + (-a + ctb*(1 - 2*nu)*y1)/rb + (a*y2**2)/rb**3 +
     -       ((2*nu + a/rb)*y2**2)/(rb*(rb + yb3))))/(rb + yb3)**2 -
     -  ((1 - 2*nu)*(nu + (a*ctb*y1*yb3)/rb**3 -
     -       (a*y2**2*yb3)/(rb**3*(rb + yb3)) -
     -       ((nu + a/rb)*y2**2*(1 + yb3/rb))/(rb + yb3)**2))/(rb + yb3)
     -    + ((-a + yb3)*(-(((-a + ctb*(1 - 2*nu)*y1)*yb3)/rb**3) -
     -       (3*a*y2**2*yb3)/rb**5 - (a*y2**2*yb3)/(rb**4*(rb + yb3)) -
     -       ((2*nu + a/rb)*y2**2*yb3)/(rb**3*(rb + yb3)) -
     -       ((2*nu + a/rb)*y2**2*(1 + yb3/rb))/(rb*(rb + yb3)**2)))/
     -   (rb + yb3) + (ctb*(1 - 2*nu)*(cb + a/rb)*(cb + yb3/rb)*zb1)/
     -   (rb + zb3)**2 - (ctb*(1 - 2*nu)*(cb + a/rb)*sb)/(rb + zb3) +
     -  (a*ctb*(1 - 2*nu)*yb3*zb1)/(rb**3*(rb + zb3)) +
     -  (1 - 2*nu)*(((2*ctb**2*(1 - nu) - nu)*(1 + yb3/rb))/
     -      (rb + yb3) - (cb*(1 + 2*ctb**2*(1 - nu) - 2*nu)*
     -        (cb + yb3/rb))/(rb + zb3)) -
     -  ((-a + yb3)*(cb + yb3/rb)*
     -     (cb**2 + (a*ctb*yb3*zb1)/rb**3 -
     -       (a*cb + ctb*(1 - 2*nu)*zb1)/rb -
     -       (cb**2*y2**2 - (a*ctb*(cb*rb + yb3)*zb1)/rb)/
     -        (rb*(rb + zb3))))/(rb + zb3)**2 +
     -  (cb**2 + (a*ctb*yb3*zb1)/rb**3 -
     -     (a*cb + ctb*(1 - 2*nu)*zb1)/rb -
     -     (cb**2*y2**2 - (a*ctb*(cb*rb + yb3)*zb1)/rb)/(rb*(rb + zb3)))
     -    /(rb + zb3) + ((-a + yb3)*
     -     (-((ctb*(1 - 2*nu)*sb)/rb) + (a*ctb*sb*yb3)/rb**3 +
     -       (a*ctb*zb1)/rb**3 - (3*a*ctb*yb3**2*zb1)/rb**5 +
     -       (yb3*(a*cb + ctb*(1 - 2*nu)*zb1))/rb**3 +
     -       ((cb + yb3/rb)*(cb**2*y2**2 -
     -            (a*ctb*(cb*rb + yb3)*zb1)/rb))/(rb*(rb + zb3)**2) +
     -       (yb3*(cb**2*y2**2 - (a*ctb*(cb*rb + yb3)*zb1)/rb))/
     -        (rb**3*(rb + zb3)) -
     -       (-((a*ctb*sb*(cb*rb + yb3))/rb) +
     -          (a*ctb*yb3*(cb*rb + yb3)*zb1)/rb**3 -
     -          (a*ctb*(1 + (cb*yb3)/rb)*zb1)/rb)/(rb*(rb + zb3))))/
     -   (rb + zb3)



      dvc(1,3,1) =
     .  (y2*(-a + yb3)*((-2*a*y1)/rb**4 -
     -       (2*nu*y1)/(rb*(rb + yb3)**2)))/rb -
     -  (y1*y2*(-a + yb3)*(a/rb**2 + (2*nu)/(rb + yb3)))/rb**3 +
     -  2*(1 - nu)*(-(((2*nu + a/rb)*y1*y2)/(rb*(rb + yb3)**2)) -
     -     (a*y1*y2)/(rb**3*(rb + yb3)) +
     -     ctb*(1 - 2*nu)*y2*
     -      (1/(y1**2 + y2**2) - cb/(y2**2 + zb1**2) +
     -        (sb*(cb*y1*(-rb**2 + y2**2) + (-rb**2 + y1**2)*zb1))/
     -         (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2))) +
     -     (cb*(cb + a/rb)*(-sb + y1/rb)*y2)/(rb + zb3)**2 +
     -     (a*cb*y1*y2)/(rb**3*(rb + zb3))) -
     -  (cb*(-sb + y1/rb)*y2*(-a + yb3)*
     -     (1 - 2*nu - (a*yb3)/rb**2 -
     -       ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/(rb*(rb + zb3)**2)
     -    - (cb*y1*y2*(-a + yb3)*
     -     (1 - 2*nu - (a*yb3)/rb**2 -
     -       ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/(rb**3*(rb + zb3))
     -    + (cb*y2*(-a + yb3)*
     -     ((2*a*y1*yb3)/rb**4 +
     -       ((cb + a/rb)*(-sb + y1/rb)*(cb*rb + yb3))/(rb + zb3)**2 -
     -       (cb*(cb + a/rb)*y1)/(rb*(rb + zb3)) +
     -       (a*y1*(cb*rb + yb3))/(rb**3*(rb + zb3))))/(rb*(rb + zb3))

      dvc(1,3,2) =
     .  (y2*(-a + yb3)*((-2*a*y2)/rb**4 -
     -       (2*nu*y2)/(rb*(rb + yb3)**2)))/rb +
     -  ((-a + yb3)*(a/rb**2 + (2*nu)/(rb + yb3)))/rb -
     -  (y2**2*(-a + yb3)*(a/rb**2 + (2*nu)/(rb + yb3)))/rb**3 +
     -  2*(1 - nu)*(-(((2*nu + a/rb)*y2**2)/(rb*(rb + yb3)**2)) +
     -     (2*nu + a/rb)/(rb + yb3) - (a*y2**2)/(rb**3*(rb + yb3)) +
     -     ctb*(1 - 2*nu)*(-(y1/(y1**2 + y2**2)) +
     -        zb1/(y2**2 + zb1**2) +
     -        (sb*(cb*(-(rb**2*y2**2) + y2**4) +
     -             y1*(rb**2 + y2**2)*zb1))/
     -         (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2))) +
     -     (cb*(cb + a/rb)*y2**2)/(rb*(rb + zb3)**2) -
     -     (cb*(cb + a/rb))/(rb + zb3) + (a*cb*y2**2)/(rb**3*(rb + zb3))
     -     ) - (cb*y2**2*(-a + yb3)*
     -     (1 - 2*nu - (a*yb3)/rb**2 -
     -       ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb**2*(rb + zb3)**2) +
     -  (cb*(-a + yb3)*(1 - 2*nu - (a*yb3)/rb**2 -
     -       ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/(rb*(rb + zb3)) -
     -  (cb*y2**2*(-a + yb3)*
     -     (1 - 2*nu - (a*yb3)/rb**2 -
     -       ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/(rb**3*(rb + zb3))
     -    + (cb*y2*(-a + yb3)*
     -     ((2*a*y2*yb3)/rb**4 +
     -       ((cb + a/rb)*y2*(cb*rb + yb3))/(rb*(rb + zb3)**2) -
     -       (cb*(cb + a/rb)*y2)/(rb*(rb + zb3)) +
     -       (a*y2*(cb*rb + yb3))/(rb**3*(rb + zb3))))/(rb*(rb + zb3))



      dvc(1,3,3) =
     .  (y2*(a/rb**2 + (2*nu)/(rb + yb3)))/rb -
     -  (y2*yb3*(-a + yb3)*(a/rb**2 + (2*nu)/(rb + yb3)))/rb**3 +
     -  (y2*(-a + yb3)*((-2*a*yb3)/rb**4 -
     -       (2*nu*(1 + yb3/rb))/(rb + yb3)**2))/rb +
     -  2*(1 - nu)*(-((a*y2*yb3)/(rb**3*(rb + yb3))) -
     -     ((2*nu + a/rb)*y2*(1 + yb3/rb))/(rb + yb3)**2 +
     -     ctb*(1 - 2*nu)*sb*y2*
     -      (-(1/(y2**2 + zb1**2)) +
     -        (-(rb**2*sb*y1) + cb*y2**2*yb3 + y1*yb3*zb1)/
     -         (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2))) +
     -     (cb*(cb + a/rb)*y2*(cb + yb3/rb))/(rb + zb3)**2 +
     -     (a*cb*y2*yb3)/(rb**3*(rb + zb3))) -
     -  (cb*y2*(-a + yb3)*(cb + yb3/rb)*
     -     (1 - 2*nu - (a*yb3)/rb**2 -
     -       ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/(rb*(rb + zb3)**2)
     -    + (cb*y2*(1 - 2*nu - (a*yb3)/rb**2 -
     -       ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/(rb*(rb + zb3)) -
     -  (cb*y2*yb3*(-a + yb3)*
     -     (1 - 2*nu - (a*yb3)/rb**2 -
     -       ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/(rb**3*(rb + zb3))
     -    + (cb*y2*(-a + yb3)*
     -     (-(a/rb**2) + (2*a*yb3**2)/rb**4 +
     -       ((cb + a/rb)*(cb*rb + yb3)*(cb + yb3/rb))/(rb + zb3)**2 +
     -       (a*yb3*(cb*rb + yb3))/(rb**3*(rb + zb3)) -
     -       ((cb + a/rb)*(1 + (cb*yb3)/rb))/(rb + zb3)))/
     -   (rb*(rb + zb3))


      dvc(2,1,1) =
     .  -((a*ctb*(-a + yb3))/rb**3) +
     -  (3*a*ctb*y1**2*(-a + yb3))/rb**5 -
     -  ((1 - 2*nu)*y1*(-a + ctb*(-1 + 2*nu)*y1 + (a*ctb*y1)/rb +
     -       nu*yb3 + ((nu + a/rb)*y1**2)/(rb + yb3)))/
     -   (rb*(rb + yb3)**2) -
     -  (y1*(-a + yb3)*(2*nu - (a*y1**2)/rb**3 +
     -       (a + ctb*(1 - 2*nu)*y1)/rb -
     -  ((2*nu + a/rb)*y1**2)/(rb*(rb + yb3))))/(rb*(rb + yb3)**2)
     -   + ((1 - 2*nu)*(ctb*(-1 + 2*nu) + (a*ctb)/rb -
     -       (a*ctb*y1**2)/rb**3 -
     -       ((nu + a/rb)*y1**3)/(rb*(rb + yb3)**2) +
     -       (2*(nu + a/rb)*y1)/(rb + yb3) -
     -       (a*y1**3)/(rb**3*(rb + yb3))))/(rb + yb3) +
     -  ((-a + yb3)*((ctb*(1 - 2*nu))/rb - (2*a*y1)/rb**3 +
     -       (3*a*y1**3)/rb**5 - (y1*(a + ctb*(1 - 2*nu)*y1))/rb**3 +
     -       ((2*nu + a/rb)*y1**3)/(rb**2*(rb + yb3)**2) -
     -       (2*(2*nu + a/rb)*y1)/(rb*(rb + yb3)) +
     -       (a*y1**3)/(rb**4*(rb + yb3)) +
     -       ((2*nu + a/rb)*y1**3)/(rb**3*(rb + yb3))))/(rb + yb3) +
     -  (ctb*(1 - 2*nu)*(-sb + y1/rb)*
     -     (-((a*(rb*sb - y1))/(cb*rb)) + cb*zb1))/(rb + zb3)**2 -
     -  (ctb*(1 - 2*nu)*(cb**2 + (a*(rb*sb - y1)*y1)/(cb*rb**3) -
     -       (a*(-1 + (sb*y1)/rb))/(cb*rb)))/(rb + zb3) +
     -  (1 - 2*nu)*(((2*ctb**2*(1 - nu) + nu)*y1)/(rb*(rb + yb3)) -
     -     (cb*(1 + 2*ctb**2*(1 - nu))*(-sb + y1/rb))/(rb + zb3)) -
     -  (ctb*(-sb + y1/rb)*(-a + yb3)*
     -     (-(cb*sb) + (a*y1*yb3)/(cb*rb**3) +
     -       ((rb*sb - y1)*(2*cb*(1 - nu) -
     -            ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/rb))/
     -   (rb + zb3)**2 + (ctb*(-a + yb3)*
     -     ((a*yb3)/(cb*rb**3) - (3*a*y1**2*yb3)/(cb*rb**5) -
     -       ((rb*sb - y1)*y1*
     -          (2*cb*(1 - nu) -
     -  ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/rb**3 +
     -       ((-1 + (sb*y1)/rb)*
     -          (2*cb*(1 - nu) -
     -            ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/rb +
     -       ((rb*sb - y1)*(((1 + a/(cb*rb))*(-sb + y1/rb)*
     -               (cb*rb + yb3))/(rb + zb3)**2 -
     -            (cb*(1 + a/(cb*rb))*y1)/(rb*(rb + zb3)) +
     -            (a*y1*(cb*rb + yb3))/(cb*rb**3*(rb + zb3))))/rb))/
     -   (rb + zb3)


      dvc(2,1,2) =
     .  (3*a*ctb*y1*y2*(-a + yb3))/rb**5 -
     -  ((1 - 2*nu)*y2*(-a + ctb*(-1 + 2*nu)*y1 + (a*ctb*y1)/rb +
     -       nu*yb3 + ((nu + a/rb)*y1**2)/(rb + yb3)))/
     -   (rb*(rb + yb3)**2) -
     -  (y2*(-a + yb3)*(2*nu - (a*y1**2)/rb**3 +
     -       (a + ctb*(1 - 2*nu)*y1)/rb -
     - ((2*nu + a/rb)*y1**2)/(rb*(rb + yb3))))/(rb*(rb + yb3)**2)
     -   + ((1 - 2*nu)*(-((a*ctb*y1*y2)/rb**3) -
     -       ((nu + a/rb)*y1**2*y2)/(rb*(rb + yb3)**2) -
     -       (a*y1**2*y2)/(rb**3*(rb + yb3))))/(rb + yb3) +
     -  ((-a + yb3)*((3*a*y1**2*y2)/rb**5 -
     -       ((a + ctb*(1 - 2*nu)*y1)*y2)/rb**3 +
     -       ((2*nu + a/rb)*y1**2*y2)/(rb**2*(rb + yb3)**2) +
     -       (a*y1**2*y2)/(rb**4*(rb + yb3)) +
     - ((2*nu + a/rb)*y1**2*y2)/(rb**3*(rb + yb3))))/(rb + yb3) +
     -  (ctb*(1 - 2*nu)*y2*(-((a*(rb*sb - y1))/(cb*rb)) + cb*zb1))/
     -   (rb*(rb + zb3)**2) -
     -  (ctb*(1 - 2*nu)*(-((a*sb*y2)/(cb*rb**2)) +
     -       (a*(rb*sb - y1)*y2)/(cb*rb**3)))/(rb + zb3) +
     -  (1 - 2*nu)*(((2*ctb**2*(1 - nu) + nu)*y2)/(rb*(rb + yb3)) -
     -     (cb*(1 + 2*ctb**2*(1 - nu))*y2)/(rb*(rb + zb3))) -
     -  (ctb*y2*(-a + yb3)*(-(cb*sb) + (a*y1*yb3)/(cb*rb**3) +
     -       ((rb*sb - y1)*(2*cb*(1 - nu) -
     -            ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/rb))/
     -   (rb*(rb + zb3)**2) +
     -  (ctb*(-a + yb3)*((-3*a*y1*y2*yb3)/(cb*rb**5) +
     -       (sb*y2*(2*cb*(1 - nu) -
     -            ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/rb**2 -
     -       ((rb*sb - y1)*y2*
     -          (2*cb*(1 - nu) -
     -            ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/rb**3 +
     -       ((rb*sb - y1)*(((1 + a/(cb*rb))*y2*(cb*rb + yb3))/
     -             (rb*(rb + zb3)**2) -
     -            (cb*(1 + a/(cb*rb))*y2)/(rb*(rb + zb3)) +
     -            (a*y2*(cb*rb + yb3))/(cb*rb**3*(rb + zb3))))/rb))/
     -   (rb + zb3)

      dvc(2,1,3) =
     .  -((a*ctb*y1)/rb**3) + (3*a*ctb*y1*yb3*(-a + yb3))/rb**5 -
     -  ((1 - 2*nu)*(1 + yb3/rb)*
     -     (-a + ctb*(-1 + 2*nu)*y1 + (a*ctb*y1)/rb + nu*yb3 +
     -       ((nu + a/rb)*y1**2)/(rb + yb3)))/(rb + yb3)**2 +
     -  (2*nu - (a*y1**2)/rb**3 + (a + ctb*(1 - 2*nu)*y1)/rb -
     -     ((2*nu + a/rb)*y1**2)/(rb*(rb + yb3)))/(rb + yb3) -
     -  ((-a + yb3)*(1 + yb3/rb)*
     -     (2*nu - (a*y1**2)/rb**3 + (a + ctb*(1 - 2*nu)*y1)/rb -
     -       ((2*nu + a/rb)*y1**2)/(rb*(rb + yb3))))/(rb + yb3)**2 +
     -  ((1 - 2*nu)*(nu - (a*ctb*y1*yb3)/rb**3 -
     -       (a*y1**2*yb3)/(rb**3*(rb + yb3)) -
     -       ((nu + a/rb)*y1**2*(1 + yb3/rb))/(rb + yb3)**2))/(rb + yb3)
     -    + ((-a + yb3)*((3*a*y1**2*yb3)/rb**5 -
     -       ((a + ctb*(1 - 2*nu)*y1)*yb3)/rb**3 +
     -       (a*y1**2*yb3)/(rb**4*(rb + yb3)) +
     -       ((2*nu + a/rb)*y1**2*yb3)/(rb**3*(rb + yb3)) +
     -       ((2*nu + a/rb)*y1**2*(1 + yb3/rb))/(rb*(rb + yb3)**2)))/
     -   (rb + yb3) + (ctb*(1 - 2*nu)*(cb + yb3/rb)*
     -     (-((a*(rb*sb - y1))/(cb*rb)) + cb*zb1))/(rb + zb3)**2 -
     -  (ctb*(1 - 2*nu)*(cb*sb - (a*sb*yb3)/(cb*rb**2) +
     -       (a*(rb*sb - y1)*yb3)/(cb*rb**3)))/(rb + zb3) +
     -  (1 - 2*nu)*(((2*ctb**2*(1 - nu) + nu)*(1 + yb3/rb))/
     -      (rb + yb3) - (cb*(1 + 2*ctb**2*(1 - nu))*(cb + yb3/rb))/
     -      (rb + zb3)) - (ctb*(-a + yb3)*(cb + yb3/rb)*
     -     (-(cb*sb) + (a*y1*yb3)/(cb*rb**3) +
     -       ((rb*sb - y1)*(2*cb*(1 - nu) -
     -            ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/rb))/
     -   (rb + zb3)**2 + (ctb*
     -     (-(cb*sb) + (a*y1*yb3)/(cb*rb**3) +
     -       ((rb*sb - y1)*(2*cb*(1 - nu) -
     -            ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/rb))/
     -   (rb + zb3) + (ctb*(-a + yb3)*
     -     ((a*y1)/(cb*rb**3) - (3*a*y1*yb3**2)/(cb*rb**5) +
     -       (sb*yb3*(2*cb*(1 - nu) -
     -            ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/rb**2 -
     -       ((rb*sb - y1)*yb3*
     -          (2*cb*(1 - nu) -
     -            ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/rb**3 +
     -       ((rb*sb - y1)*(((1 + a/(cb*rb))*(cb*rb + yb3)*
     -               (cb + yb3/rb))/(rb + zb3)**2 +
     -            (a*yb3*(cb*rb + yb3))/(cb*rb**3*(rb + zb3)) -
     -            ((1 + a/(cb*rb))*(1 + (cb*yb3)/rb))/(rb + zb3)))/rb))/
     -   (rb + zb3)

      dvc(2,2,1) =
     .  (3*a*ctb*y1*y2*(-a + yb3))/rb**5 -
     -  ((1 - 2*nu)*y1*y2*(ctb*(-1 + 2*nu + a/rb) +
     -       ((nu + a/rb)*y1)/(rb + yb3)))/(rb*(rb + yb3)**2) +
     -  ((1 - 2*nu)*y2*(-((a*ctb*y1)/rb**3) -
     -       ((nu + a/rb)*y1**2)/(rb*(rb + yb3)**2) +
     -       (nu + a/rb)/(rb + yb3) - (a*y1**2)/(rb**3*(rb + yb3))))/
     -   (rb + yb3) - (y1*y2*(-a + yb3)*
     -     (ctb*(1 - 2*nu) - (2*nu*y1)/(rb + yb3) -
     -       (a*y1*(1/rb + 1/(rb + yb3)))/rb))/(rb**2*(rb + yb3)**2) -
     -  (y1*y2*(-a + yb3)*(ctb*(1 - 2*nu) - (2*nu*y1)/(rb + yb3) -
     -       (a*y1*(1/rb + 1/(rb + yb3)))/rb))/(rb**3*(rb + yb3)) +
     -  (y2*(-a + yb3)*((2*nu*y1**2)/(rb*(rb + yb3)**2) -
     -       (2*nu)/(rb + yb3) -
     -       (a*y1*(-(y1/rb**3) - y1/(rb*(rb + yb3)**2)))/rb -
     -       (a*(1/rb + 1/(rb + yb3)))/rb +
     -       (a*y1**2*(1/rb + 1/(rb + yb3)))/rb**3))/(rb*(rb + yb3)) +
     -  2*ctb**2*(1 - 2*nu)*(1 - nu)*y2*
     -   (1/(y1**2 + y2**2) - cb/(y2**2 + zb1**2) +
     -     (sb*(cb*y1*(-rb**2 + y2**2) + (-rb**2 + y1**2)*zb1))/
     -      (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2))) +
     -  (ctb*(1 - 2*nu)*(1 + a/(cb*rb))*(-sb + y1/rb)*y2)/
     -   (rb + zb3)**2 + (a*ctb*(1 - 2*nu)*y1*y2)/
     -   (cb*rb**3*(rb + zb3)) -
     -  (ctb*(-sb + y1/rb)*y2*(-a + yb3)*
     -     (-2*cb*(1 - nu) + (a*yb3)/(cb*rb**2) +
     -       ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb*(rb + zb3)**2) -
     -  (ctb*y1*y2*(-a + yb3)*
     -     (-2*cb*(1 - nu) + (a*yb3)/(cb*rb**2) +
     -       ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb**3*(rb + zb3)) +
     -  (ctb*y2*(-a + yb3)*((-2*a*y1*yb3)/(cb*rb**4) -
     -       ((1 + a/(cb*rb))*(-sb + y1/rb)*(cb*rb + yb3))/
     -        (rb + zb3)**2 + (cb*(1 + a/(cb*rb))*y1)/(rb*(rb + zb3)) -
     -       (a*y1*(cb*rb + yb3))/(cb*rb**3*(rb + zb3))))/
     -   (rb*(rb + zb3))

      dvc(2,2,2) =
     .  -((a*ctb*(-a + yb3))/rb**3) +
     -  (3*a*ctb*y2**2*(-a + yb3))/rb**5 -
     -  ((1 - 2*nu)*y2**2*(ctb*(-1 + 2*nu + a/rb) +
     -       ((nu + a/rb)*y1)/(rb + yb3)))/(rb*(rb + yb3)**2) +
     -  ((1 - 2*nu)*(ctb*(-1 + 2*nu + a/rb) +
     -       ((nu + a/rb)*y1)/(rb + yb3)))/(rb + yb3) +
     -  ((1 - 2*nu)*y2*(-((a*ctb*y2)/rb**3) -
     -       ((nu + a/rb)*y1*y2)/(rb*(rb + yb3)**2) -
     -       (a*y1*y2)/(rb**3*(rb + yb3))))/(rb + yb3) -
     -  (y2**2*(-a + yb3)*(ctb*(1 - 2*nu) - (2*nu*y1)/(rb + yb3) -
     -       (a*y1*(1/rb + 1/(rb + yb3)))/rb))/(rb**2*(rb + yb3)**2) +
     -  ((-a + yb3)*(ctb*(1 - 2*nu) - (2*nu*y1)/(rb + yb3) -
     -       (a*y1*(1/rb + 1/(rb + yb3)))/rb))/(rb*(rb + yb3)) -
     -  (y2**2*(-a + yb3)*(ctb*(1 - 2*nu) - (2*nu*y1)/(rb + yb3) -
     -       (a*y1*(1/rb + 1/(rb + yb3)))/rb))/(rb**3*(rb + yb3)) +
     -  (y2*(-a + yb3)*((2*nu*y1*y2)/(rb*(rb + yb3)**2) -
     -       (a*y1*(-(y2/rb**3) - y2/(rb*(rb + yb3)**2)))/rb +
     -       (a*y1*y2*(1/rb + 1/(rb + yb3)))/rb**3))/(rb*(rb + yb3)) +
     -  2*ctb**2*(1 - 2*nu)*(1 - nu)*
     -   (-(y1/(y1**2 + y2**2)) + zb1/(y2**2 + zb1**2) +
     -     (sb*(cb*(-(rb**2*y2**2) + y2**4) + y1*(rb**2 + y2**2)*zb1))/
     -      (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2))) +
     -  (ctb*(1 - 2*nu)*(1 + a/(cb*rb))*y2**2)/(rb*(rb + zb3)**2) -
     -  (ctb*(1 - 2*nu)*(1 + a/(cb*rb)))/(rb + zb3) +
     -  (a*ctb*(1 - 2*nu)*y2**2)/(cb*rb**3*(rb + zb3)) -
     -  (ctb*y2**2*(-a + yb3)*
     -     (-2*cb*(1 - nu) + (a*yb3)/(cb*rb**2) +
     -       ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb**2*(rb + zb3)**2) +
     -  (ctb*(-a + yb3)*(-2*cb*(1 - nu) + (a*yb3)/(cb*rb**2) +
     -       ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb*(rb + zb3)) - (ctb*y2**2*(-a + yb3)*
     -     (-2*cb*(1 - nu) + (a*yb3)/(cb*rb**2) +
     -       ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb**3*(rb + zb3)) +
     -  (ctb*y2*(-a + yb3)*((-2*a*y2*yb3)/(cb*rb**4) -
     -       ((1 + a/(cb*rb))*y2*(cb*rb + yb3))/(rb*(rb + zb3)**2) +
     -       (cb*(1 + a/(cb*rb))*y2)/(rb*(rb + zb3)) -
     -       (a*y2*(cb*rb + yb3))/(cb*rb**3*(rb + zb3))))/
     -   (rb*(rb + zb3))

      dvc(2,2,3) =
     .  -((a*ctb*y2)/rb**3) + (3*a*ctb*y2*yb3*(-a + yb3))/rb**5 -
     -  ((1 - 2*nu)*y2*(1 + yb3/rb)*
     -     (ctb*(-1 + 2*nu + a/rb) + ((nu + a/rb)*y1)/(rb + yb3)))/
     -   (rb + yb3)**2 + ((1 - 2*nu)*y2*
     -     (-((a*ctb*yb3)/rb**3) - (a*y1*yb3)/(rb**3*(rb + yb3)) -
     -       ((nu + a/rb)*y1*(1 + yb3/rb))/(rb + yb3)**2))/(rb + yb3) +
     -  (y2*(ctb*(1 - 2*nu) - (2*nu*y1)/(rb + yb3) -
     -       (a*y1*(1/rb + 1/(rb + yb3)))/rb))/(rb*(rb + yb3)) -
     -  (y2*yb3*(-a + yb3)*(ctb*(1 - 2*nu) - (2*nu*y1)/(rb + yb3) -
     -       (a*y1*(1/rb + 1/(rb + yb3)))/rb))/(rb**3*(rb + yb3)) -
     -  (y2*(-a + yb3)*(1 + yb3/rb)*
     -     (ctb*(1 - 2*nu) - (2*nu*y1)/(rb + yb3) -
     -       (a*y1*(1/rb + 1/(rb + yb3)))/rb))/(rb*(rb + yb3)**2) +
     -  (y2*(-a + yb3)*((2*nu*y1*(1 + yb3/rb))/(rb + yb3)**2 +
     -       (a*y1*yb3*(1/rb + 1/(rb + yb3)))/rb**3 -
     -       (a*y1*(-(yb3/rb**3) - (1 + yb3/rb)/(rb + yb3)**2))/rb))/
     -   (rb*(rb + yb3)) + 2*ctb**2*(1 - 2*nu)*(1 - nu)*sb*y2*
     -   (-(1/(y2**2 + zb1**2)) +
     -     (-(rb**2*sb*y1) + cb*y2**2*yb3 + y1*yb3*zb1)/
     -      (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2))) +
     -  (ctb*(1 - 2*nu)*(1 + a/(cb*rb))*y2*(cb + yb3/rb))/
     -   (rb + zb3)**2 + (a*ctb*(1 - 2*nu)*y2*yb3)/
     -   (cb*rb**3*(rb + zb3)) -
     -  (ctb*y2*(-a + yb3)*(cb + yb3/rb)*
     -     (-2*cb*(1 - nu) + (a*yb3)/(cb*rb**2) +
     -       ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb*(rb + zb3)**2) +
     -  (ctb*y2*(-2*cb*(1 - nu) + (a*yb3)/(cb*rb**2) +
     -       ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb*(rb + zb3)) - (ctb*y2*yb3*(-a + yb3)*
     -     (-2*cb*(1 - nu) + (a*yb3)/(cb*rb**2) +
     -       ((1 + a/(cb*rb))*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb**3*(rb + zb3)) +
     -  (ctb*y2*(-a + yb3)*(a/(cb*rb**2) - (2*a*yb3**2)/(cb*rb**4) -
     -       ((1 + a/(cb*rb))*(cb*rb + yb3)*(cb + yb3/rb))/
     -        (rb + zb3)**2 -
     -       (a*yb3*(cb*rb + yb3))/(cb*rb**3*(rb + zb3)) +
     -       ((1 + a/(cb*rb))*(1 + (cb*yb3)/rb))/(rb + zb3)))/
     -   (rb*(rb + zb3))


      dvc(2,3,1) =
     .  (2*(1 - nu)*(2*nu + a/rb)*y1**2)/(rb*(rb + yb3)**2) -
     -  (2*(1 - nu)*(2*nu + a/rb))/(rb + yb3) +
     -  (2*a*(1 - nu)*y1**2)/(rb**3*(rb + yb3)) +
     -  ((-a + yb3)*(-(a/rb**2) + (2*a*y1**2)/rb**4 +
     -       (2*nu*y1**2)/(rb*(rb + yb3)**2) - (2*nu)/(rb + yb3)))/rb -
     -  (y1*(-a + yb3)*(ctb*(1 - 2*nu) - (a*y1)/rb**2 -
     -       (2*nu*y1)/(rb + yb3)))/rb**3 -
     -  (2*(1 - nu)*(cb + a/rb)*(-sb + y1/rb)*zb1)/(rb + zb3)**2 +
     -  (2*cb*(1 - nu)*(cb + a/rb))/(rb + zb3) -
     -  (2*a*(1 - nu)*y1*zb1)/(rb**3*(rb + zb3)) -
     -  2*ctb*(1 - 2*nu)*(1 - nu)*
     -   (y1/(rb*(rb + yb3)) - (cb*(-sb + y1/rb))/(rb + zb3)) +
     -  ((-sb + y1/rb)*(-a + yb3)*
     -     (cb*sb + (ctb*(cb*rb + yb3)*
     -          (2*cb*(1 - nu) - (cb*rb + yb3)/(rb + zb3)))/rb +
     -       (a*(sb - (yb3*zb1)/rb**2 -
     -            ((cb*rb + yb3)*zb1)/(rb*(rb + zb3))))/rb))/
     -   (rb + zb3)**2 - ((-a + yb3)*
     -     ((ctb*(cb*rb + yb3)*
     -          (((-sb + y1/rb)*(cb*rb + yb3))/(rb + zb3)**2 -
     -            (cb*y1)/(rb*(rb + zb3))))/rb +
     -       (cb*ctb*y1*(2*cb*(1 - nu) - (cb*rb + yb3)/(rb + zb3)))/
     -        rb**2 - (ctb*y1*(cb*rb + yb3)*
     -          (2*cb*(1 - nu) - (cb*rb + yb3)/(rb + zb3)))/rb**3 -
     -       (a*y1*(sb - (yb3*zb1)/rb**2 -
     -            ((cb*rb + yb3)*zb1)/(rb*(rb + zb3))))/rb**3 +
     -       (a*(-((cb*yb3)/rb**2) + (2*y1*yb3*zb1)/rb**4 +
     -            ((-sb + y1/rb)*(cb*rb + yb3)*zb1)/
     -             (rb*(rb + zb3)**2) -
     -            (cb*(cb*rb + yb3))/(rb*(rb + zb3)) -
     -            (cb*y1*zb1)/(rb**2*(rb + zb3)) +
     -            (y1*(cb*rb + yb3)*zb1)/(rb**3*(rb + zb3))))/rb))/
     -   (rb + zb3)

      dvc(2,3,2) =
     .  (2*(1 - nu)*(2*nu + a/rb)*y1*y2)/(rb*(rb + yb3)**2) +
     -  (2*a*(1 - nu)*y1*y2)/(rb**3*(rb + yb3)) +
     -  ((-a + yb3)*((2*a*y1*y2)/rb**4 +
     -       (2*nu*y1*y2)/(rb*(rb + yb3)**2)))/rb -
     -  (y2*(-a + yb3)*(ctb*(1 - 2*nu) - (a*y1)/rb**2 -
     -       (2*nu*y1)/(rb + yb3)))/rb**3 -
     -  (2*(1 - nu)*(cb + a/rb)*y2*zb1)/(rb*(rb + zb3)**2) -
     -  (2*a*(1 - nu)*y2*zb1)/(rb**3*(rb + zb3)) -
     -  2*ctb*(1 - 2*nu)*(1 - nu)*
     -   (y2/(rb*(rb + yb3)) - (cb*y2)/(rb*(rb + zb3))) +
     -  (y2*(-a + yb3)*(cb*sb +
     -       (ctb*(cb*rb + yb3)*
     -          (2*cb*(1 - nu) - (cb*rb + yb3)/(rb + zb3)))/rb +
     -       (a*(sb - (yb3*zb1)/rb**2 -
     -            ((cb*rb + yb3)*zb1)/(rb*(rb + zb3))))/rb))/
     -   (rb*(rb + zb3)**2) -
     -  ((-a + yb3)*((ctb*(cb*rb + yb3)*
     -          ((y2*(cb*rb + yb3))/(rb*(rb + zb3)**2) -
     -            (cb*y2)/(rb*(rb + zb3))))/rb +
     -       (cb*ctb*y2*(2*cb*(1 - nu) - (cb*rb + yb3)/(rb + zb3)))/
     -        rb**2 - (ctb*y2*(cb*rb + yb3)*
     -          (2*cb*(1 - nu) - (cb*rb + yb3)/(rb + zb3)))/rb**3 -
     -       (a*y2*(sb - (yb3*zb1)/rb**2 -
     -            ((cb*rb + yb3)*zb1)/(rb*(rb + zb3))))/rb**3 +
     -       (a*((2*y2*yb3*zb1)/rb**4 +
     -            (y2*(cb*rb + yb3)*zb1)/(rb**2*(rb + zb3)**2) -
     -            (cb*y2*zb1)/(rb**2*(rb + zb3)) +
     -            (y2*(cb*rb + yb3)*zb1)/(rb**3*(rb + zb3))))/rb))/
     -   (rb + zb3)



      dvc(2,3,3) =
     .  (2*a*(1 - nu)*y1*yb3)/(rb**3*(rb + yb3)) +
     -  (2*(1 - nu)*(2*nu + a/rb)*y1*(1 + yb3/rb))/(rb + yb3)**2 +
     -  (ctb*(1 - 2*nu) - (a*y1)/rb**2 - (2*nu*y1)/(rb + yb3))/rb -
     -  (yb3*(-a + yb3)*(ctb*(1 - 2*nu) - (a*y1)/rb**2 -
     -       (2*nu*y1)/(rb + yb3)))/rb**3 +
     -  ((-a + yb3)*((2*a*y1*yb3)/rb**4 +
     -       (2*nu*y1*(1 + yb3/rb))/(rb + yb3)**2))/rb -
     -  (2*(1 - nu)*(cb + a/rb)*(cb + yb3/rb)*zb1)/(rb + zb3)**2 +
     -  (2*(1 - nu)*(cb + a/rb)*sb)/(rb + zb3) -
     -  (2*a*(1 - nu)*yb3*zb1)/(rb**3*(rb + zb3)) -
     -  2*ctb*(1 - 2*nu)*(1 - nu)*
     -   ((1 + yb3/rb)/(rb + yb3) - (cb*(cb + yb3/rb))/(rb + zb3)) +
     -  ((-a + yb3)*(cb + yb3/rb)*
     -     (cb*sb + (ctb*(cb*rb + yb3)*
     -          (2*cb*(1 - nu) - (cb*rb + yb3)/(rb + zb3)))/rb +
     -       (a*(sb - (yb3*zb1)/rb**2 -
     -            ((cb*rb + yb3)*zb1)/(rb*(rb + zb3))))/rb))/
     -   (rb + zb3)**2 - (cb*sb +
     -     (ctb*(cb*rb + yb3)*
     -        (2*cb*(1 - nu) - (cb*rb + yb3)/(rb + zb3)))/rb +
     -     (a*(sb - (yb3*zb1)/rb**2 -
     -          ((cb*rb + yb3)*zb1)/(rb*(rb + zb3))))/rb)/(rb + zb3) -
     -  ((-a + yb3)*(-((ctb*yb3*(cb*rb + yb3)*
     -            (2*cb*(1 - nu) - (cb*rb + yb3)/(rb + zb3)))/rb**3) +
     -       (ctb*(1 + (cb*yb3)/rb)*
     -          (2*cb*(1 - nu) - (cb*rb + yb3)/(rb + zb3)))/rb +
     -       (ctb*(cb*rb + yb3)*
     -          (((cb*rb + yb3)*(cb + yb3/rb))/(rb + zb3)**2 -
     -            (1 + (cb*yb3)/rb)/(rb + zb3)))/rb -
     -       (a*yb3*(sb - (yb3*zb1)/rb**2 -
     -            ((cb*rb + yb3)*zb1)/(rb*(rb + zb3))))/rb**3 +
     -       (a*(-((sb*yb3)/rb**2) - zb1/rb**2 + (2*yb3**2*zb1)/rb**4 +
     -            ((cb*rb + yb3)*(cb + yb3/rb)*zb1)/
     -             (rb*(rb + zb3)**2) -
     -            (sb*(cb*rb + yb3))/(rb*(rb + zb3)) +
     -            (yb3*(cb*rb + yb3)*zb1)/(rb**3*(rb + zb3)) -
     -            ((1 + (cb*yb3)/rb)*zb1)/(rb*(rb + zb3))))/rb))/
     -   (rb + zb3)


      dvc(3,1,1) =
     .  -((y2*(-a + yb3)*((-2*a*y1)/rb**4 - y1/(rb*(rb + yb3)**2)))/
     -     rb) + (y1*y2*(-a + yb3)*(a/rb**2 + 1/(rb + yb3)))/rb**3 +
     -  (1 - 2*nu)*(-(((1 + a/rb)*y1*y2)/(rb*(rb + yb3)**2)) -
     -     (a*y1*y2)/(rb**3*(rb + yb3)) +
     -     (cb*(cb + a/rb)*(-sb + y1/rb)*y2)/(rb + zb3)**2 +
     -     (a*cb*y1*y2)/(rb**3*(rb + zb3))) -
     -  (cb*(-sb + y1/rb)*y2*(-a + yb3)*
     -     ((a*yb3)/rb**2 + ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb*(rb + zb3)**2) -
     -  (cb*y1*y2*(-a + yb3)*
     -     ((a*yb3)/rb**2 + ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb**3*(rb + zb3)) +
     -  (cb*y2*(-a + yb3)*((-2*a*y1*yb3)/rb**4 -
     -       ((cb + a/rb)*(-sb + y1/rb)*(cb*rb + yb3))/(rb + zb3)**2 +
     -       (cb*(cb + a/rb)*y1)/(rb*(rb + zb3)) -
     -       (a*y1*(cb*rb + yb3))/(rb**3*(rb + zb3))))/(rb*(rb + zb3))


      dvc(3,1,2) =
     .  -((y2*(-a + yb3)*((-2*a*y2)/rb**4 - y2/(rb*(rb + yb3)**2)))/
     -     rb) - ((-a + yb3)*(a/rb**2 + 1/(rb + yb3)))/rb +
     -  (y2**2*(-a + yb3)*(a/rb**2 + 1/(rb + yb3)))/rb**3 +
     -  (1 - 2*nu)*(-(((1 + a/rb)*y2**2)/(rb*(rb + yb3)**2)) +
     -     (1 + a/rb)/(rb + yb3) - (a*y2**2)/(rb**3*(rb + yb3)) +
     -     (cb*(cb + a/rb)*y2**2)/(rb*(rb + zb3)**2) -
     -     (cb*(cb + a/rb))/(rb + zb3) + (a*cb*y2**2)/(rb**3*(rb + zb3))
     -     ) - (cb*y2**2*(-a + yb3)*
     -     ((a*yb3)/rb**2 + ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb**2*(rb + zb3)**2) +
     -  (cb*(-a + yb3)*((a*yb3)/rb**2 +
     -       ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/(rb*(rb + zb3)) -
     -  (cb*y2**2*(-a + yb3)*
     -     ((a*yb3)/rb**2 + ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb**3*(rb + zb3)) +
     -  (cb*y2*(-a + yb3)*((-2*a*y2*yb3)/rb**4 -
     -       ((cb + a/rb)*y2*(cb*rb + yb3))/(rb*(rb + zb3)**2) +
     -       (cb*(cb + a/rb)*y2)/(rb*(rb + zb3)) -
     -       (a*y2*(cb*rb + yb3))/(rb**3*(rb + zb3))))/(rb*(rb + zb3))


      dvc(3,1,3) =
     .  -((y2*(a/rb**2 + 1/(rb + yb3)))/rb) +
     -  (y2*yb3*(-a + yb3)*(a/rb**2 + 1/(rb + yb3)))/rb**3 -
     -  (y2*(-a + yb3)*((-2*a*yb3)/rb**4 - (1 + yb3/rb)/(rb + yb3)**2))/
     -   rb + (1 - 2*nu)*(-((a*y2*yb3)/(rb**3*(rb + yb3))) -
     -     ((1 + a/rb)*y2*(1 + yb3/rb))/(rb + yb3)**2 +
     -     (cb*(cb + a/rb)*y2*(cb + yb3/rb))/(rb + zb3)**2 +
     -     (a*cb*y2*yb3)/(rb**3*(rb + zb3))) -
     -  (cb*y2*(-a + yb3)*(cb + yb3/rb)*
     -     ((a*yb3)/rb**2 + ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb*(rb + zb3)**2) +
     -  (cb*y2*((a*yb3)/rb**2 +
     -       ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/(rb*(rb + zb3)) -
     -  (cb*y2*yb3*(-a + yb3)*
     -     ((a*yb3)/rb**2 + ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/
     -   (rb**3*(rb + zb3)) +
     -  (cb*y2*(-a + yb3)*(a/rb**2 - (2*a*yb3**2)/rb**4 -
     -       ((cb + a/rb)*(cb*rb + yb3)*(cb + yb3/rb))/(rb + zb3)**2 -
     -       (a*yb3*(cb*rb + yb3))/(rb**3*(rb + zb3)) +
     -       ((cb + a/rb)*(1 + (cb*yb3)/rb))/(rb + zb3)))/
     -   (rb*(rb + zb3))


      dvc(3,2,1) =
     .  (y1*(-a + yb3)*((-2*a*y1)/rb**4 - y1/(rb*(rb + yb3)**2)))/rb +
     -  ((-a + yb3)*(a/rb**2 + 1/(rb + yb3)))/rb -
     -  (y1**2*(-a + yb3)*(a/rb**2 + 1/(rb + yb3)))/rb**3 +
     -  (1 - 2*nu)*(((1 + a/rb)*y1**2)/(rb*(rb + yb3)**2) -
     -     (1 + a/rb)/(rb + yb3) + (a*y1**2)/(rb**3*(rb + yb3)) -
     -     ((cb + a/rb)*(-sb + y1/rb)*zb1)/(rb + zb3)**2 +
     -     (cb*(cb + a/rb))/(rb + zb3) -
     -     (sb*(-sb + y1/rb))/(rb + zb3) - (a*y1*zb1)/(rb**3*(rb + zb3))
     -     ) + ((-sb + y1/rb)*(-a + yb3)*
     -     ((cb - a/rb)*sb + ((1 + (a*yb3)/rb**2)*zb1)/rb -
     -       (cb*sb*y2**2 - (a*(cb*rb + yb3)*zb1)/rb)/(rb*(rb + zb3))))/
     -   (rb + zb3)**2 - ((-a + yb3)*
     -     ((a*sb*y1)/rb**3 + (cb*(1 + (a*yb3)/rb**2))/rb -
     -       (2*a*y1*yb3*zb1)/rb**5 -
     -       (y1*(1 + (a*yb3)/rb**2)*zb1)/rb**3 +
     -       ((-sb + y1/rb)*(cb*sb*y2**2 - (a*(cb*rb + yb3)*zb1)/rb))/
     -        (rb*(rb + zb3)**2) +
     -       (y1*(cb*sb*y2**2 - (a*(cb*rb + yb3)*zb1)/rb))/
     -        (rb**3*(rb + zb3)) -
     -       (-((a*cb*(cb*rb + yb3))/rb) - (a*cb*y1*zb1)/rb**2 +
     -          (a*y1*(cb*rb + yb3)*zb1)/rb**3)/(rb*(rb + zb3))))/
     -   (rb + zb3)


      dvc(3,2,2) =
     .  (y1*(-a + yb3)*((-2*a*y2)/rb**4 - y2/(rb*(rb + yb3)**2)))/rb -
     -  (y1*y2*(-a + yb3)*(a/rb**2 + 1/(rb + yb3)))/rb**3 +
     -  (1 - 2*nu)*(((1 + a/rb)*y1*y2)/(rb*(rb + yb3)**2) +
     -     (a*y1*y2)/(rb**3*(rb + yb3)) -
     -     ((cb + a/rb)*y2*zb1)/(rb*(rb + zb3)**2) -
     -     (sb*y2)/(rb*(rb + zb3)) - (a*y2*zb1)/(rb**3*(rb + zb3))) +
     -  (y2*(-a + yb3)*((cb - a/rb)*sb + ((1 + (a*yb3)/rb**2)*zb1)/rb -
     -       (cb*sb*y2**2 - (a*(cb*rb + yb3)*zb1)/rb)/(rb*(rb + zb3))))/
     -   (rb*(rb + zb3)**2) -
     -  ((-a + yb3)*((a*sb*y2)/rb**3 - (2*a*y2*yb3*zb1)/rb**5 -
     -       (y2*(1 + (a*yb3)/rb**2)*zb1)/rb**3 +
     -       (y2*(cb*sb*y2**2 - (a*(cb*rb + yb3)*zb1)/rb))/
     -        (rb**2*(rb + zb3)**2) +
     -       (y2*(cb*sb*y2**2 - (a*(cb*rb + yb3)*zb1)/rb))/
     -        (rb**3*(rb + zb3)) -
     -       (2*cb*sb*y2 - (a*cb*y2*zb1)/rb**2 +
     -          (a*y2*(cb*rb + yb3)*zb1)/rb**3)/(rb*(rb + zb3))))/
     -   (rb + zb3)

      dvc(3,2,3) =
     .  (y1*(a/rb**2 + 1/(rb + yb3)))/rb -
     -  (y1*yb3*(-a + yb3)*(a/rb**2 + 1/(rb + yb3)))/rb**3 +
     -  (y1*(-a + yb3)*((-2*a*yb3)/rb**4 - (1 + yb3/rb)/(rb + yb3)**2))/
     -   rb + (1 - 2*nu)*((a*y1*yb3)/(rb**3*(rb + yb3)) +
     -     ((1 + a/rb)*y1*(1 + yb3/rb))/(rb + yb3)**2 -
     -     ((cb + a/rb)*(cb + yb3/rb)*zb1)/(rb + zb3)**2 +
     -     ((cb + a/rb)*sb)/(rb + zb3) -
     -     (sb*(cb + yb3/rb))/(rb + zb3) -
     -     (a*yb3*zb1)/(rb**3*(rb + zb3))) +
     -  ((-a + yb3)*(cb + yb3/rb)*
     -     ((cb - a/rb)*sb + ((1 + (a*yb3)/rb**2)*zb1)/rb -
     -       (cb*sb*y2**2 - (a*(cb*rb + yb3)*zb1)/rb)/(rb*(rb + zb3))))/
     -   (rb + zb3)**2 - ((cb - a/rb)*sb +
     -     ((1 + (a*yb3)/rb**2)*zb1)/rb -
     -     (cb*sb*y2**2 - (a*(cb*rb + yb3)*zb1)/rb)/(rb*(rb + zb3)))/
     -   (rb + zb3) - ((-a + yb3)*
     -     ((a*sb*yb3)/rb**3 + (sb*(1 + (a*yb3)/rb**2))/rb -
     -       (yb3*(1 + (a*yb3)/rb**2)*zb1)/rb**3 +
     -       ((a/rb**2 - (2*a*yb3**2)/rb**4)*zb1)/rb +
     -       ((cb + yb3/rb)*(cb*sb*y2**2 - (a*(cb*rb + yb3)*zb1)/rb))/
     -        (rb*(rb + zb3)**2) +
     -       (yb3*(cb*sb*y2**2 - (a*(cb*rb + yb3)*zb1)/rb))/
     -        (rb**3*(rb + zb3)) -
     -       (-((a*sb*(cb*rb + yb3))/rb) +
     -          (a*yb3*(cb*rb + yb3)*zb1)/rb**3 -
     -          (a*(1 + (cb*yb3)/rb)*zb1)/rb)/(rb*(rb + zb3))))/
     -   (rb + zb3)

      dvc(3,3,1) =
     .  2*(1 - nu)*(y2*(1/(y1**2 + y2**2) - cb/(y2**2 + zb1**2) +
     -        (sb*(cb*y1*(-rb**2 + y2**2) + (-rb**2 + y1**2)*zb1))/
     -         (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2))) -
     -     ((cb + a/rb)*sb*(-sb + y1/rb)*y2)/(rb + zb3)**2 -
     -     (a*sb*y1*y2)/(rb**3*(rb + zb3))) -
     -  (sb*(-sb + y1/rb)*y2*(-a + yb3)*
     -     (1 + (a*yb3)/rb**2 + ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))
     -    /(rb*(rb + zb3)**2) -
     -  (sb*y1*y2*(-a + yb3)*
     -     (1 + (a*yb3)/rb**2 + ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))
     -    /(rb**3*(rb + zb3)) +
     -  (sb*y2*(-a + yb3)*((-2*a*y1*yb3)/rb**4 -
     -       ((cb + a/rb)*(-sb + y1/rb)*(cb*rb + yb3))/(rb + zb3)**2 +
     -       (cb*(cb + a/rb)*y1)/(rb*(rb + zb3)) -
     -       (a*y1*(cb*rb + yb3))/(rb**3*(rb + zb3))))/(rb*(rb + zb3))


      dvc(3,3,2) =
     .  2*(1 - nu)*(-(y1/(y1**2 + y2**2)) + zb1/(y2**2 + zb1**2) +
     -     (sb*(cb*(-(rb**2*y2**2) + y2**4) + y1*(rb**2 + y2**2)*zb1))/
     -      (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2)) -
     -     ((cb + a/rb)*sb*y2**2)/(rb*(rb + zb3)**2) +
     -     ((cb + a/rb)*sb)/(rb + zb3) - (a*sb*y2**2)/(rb**3*(rb + zb3))
     -     ) - (sb*y2**2*(-a + yb3)*
     -     (1 + (a*yb3)/rb**2 + ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))
     -    /(rb**2*(rb + zb3)**2) +
     -  (sb*(-a + yb3)*(1 + (a*yb3)/rb**2 +
     -       ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/(rb*(rb + zb3)) -
     -  (sb*y2**2*(-a + yb3)*
     -     (1 + (a*yb3)/rb**2 + ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))
     -    /(rb**3*(rb + zb3)) +
     -  (sb*y2*(-a + yb3)*((-2*a*y2*yb3)/rb**4 -
     -       ((cb + a/rb)*y2*(cb*rb + yb3))/(rb*(rb + zb3)**2) +
     -       (cb*(cb + a/rb)*y2)/(rb*(rb + zb3)) -
     -       (a*y2*(cb*rb + yb3))/(rb**3*(rb + zb3))))/(rb*(rb + zb3))

      dvc(3,3,3) =
     .  2*(1 - nu)*(sb*y2*(-(1/(y2**2 + zb1**2)) +
     -        (-(rb**2*sb*y1) + cb*y2**2*yb3 + y1*yb3*zb1)/
     -         (rb*(rb**2*sb**2*y2**2 + (cb*y2**2 + y1*zb1)**2))) -
     -     ((cb + a/rb)*sb*y2*(cb + yb3/rb))/(rb + zb3)**2 -
     -     (a*sb*y2*yb3)/(rb**3*(rb + zb3))) -
     -  (sb*y2*(-a + yb3)*(cb + yb3/rb)*
     -     (1 + (a*yb3)/rb**2 + ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))
     -    /(rb*(rb + zb3)**2) +
     -  (sb*y2*(1 + (a*yb3)/rb**2 +
     -       ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))/(rb*(rb + zb3)) -
     -  (sb*y2*yb3*(-a + yb3)*
     -     (1 + (a*yb3)/rb**2 + ((cb + a/rb)*(cb*rb + yb3))/(rb + zb3)))
     -    /(rb**3*(rb + zb3)) +
     -  (sb*y2*(-a + yb3)*(a/rb**2 - (2*a*yb3**2)/rb**4 -
     -       ((cb + a/rb)*(cb*rb + yb3)*(cb + yb3/rb))/(rb + zb3)**2 -
     -       (a*yb3*(cb*rb + yb3))/(rb**3*(rb + zb3)) +
     -       ((cb + a/rb)*(1 + (cb*yb3)/rb))/(rb + zb3)))/
     -   (rb*(rb + zb3))


      pi   = 4.d0*atan2(1.d0,1.d0)
      vfac = 1.d0/(8.d0*pi*(1.d0-nu))

      do 20  i=1,3
      do 20  j=1,3
      v(i,j)    = vfac*( v(i,j)    + 2.d0*vc(i,j))
      do 10 k=1,3
      dv(i,j,k) = vfac*( dv(i,j,k) + 2.d0*dvc(i,j,k))
   10 continue
   20 continue

      return
      end

c---+-|--1----+----2----+----3----+----4----+----5----+----6----+----7--
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
