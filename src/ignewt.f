c
c  Original date of RCS archived version is  May 21  1998
c
c
c     $Id: ignewt.f,v 1.3 2012/03/01 21:36:12 agnew Exp agnew $
c
c     $Log: ignewt.f,v $
c     Revision 1.3  2012/03/01 21:36:12  agnew
c     for gravity, added case when height is zero
c
c     Revision 1.2  2012/03/01 21:32:23  agnew
c     gravity calculation modified to use a single formula
c     both gravity and tilt calculations converted to double precision
c
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine ignewt(del,stp,grav,tilt,pot)
c  returns the integrated green function for newtonian gravity,
c tilt, and potential.
c  the function is the integral over the interval cntered at del with
c width stp (ie, the interval [del-stp/2,del+stp/2],
c del and stp both being in radians
c  the height correction is included in the green functions,
c the station height in meters being passed as ht in the common block
c stloc
c
c  $$$$$$calls only system routines
c
      double precision eps,eps1,eps2,s,gt,c1
      save iflg,eps,eps1,eps2,g2,em,plc
      common/stloc/ct,st,rlam,ht
      data gn/6.67e-11/,a/6.371e6/,iflg/0/
c  on first call do part of height correction, and find constant parts
      if(iflg.eq.0) then
         iflg = 1
         eps = ht/a
         eps1=1.+eps
         eps2=eps*eps
         g2 = gn/(eps1*eps1)
	 g = 9.7803327*(1+.005279*ct*ct) - 3.08e-6*ht
	 em = gn/g
	 plc = 4*a*em
      endif
c sign for gravity makes + acceleration up
      if(eps.ne.0) then
        s=dsin((del+stp/2)/2.d0)
        gt=(2.d0*eps1*s**2-eps)/dsqrt(4*eps1*s**2+eps2)
        s=dsin((del-stp/2)/2.d0)
        grav=gt-(2.d0*eps1*s**2-eps)/dsqrt(4*eps1*s**2+eps2)
        grav=-g2*grav
      endif
      if(eps.eq.0) then
        grav=-g2*(dsin((del+stp/2)/2.d0)-dsin((del-stp/2)/2.d0))
      endif
c for tilt, use different approximations if < or > 0.05 rad (2.9 deg)
      c1 = dsin(0.25d0*stp)
      if(del.ge.0.05) then
        tilt = em*(-2.*dsin(del/2.d0)*c1 +
     1    dlog(dtan((del+0.5*stp)/4.d0)/dtan((del-0.5*stp)/4.d0)))
      endif
      if(del.lt.0.05) then
         d1 = del + stp/2
         d2 = del - stp/2
         tilt = em*(d2/dsqrt(eps2+d2*d2) - d1/dsqrt(eps2+d1*d1) + 
     1      dlog( (d1+dsqrt(eps2+d1*d1))/(d2+dsqrt(eps2+d2*d2)) ) )
      endif
      pot=plc*cos(del/2.)*c1
      return
      end
