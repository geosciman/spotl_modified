c
c  Original date of RCS archived version is  Jan  4  1995
c
c
c     $Id: invspt.f,v 1.1 2011/09/23 22:39:01 agnew Exp agnew $
c
c     $Log: invspt.f,v $
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine invspt(alp,del,b,rlong)
c  solves the inverse spherical triangle problem - given a station whose
c  colatitude is t, and east longitude rlam (in degrees), returns the
c  colatitude b and east longitude rlong of the place which is at a distance
c  of delta degrees and at an azimuth of alp (clockwise from north).
c  the common block stloc holds the cosine and sine of the station colatitude,
c  and its east longitude.
      common/stloc/ct,st,rlam,ht
      data dtr/.01745329251/
      ca = cos(alp*dtr)
      sa = sin(alp*dtr)
      cd = cos(del*dtr)
      sd = sin(del*dtr)
      cb = cd*ct + sd*st*ca
      sb = sqrt(1.-cb*cb)
      b = acos(cb)/dtr
      if(sb.le.1.e-3) then
c  special case - the point is at the poles
	 rlong = 0
	 return
      endif
      sg = sd*sa/sb
      cg = (st*cd-sd*ct*ca)/sb
      g = atan2(sg,cg)/dtr
      rlong = rlam + g
      return
      end
