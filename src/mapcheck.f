c
c  Original date of RCS archived version is  Aug  3  2001
c
c
c     $Id: mapcheck.f,v 1.1 2011/09/23 22:39:01 agnew Exp agnew $
c
c     $Log: mapcheck.f,v $
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
c   program mapmak
c  reads from the random-access land-sea data file and outputs:
c   1. a printer-plot of the land and sea within 300 km of a specified
c      point
c   2. a printer-plot of a 2 by 2 degree area, 48 by 48 points, to
c      show the local structure of the land-sea file.
c
c $$$$calls lndsea, mrsden, invspt
c
      character*1 lin,opt
      character*50 dumm
      dimension lin(120)
      common/stloc/ct,st,rlam,ht
      data ht/0./,degk/111./,dr/.017453293/,
     1 del/0./,minaz/60/,delmax/12./
      if(iargc().lt.3) then
	 write(6,100)
 100     format('usage: mapcheck r lat long dist',/,
     1          ' or    mapcheck h lat long ',/,
     2    ' or    mapcheck o lat long')
	 stop
      endif
      call getarg(1,opt)
      call getarg(2,dumm)
      read(dumm,102) rlat
 102  format(f10.0)
      call getarg(3,dumm)
      read(dumm,102) rlam
      if(opt.eq.'r') then
         call getarg(4,dumm)
	 read(dumm,102) stp
      endif
      cl = dr*(90.-rlat)
      ct = cos(cl)
      st = sin(cl)
      if(opt.eq.'r') then
         stp = stp/33.
c 43 lines of 67 characters gives a square map
         fac = 33./21.
         do 9 jj=1,43
         y = fac*stp*(22-jj)
         maxp = 0
         do 7 kk=1,67
         lin(kk) = '-'
         x = stp*(kk-34)
         del = sqrt(x**2+y**2)/degk
         asz = 0.
         if(del.ne.0) asz = atan2(y,x)/dr
         asz = 90-asz
         call invspt(asz,del,b,rlong)
         rlat = 90 - b
         call lndsea(rlat,rlong,lnd)
         if(lnd.eq.2) lin(kk)=' '
         if(del.eq.0) lin(kk)='X'
         if(lin(kk).ne.' ') maxp = kk
 7       continue
         if(maxp.eq.0) write(6,110)
 110     format(1x)
         if(maxp.ne.0) write(6,112) (lin(i),i=1,maxp)
 112     format(1x,67a1)
 9       continue
      elseif(opt.eq.'h') then
         lat=nint(rlat)
         long=nint(rlam)
         off = 1./128.
         dt = 1./64.
         do 15 i=1,64
         rlat = lat + off + (64-i)*dt
         do 13 j=1,64
         lin(j) = '-'
         rlong = long + off + (j-1.)*dt
         if(i.eq.1.and.j.eq.1) write(6,122) rlat,rlong
 122     format(' corner coordinates: ',2f10.3)
         call lndsea(rlat,rlong,lnd)
 13      if(lnd.eq.2) lin(j) = '  '
         write(6,112) (lin(k),k=1,64)
         if(i.eq.64) write(6,122) rlat,rlong
 15      continue
      elseif(opt.eq.'o') then
 23      if(del.lt.0.6) del = del + .01
         if(del.ge.0.6.and.del.lt.1.2) del = del + .05
         if(del.ge.1.2) del = del + .2
         naz = 360.*sin(dr*del)
         naz = max0(naz,minaz)
         azstp = 360./naz
         alp = azstp/2.
         do 25 jj=1,naz
         call invspt(alp,del,b,rlong)
         call lndsea(90.-b,rlong,lnd)
         if(lnd.eq.1) go to 25
c  have an ocean cell - give distance and azimuth
         write(6,132) del,alp
 132     format(' Nearest ocean at distance ',f6.2,' deg, azimuth ',
     1           f6.1)
         stop
 25      alp = alp + azstp
         if(del.le.delmax) go to 23
         write(6,134) delmax
 134     format(' found no water within ',f6.1,' degrees')
      endif
      stop
      end
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
      if(sb.le.1.e-3) go to 1
      sg = sd*sa/sb
      cg = (st*cd-sd*ct*ca)/sb
      g = atan2(sg,cg)/dtr
      rlong = rlam + g
      return
c  special case - the point is at the poles
 1    rlong = 0
      return
      end
