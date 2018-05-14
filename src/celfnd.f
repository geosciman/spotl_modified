c
c  Original date of RCS archived version is  Mar  8  1995
c
c
c     $Id: celfnd.f,v 1.1 2011/09/23 22:39:01 agnew Exp agnew $
c
c     $Log: celfnd.f,v $
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine celfnd(rlong,rlato,x,y,i,j,ind)
c  returns the cell index (as i and j also ind), and
c  fractional position within a cell, for a given long
c  and lat (in degrees)
c
      common/modlim/tlat,blat,wlong,elong,latc,longc
c  check to see that we are within model
      if(rlato.gt.tlat.or.rlato.lt.blat) then
	 ind = 0
	 return
      endif
      dlong = amod(rlong-wlong,360.)
      if(dlong.lt.0) dlong = dlong + 360
      if(dlong.gt.elong-wlong) then
	 ind = 0
	 return
      endif
c  put rlong into the range of longitudes that are specified.
      x = (longc*dlong)/(elong-wlong)
      y = (latc*(tlat-rlato))/(tlat-blat)
      i= int(x) + 1
      j= int(y)
      x = x - i + .5
      y = j - y + .5
      ind = i + longc*j
      return
      end
