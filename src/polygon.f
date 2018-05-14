c
c  Original date of RCS archived version is  Mar 25  1996
c
c
c     $Id: polygon.f,v 1.3 2012/03/13 21:18:10 agnew Exp agnew $
c
c     $Log: polygon.f,v $
c     Revision 1.3  2012/03/13 21:18:10  agnew
c      Modified polygon common block to include polygon names
c
c     Revision 1.2  2011/11/27 04:32:37  agnew
c     larger dimension for polygon information
c
c     Revision 1.1  2011/09/23 22:39:02  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine plstart
c
c    Reads in, if requested, information about
c  polygonal areas (up to 10) within which the point either must or
c  must not fall.
c
      logical*1 ispoly,use
      character*1 iuse,alpoly
      character*80 polyf,polnam,polynm
c  common block holds polygon file information
      common/polinf/xpoly(500,10),ypoly(500,10),np,npoly(10),
     1              polyf,polynm(10),use(10),ispoly
      common/polmor/polnam,alpoly
      np=0
      llu = 4
      open(llu,file=polyf,status='old')
c first read name for whole file, and number of polygons; for each one,
c then read name (not used), number of points, + for point must be in,
c - for may not be, and then the points (long and lat) themselves.
c   The polygon coordinates may not cross a 360-degree discontinuity
c and must be in the range -180 to 180 or 0-360 (E +)
      read(llu,102) polnam
 102  format(a)
      read(llu,*) np
      do i=1,np
        read(llu,102) polynm(i)
        read(llu,*) npoly(i)
        read(llu,104) iuse
 104    format(a1)
        if(iuse.eq.'+') use(i) = .true.
        if(iuse.eq.'-') use(i) = .false.
c global override from variable alpoly, read from command line
        if(alpoly.eq.'+') use(i) = .true.
        if(alpoly.eq.'-') use(i) = .false.
        do j=1,npoly(i)
          read(llu,*) xpoly(j,i),ypoly(j,i)
        enddo
      enddo
      close(llu)
      return
      end
      subroutine chkgon(rlong,rlat,iok)
c
c  returns iok=0 if
c     1. there is any polygon (of all those read in) in which the
c        coordinate should not fall, and it does
c             or
c     2. the coordinate should fall in at least one polygon
c        (of those read in) and it does not
c
c    otherwise returns iok=1
c
c    $$$$calls inpoly
c
      logical*1 ispoly,use
      character*80 polyf,polynm
      dimension xu(500),yu(500)
c  common block holds polygon file information
      common/polinf/xpoly(500,10),ypoly(500,10),np,npoly(10),
     1              polyf,polynm(10),use(10),ispoly
c first convert rlong to -180/180 (or leave in this, or 0/360)
      if(rlong.lt.-180) rlong = rlong + 360
c now make rlong2 fall into the other system than rlong (unless rlong is
c between 0 and 180)
      if(rlong.lt.0) rlong2 = rlong + 360
      if(rlong.gt.180) rlong2 = rlong - 360
      do 5 ip=1,np
      if(.not.use(ip)) then
c polygon is one we should not be in; test to see if we are, and if so set
c iok to 0 and return - check all such before seeing if we might also be
c in a polygon we should be
         do 3 i=1,npoly(ip)
            xu(i) = xpoly(i,ip)
            yu(i) = ypoly(i,ip)
 3       continue
	 if(inpoly(rlong,rlat,xu,yu,npoly(ip)).ne.0.or.
     1      inpoly(rlong2,rlat,xu,yu,npoly(ip)).ne.0) then
	    iok=0
	    return
	 endif
      endif
 5    continue
      ianyok=0
      do 15 ip=1,np
      if(use(ip)) then
c polygon is one we should be in; test to see if we are, and if so set
c iok to 1 and return
         ianyok = ianyok+1
         do 13 i=1,npoly(ip)
            xu(i) = xpoly(i,ip)
            yu(i) = ypoly(i,ip)
 13      continue
         if(inpoly(rlong,rlat,xu,yu,npoly(ip)).ne.0.or.
     1      inpoly(rlong2,rlat,xu,yu,npoly(ip)).ne.0) then
	    iok=1
	    return
	 endif
      endif
 15   continue
c not inside any polygons; set iok to 0 if there are any we should have
c been in
      iok = 1
      if(ianyok.gt.0) iok = 0
      return
      end
