c
c  Original date of RCS archived version is  Oct 29  1996
c
c
c     $Id: getcl.f,v 1.3 2012/03/13 21:17:14 agnew Exp agnew $
c
c     $Log: getcl.f,v $
c     Revision 1.3  2012/03/13 21:17:14  agnew
c      Modified polygon common block
c
c     Revision 1.2  2011/11/27 04:30:55  agnew
c     larger dimension for polygon information
c
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine getcl(stnam,rlat,rlong,ht,mdfile,grfile,modo,distm)
c  gets parameters needed from command line; stops program if there
c  are too few
      logical*1 ispoly,use
      character*1 modo,alpoly
      character*80 stnam,mdfile,grfile,dumm
      character*80 polyf,polnam,polynm
c  common block holds polygon file information
      common/polinf/xpoly(500,10),ypoly(500,10),np,npoly(10),
     1              polyf,polynm(10),use(10),ispoly
      common/polmor/polnam,alpoly
      data ispoly/.false./,alpoly/'N'/
      if(iargc().lt.7) then
	 write(6,100)
 100     format('usage: nloadf stname lat long ht modelfile',
     1' greenfile mode [dist] [polyfile]'/,'e.g.',/,
     2'nloadf PFO 33.609 -116.455 1280 m2.tpxo2 grn.gb g',/,
     3'nloadf PFO 33.609 -116.455 1280 m2.tpxo2 grn.gb m 5',/,
     4'nloadf PFO 33.609 -116.455 1280 m2.tpxo2 grn.gb g nobaja',/,
     6'mode is g or l for loading, m to output cells (needs a max',
     7' distance in deg')
	 stop
      endif
      call getarg(1,stnam)
      call getarg(2,dumm)
      read(dumm,102) rlat
 102  format(f10.0)
      call getarg(3,dumm)
      read(dumm,102) rlong
      call getarg(4,dumm)
      read(dumm,102) ht
      call getarg(5,mdfile)
      call getarg(6,grfile)
      call getarg(7,modo)
      if(modo.ne.'g'.and.modo.ne.'l'.and.modo.ne.'m') then
         write(6,103)
 103     format('Mode must be g, l (Greenwich or local phase), or m')
         stop
      endif
      if(modo.eq.'m') then
	 if(iargc().lt.8) then
	    write(6,104)
 104        format('Distance needed for mode m')
	    stop
         endif
         call getarg(8,dumm)
         read(dumm,102) distm
	 if(iargc().ge.9) then
	    ispoly = .true.
	    alpoly = 'Y'
	    call getarg(9,polyf)
	 endif
	 if(iargc().eq.10) call getarg(10,alpoly)
	 return
      endif
      if(iargc().eq.7) return
      if(iargc().gt.9) then
         write(6,105) iargc()
 105     format('Too many arguments (',i2,')')
	 stop
      endif
      call getarg(8,polyf)
      ispoly = .true.
      alpoly = 'Y'
      if(iargc().eq.9) call getarg(9,alpoly)
      return
      end
