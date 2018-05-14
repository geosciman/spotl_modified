c
c  Original date of RCS archived version is  Jul  2  2006
c
c
c     $Id: modcon.f,v 1.1 2011/09/23 22:39:01 agnew Exp agnew $
c
c     $Log: modcon.f,v $
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
c program modcon - converts ASCII files of ocean models to binary (Fortran
c sequential unformatted), or the reverse. The ASCII is read in (or written
c out) via standard input or output; the binary file is specified on the 
c command line
c
c   calls no other routines
c
      character*80 mdfile
      character*50 mdnam
      character*4 dsym
      character*1 cdir
      integer*2 ir1,im1
      dimension ir1(5000000),im1(5000000)
      common /modpar/mdfile,dsym,icte(6),mdnam
      common /modlim/tlat,blat,wlong,elong,latc,longc
      if(iargc().ne.2) then
	 write(6,100)
 100     format('usage: modcon [f|r] binaryfile',/
     1     'f for ASCII -> binary, r the reverse')
	 stop
      endif
      call getarg(1,cdir)
      if(cdir.ne.'r'.and.cdir.ne.'f') then
         write(6,100)
         stop
      endif
      call getarg(2,mdfile)
      if(cdir.eq.'f') then
         llu =  5
         read(llu,102) dsym
 102     format(a4)
         read(llu,104) (icte(i),i=1,6)
 104     format(6i3)
         read(llu,106) latd,latf
 106     format(2i8)
         read(llu,106) latd2,latf2
         read(llu,106) latd3,latf3
         read(llu,106) latd4,latf4
         read(llu,106) latc,longc
         read(llu,108) mdnam
 108     format(a50)
         read(llu,110) (ir1(i),i=1,latc*longc)
 110     format(10i7)
         read(llu,110) (im1(i),i=1,latc*longc)
         llu =  2
         open(unit=llu,file=mdfile,status='unknown',access=
     $       'sequential',form="unformatted")
         write(llu) dsym
         write(llu) (icte(i),i=1,6)
         write(llu) latd,latf
         write(llu) latd2,latf2
         write(llu) latd3,latf3
         write(llu) latd4,latf4
         write(llu) latc,longc
         write(llu) mdnam
         write(llu) (ir1(i),i=1,latc*longc)
         write(llu) (im1(i),i=1,latc*longc)
         close(llu)
      endif
      if(cdir.eq.'r') then
         llu =  2
         open(unit=llu,file=mdfile,status='old',access=
     $       'sequential',form="unformatted")
         read(llu) dsym
         read(llu) (icte(i),i=1,6)
         read(llu) latd,latf
         read(llu) latd2,latf2
         read(llu) latd3,latf3
         read(llu) latd4,latf4
         read(llu) latc,longc
         read(llu) mdnam
         read(llu) (ir1(i),i=1,latc*longc)
         read(llu) (im1(i),i=1,latc*longc)
         close(llu)
         llu =  6
         write(llu,102) dsym
         write(llu,104) (icte(i),i=1,6)
         write(llu,106) latd,latf
         write(llu,106) latd2,latf2
         write(llu,106) latd3,latf3
         write(llu,106) latd4,latf4
         write(llu,106) latc,longc
         write(llu,108) mdnam
         write(llu,110) (ir1(i),i=1,latc*longc)
         write(llu,110) (im1(i),i=1,latc*longc)
      endif
      stop
      end
