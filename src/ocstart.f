c
c  Original date of RCS archived version is  Feb 15  1995
c
c
c     $Id: ocstart.f,v 1.1 2011/09/23 22:39:02 agnew Exp agnew $
c
c     $Log: ocstart.f,v $
c     Revision 1.1  2011/09/23 22:39:02  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine ocstart(ir1,im1)
c
c    Reads in an ocean model, from a filename passed in common.  The arrays
c  ir1 and im1 hold the model values (real and imaginary parts); the following
c  are returned also (in another common):
c     tlat   latitude (N) of top edge of model grid
c     blat   latitude ( ) of bottom edge of model grid
c     elong   longitude (E) of "right" edge of model grid
c     wlong   longitude (E) of "left" edge of model grid
c     latc   number of cells EW
c     longc   number of cells NS
c
c  other auxiliary information is passed to other routines in common block
c   modpar, which also contains the file name
c
c   This version reads a file of model values, which should contain:
c     bytes 1-4  the Darwin symbol (a4)
c     bytes 5-28 the corresponding C-T-E integers (6i)
c     bytes 29-30 the latitude of the top edge, in degrees N (integer part)
c     bytes 31-32 the latitude of the top edge, in degrees N
c                 (fract part, in thousandths)
c     bytes 33-34 the latitude of the bottom edge, in degrees N (integer part)
c     bytes 35-36 the latitude of the bottom edge, in degrees N
c                 (fract part, in thousandths)
c     bytes 37-38 the longitude of the E edge, in degrees E (integer part)
c     bytes 39-40 the longitude of the E edge, in degrees E
c                 (fract part, in thousandths)
c     bytes 41-42 the longitude of the W edge, in degrees E (integer part)
c     bytes 43-44 the longitude of the W edge, in degrees E
c                 (fract part, in thousandths)
c     bytes 45-47 the number of NS cells
c     bytes 48-49 the number of EW cells
c     bytes 50-99 the name of the model (alpha)
c
c     bytes 100--> the real values of the tide in each cell, running EW
c           and then NS (ie ordered as a matrix), followed by the
c           imaginary values. The units are mm
c
c
      character*80 mdfile
      character*50 mdnam
      character*4 dsym
      integer*2 ir1,im1
      dimension ir1(1),im1(1)
      common/modpar/mdfile,dsym,icte(6),mdnam
      common/modlim/tlat,blat,wlong,elong,latc,longc
      llu =  1
      open(unit=llu,file=mdfile,status='old',access='sequential',
     $    form="unformatted")
      read(llu) dsym
      read(llu) (icte(i),i=1,6)
      read(llu) latd,latf
      tlat = latd + latf/1000.
      read(llu) latd,latf
      blat = latd + latf/1000.
      read(llu) latd,latf
      elong = latd + latf/1000.
      read(llu) latd,latf
      wlong = latd + latf/1000.
      read(llu) latc,longc
      read(llu) mdnam
      read(llu) (ir1(i),i=1,latc*longc)
      read(llu) (im1(i),i=1,latc*longc)
      close(llu)
      return
      end
