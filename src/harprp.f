c
c  Original date of RCS archived version is  Jul 23  1999
c
c
c     $Id: harprp.f,v 1.1 2011/09/23 22:39:01 agnew Exp agnew $
c
c     $Log: harprp.f,v $
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
c  Program harprp
c
c  reads from standard input, a load file (or set of them),
c  selects the numbers associated with a particular quantity,
c  and writes out a file of harmonic constants suitable for tidal
c  prediction by hartid.
c
c  The designators (on the command line) are
c    o - ocean-tide height
c    g - gravity
c    p - potential height
c    z - vertical displacement
c    d - horizontal displacement at azimuth nnn
c    l - extension at azimuth nnn
c    v - volume strain (for poisson's ratio 0.25)
c    s - shear strain at azimuth nnn
c    t - tilt at azimuth nnn
c  where nnn is also given on the command line
c
      character*80 stnam
      character*50 mdnam,dumm
      character*4 dsym
      character*1 modo,type
      complex cload1,cload2,scale,phasor,cz
      dimension cload1(10),cload2(10)
      common /modpar/mdfile,dsym,icte(6),mdnam
c the following pair of functions convert between complex stored as
c amp and phase and complex stored as real and imaginary
      phasor(scale) =
     1 cmplx(cabs(scale),57.2958*atan2(aimag(scale),real(scale)))
      data luo/6/,lui/5/
      data scale/(1.,0.)/azim/0./,ifrst/0/cz/(0.,0.)/
      data cload1/10*(0.,0.)/
      if(iargc().lt.1.or.iargc().gt.2) then
	 write(luo,100)
 100  format('Usage: harprp o | g | p | z | v',/,
     1  'or:  harprp d | l | | s | t  azimuth ')
	 stop
      endif
      call getarg(1,type)
      if(type.eq.'d'.or.type.eq.'l'.or.type.eq.'s'.or.type.eq.'t') then
         if(iargc().ne.2) then
	    write(luo,100)
	    stop
         endif
         call getarg(2,dumm)
         read(dumm,102) azim
 102     format(f10.0)
      endif
c
c now read load file, getting coordinates and other information from
c the first two lines (plus the last "L" line for the phase)
c
 3    read(lui,201,end=33) stnam, rlat, rlam, ht
 201  format(3x,a40,3x,f10.4,f10.4,f10.0)
      read(lui,203) dsym,(icte(j),j=1,6),mdnam
 203  format(2x,a4,5x,6i2,3x,a1,4x,a50)
 5    read(lui,204) dumm
 204  format(a50)
      if(dumm(1:1).eq.'L') modo = dumm(3:3)
      if(dumm(1:1).ne.'X') go to 5
c      gravity or ocean amplitude
      read(lui,205) cload1(1)
 205  format(10x,2f10.4)
      if(type.ne.'o') then
c      potential
         read(lui,205) cload1(10)
c      displacement
         read(lui,207) (cload1(i),i=2,4)
 207     format(10x,6f10.4)
c      tilt
         read(lui,209) (cload1(i),i=5,6)
 209     format(10x,4f10.4)
c      strain
         read(lui,207) (cload1(i),i=7,9)
      endif
      call argand(cload1,10)
c
c adjust to a new azimuth, rotated azim deg (clockwise in the command line,
c converted to counterclockwise below to match the usual formulae) from the
C EW-NS one that the loads are orignally in.
c
      if(azim.ne.0) then
	 caz = cos(-azim/57.2958)
	 saz = sin(-azim/57.2958)
	 cload2(1) = cload1(1)
	 cload2(10) = cload1(10)
	 cload2(4) = cload1(4)
c           horizontal displacements (2 EW 3 NS 4 vert)
	 cload2(2) = caz*cload1(2) + saz*cload1(3)
	 cload2(3) = -saz*cload1(2) + caz*cload1(3)
c           tilts (5 EW 6 NS)
	 cload2(5) = caz*cload1(5) + saz*cload1(6)
	 cload2(6) = -saz*cload1(5) + caz*cload1(6)
c           strains (originally, 7 is EW 8 NS 9 shear (EN) )
	 cload2(7) = caz*caz*cload1(7) + saz*saz*cload1(8)
     1       + 2*caz*saz*cload1(9)
	 cload2(8) = caz*caz*cload1(8) + saz*saz*cload1(7)
     1       - 2*caz*saz*cload1(9)
	 cload2(9) = cload1(9)*(caz*caz-saz*saz)
     1       - caz*saz*(cload1(7)-cload1(8))
	 do 13 i=1,10
 13      cload1(i) = cload2(i)
      endif
c     special conversion for volume strain (do not need EW
c     after here)
      cload1(7) = 0.6667*(cload1(7) + cload1(8))
c
c convert back to amp and phase and write out the results
c
      do 15 i=1,10
      cload1(i) = phasor(cload1(i))
 15   continue
      if(ifrst.eq.0) then
	 ifrst = 1
	 write(luo,213) modo
 213     format(a1)
	 if(modo.eq.'l') write(luo,215) rlam
 215     format(f8.3)
      endif
      if(type.eq.'o') write(luo,217) icte,cload1(1)
      if(type.eq.'g') write(luo,217) icte,cload1(1)
      if(type.eq.'p') write(luo,217) icte,cload1(10)
      if(type.eq.'z') write(luo,217) icte,cload1(4)
      if(type.eq.'d') write(luo,217) icte,cload1(3)
      if(type.eq.'t') write(luo,217) icte,cload1(6)
      if(type.eq.'l') write(luo,217) icte,cload1(8)
      if(type.eq.'s') write(luo,217) icte,cload1(9)
      if(type.eq.'v') write(luo,217) icte,cload1(7)
 217  format(6i2,f9.4,f9.2)
      go to 3
 33   icte(1) = -1
      do 35 i=2,6
 35   icte(i) =  0
      write(luo,217) icte,cz
      stop
      end
      subroutine argand(comp,n)
      complex comp
      dimension comp(n)
      do 3 i=1,n
      p = aimag(comp(i))/57.2958
      a = real(comp(i))
      comp(i) = cmplx(a*cos(p),a*sin(p))
 3    continue
      return
      end
