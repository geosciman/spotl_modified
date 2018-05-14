c
c  Original date of RCS archived version is  May 19  1998
c
c
c     $Id: loadcomb.f,v 1.2 2011/12/30 22:59:49 agnew Exp agnew $
c
c     $Log: loadcomb.f,v $
c     Revision 1.2  2011/12/30 22:59:49  agnew
c     fixed error in format statement 211
c
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
      program loadcomb
c
c  combines load files (eg global and local), also with body tides
c
c  reads from standard input at least one load file (two if they
c  are to be combined).
c
c  If there are two files, the values are added, after scaling the
c  values in the second load file by an amplitude and phase (if these
c  are given on the command line). The summed values are written out,
c  with the tilts and strains being (optionally) adjusted to a new direction,
c  if the rotation to this is given on the command line. A single file may
c  be scaled and rotated; added to the body tide; or used to produce the
c  body tide.
c
      character*80 stnam,hedlin
      character*50 mdnam,dumm,dumm2
      character*4 dsym
      character*1 body,dumc
      complex cload1,cload2,scale
      dimension cload1(10),cload2(10),btlt(2),bstr(3),bdis(3)
      common /modpar/mdfile,dsym,icte(6),mdnam
c the following pair of functions convert between complex stored as
c amp and phase and complex stored as real and imaginary
      data luo/6/,lui/5/
      data scale/(1.,0.)/rotdir/0./
      data dumm2/'Body Tide computed for spherical Earth            '/
      data cload1/10*(0.,0.)/
      if(iargc().lt.1.or.iargc().gt.4) then
	 write(luo,100)
 100  format('Usage: loadcomb c [rotation] [scale(amp) scale(phase)] ',/,
     1       '        (reads two load files from std input)',/,
     2  '  or:  loadcomb b [rotation]',/,
     3  '        (reads load file to get location and harmonic, outputs
     4 body tide only)',/,
     5  '  or:  loadcomb t [rotation] [scale (amp) scale (ph)] ',/,
     4  '        (adds body tide to load read in)',/,
     5  '  or:  loadcomb r [rotation] [scale (amp) scale (ph)] ',/,
     6  '        (scales and rotates single load tide read in)')
	 stop
      endif
      call getarg(1,body)
      if(body.ne.'b'.and.body.ne.'t'.and.body.ne.'r'.
     1 and.body.ne.'c'.or.(body.eq.'b'.and.iargc().gt.2)) then
	 write(luo,100)
         stop
      endif
      if(body.eq.'b'.and.iargc().eq.2) then
         call getarg(2,dumm)
         read(dumm,102) rotdir
 102     format(f10.0)
      endif
      if(body.eq.'c'.or.body.eq.'r'.or.body.eq.'t') then
         if(iargc().eq.2) then
            call getarg(2,dumm)
            read(dumm,102) rotdir
         endif
         if(iargc().eq.3) then
            call getarg(2,dumm)
            read(dumm,102) amp
            call getarg(3,dumm)
            read(dumm,102) phase
	    scale = cmplx(amp,phase)
	    call argand(scale,1)
         endif
         if(iargc().eq.4) then
            call getarg(4,dumm)
            read(dumm,102) rotdir
         endif
      endif
c
c done with reading command line
c
      if(body.eq.'b'.or.body.eq.'t') then
c  read in needed information about place from a load file, to compute
c  body tides
         read(lui,201) dumc,stnam, rlat, rlam, ht
 201     format(a1,2x,a40,3x,f10.4,f10.4,f10.0)
         read(lui,203) dumc,dsym,(icte(j),j=1,6),mdnam
 203     format(a1,1x,a4,5x,6i2,5x,a50)
c  compute the body tides, with sign adjusted to make the potential positive
	 call sph(rlat,rlam,ht)
         call cteput(icte,amp,n,m)
	 amp=abs(amp)
         call etharm(amp,n,m,pot,grav,bdis,btlt,bstr)
	 ips = sign(1.0,pot)
	 cload1(1) = ips*cmplx(grav,0.)
	 cload1(10) = ips*cmplx(pot,0.)
c  displacements rearranged to EW NS Z
	 cload1(2) = ips*cmplx(bdis(2),0.)
	 cload1(3) = ips*cmplx(0.,bdis(3))
	 cload1(4) = ips*cmplx(bdis(1),0.)
c  tilts rearranged to EW NS
	 cload1(5) = ips*cmplx(0.,btlt(2))
	 cload1(6) = ips*cmplx(btlt(1),0.)
c  strains rearranged to EW NS EN (shear, which needs a change of sign)
	 cload1(7) = ips*cmplx(bstr(2),0.)
	 cload1(8) = ips*cmplx(bstr(1),0.)
	 cload1(9) = -ips*cmplx(0.,bstr(3))
         write(luo,201) 'S',stnam, rlat, rlam, ht
         if(rotdir.ne.0.0) write(luo,211) rotdir
 211     format('R   vector and tensor quantities in frame rotated by',
     1          f8.2)
         write(luo,203) dumc,dsym,(icte(j),j=1,6),dumm2
	 if(body.eq.'b') then
            write(luo,212)
 212        format('L l          Phases are local, lags negative'/'X')
c body tides only: rotate if need be, write out, and stop
            if(rotdir.ne.0) call rotlod(cload1,rotdir)
	    call phasor(cload1,10)
            write(luo,225) cload1(1)
            write(luo,226) cload1(10)
            write(luo,227) (cload1(i),i=2,4)
            write(luo,229) (cload1(i),i=5,6)
            write(luo,231) (cload1(i),i=7,9)
 225        format('g',9x,2f10.4)
 226        format('p',9x,2f10.4)
 227        format('d',9x,6f10.4)
 229        format('t',9x,4f10.4)
 231        format('s',9x,6f10.4)
	    stop
         endif
	 if(body.eq.'t') then
            write(luo,203) dumc,dsym,(icte(j),j=1,6),mdnam
c  skip to reading loads into array cload2, since this is the one
c  that is scaled
	    go to  3
         endif
      endif
      if(body.eq.'r') then
c  only reading in one load; read a few lines, and then
c  skip to reading loads into array cload2, since this is the one
c  that is scaled
         read(lui,201) dumc,stnam, rlat, rlam, ht
         write(luo,201) 'S',stnam, rlat, rlam, ht
         read(lui,203) dumc,dsym,(icte(j),j=1,6),mdnam
         write(luo,203) dumc,dsym,(icte(j),j=1,6),mdnam
         if(rotdir.ne.0.0) write(luo,211) rotdir
	 go to 3
      endif
c
c in this case we are combining two load files; read and write the
c header lines for the first one
c
      read(lui,201) dumc,stnam, rlat, rlam, ht
      write(luo,201) dumc,stnam, rlat, rlam, ht
      if(rotdir.ne.0.0) write(luo,211) rotdir
      read(lui,203) dumc,dsym,(icte(j),j=1,6),mdnam
      write(luo,203) dumc,dsym,(icte(j),j=1,6),mdnam
 1    read(lui,200) hedlin
 200  format(a)
      if(hedlin(1:1).ne.'X') write(luo,200) hedlin
      if(hedlin(1:1).ne.'X') go to 1
c      gravity
      read(lui,205) cload1(1)
 205  format(10x,2f10.4)
c      potential
      read(lui,205) cload1(10)
c      displacement
      read(lui,207) (cload1(i),i=2,4)
 207  format(10x,6f10.4)
c      tilt
      read(lui,209) (cload1(i),i=5,6)
 209  format(10x,4f10.4)
c      strain
      read(lui,207) (cload1(i),i=7,9)
      call argand(cload1,10)
c
c read second load file (the names etc will be from the first one, so we
c have to read and ignore the first line) and add them
c
      read(lui,200) dumm
 3    if(body.eq.'r'.and.scale.ne.(1.,0.)) write(luo,213) amp,phase
 213     format('R   scaling of',g13.5,' (amp) and ',f9.2,' (phase)',
     1      ' applied')
      if(body.eq.'t'.and.scale.ne.(1.,0.)) write(luo,215) amp,phase
 215  format('R   scaling of',g13.5,' (amp) and ',f9.2,' (phase)',
     1      ' applied to load')
      if(body.eq.'c'.and.scale.ne.(1.,0.)) write(luo,217) amp,phase
 217  format('R   scaling of',g13.5,' (amp) and ',f9.2,' (phase)',
     1      ' applied to second load')
 5    read(lui,200) hedlin
      write(luo,200) hedlin
      if(hedlin(1:1).ne.'X') go to 5
      read(lui,205) cload2(1)
      read(lui,205) cload2(10)
      read(lui,207) (cload2(i),i=2,4)
      read(lui,209) (cload2(i),i=5,6)
      read(lui,207) (cload2(i),i=7,9)
      call argand(cload2,10)
      do 7 i=1,10
      cload1(i) = cload1(i) + scale*cload2(i)
 7    continue
c
c adjust to a new direction, rotated rotdir deg (clockwise in the command line,
c converted to counterclockwise below to match the usual formulae) from the
C EW-NS one that the loads are orignally in.
c
      if(rotdir.ne.0) call rotlod(cload1,rotdir)
c
c convert back to amp and phase and write out the results
c
      call phasor(cload1,10)
      write(luo,225) cload1(1)
      write(luo,226) cload1(10)
      write(luo,227) (cload1(i),i=2,4)
      write(luo,229) (cload1(i),i=5,6)
      write(luo,231) (cload1(i),i=7,9)
      stop
      end
c
c
      subroutine rotlod(cl,az)
c rotates loads by amount az
      complex cl,cl2
      dimension cl(10),cl2(10)
      caz = cos(-az/57.2958)
      saz = sin(-az/57.2958)
      cl2(1) = cl(1)
      cl2(10) = cl(10)
      cl2(4) = cl(4)
c           horizontal displacements (2 EW 3 NS 4 vert)
      cl2(2) = caz*cl(2) + saz*cl(3)
      cl2(3) = -saz*cl(2) + caz*cl(3)
c           tilts (5 EW 6 NS)
      cl2(5) = caz*cl(5) + saz*cl(6)
      cl2(6) = -saz*cl(5) + caz*cl(6)
c           strains (originally, 7 is EW 8 NS 9 shear (EN) )
      cl2(7) = caz*caz*cl(7) + saz*saz*cl(8) + 2*caz*saz*cl(9)
      cl2(8) = caz*caz*cl(8) + saz*saz*cl(7) - 2*caz*saz*cl(9)
      cl2(9) = cl(9)*(caz*caz-saz*saz) - caz*saz*(cl(7)-cl(8))
      do 13 i=1,10
 13   cl(i) = cl2(i)
      return
      end
c
c
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
c
c
      subroutine phasor(comp,n)
      complex comp
      dimension comp(n)
      do 3 i=1,n
      p = 57.2958*atan2(aimag(comp(i)),real(comp(i)))
      a = cabs(comp(i))
      comp(i) = cmplx(a,p)
 3    continue
      return
      end
