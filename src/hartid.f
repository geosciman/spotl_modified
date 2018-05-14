c
c  Original date of RCS archived version is  Jul 30  2009
c
c
c     $Id: hartidnew.f,v 1.2 2011/11/18 20:51:13 agnew Exp agnew $
c
c     $Log: hartidnew.f,v $
c     Revision 1.2  2011/11/18 20:51:13  agnew
c     modified for larger number of harmonics (342)
c
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
c  program hartid
c    finds tides by summing harmonics; the amp and phase of these are
c  specified in subroutine admint
c
c    **calls admint,tdfrph,recurs,spline,shells
c
      parameter (nt=342)
      parameter (ncon=20)
      character*1 ipht
      character*40 dumm
      double precision f,p,dr,pi,scr
      dimension a(nt),f(nt),p(nt),hc(2*nt),scr(3*nt),wf(nt)
      dimension idt(6,ncon),amp(ncon),phase(ncon)
      dimension x(600)
      common/date/it(5)
      data dr/.01745329252d0/,pi/3.1415926535898d0/,irli/1/
      data lui/5/,luo/6/,nin/1/
      if(iargc().lt.7.or.iargc().gt.8) then
	 write(luo,100)
 100  format('Usage: hartid y [d | m d] h m s num samp',/,
     1  '  harmonics file read from standard input',/,
     2  '  results written to standard output')
	 stop
      endif
      call getarg(1,dumm)
      read(dumm,102) it(1)
 102  format(i4)
      if(iargc().eq.7) then
        call getarg(2,dumm)
        read(dumm,102) it(2)
        nb=0
      endif
      if(iargc().eq.8) then
        call getarg(2,dumm)
        read(dumm,102) imonth
        call getarg(3,dumm)
        read(dumm,102) iday
        nb=1
        it(2) = iday + mday(it(1),imonth)
      endif
      call getarg(nb+3,dumm)
      read(dumm,102) it(3)
      call getarg(nb+4,dumm)
      read(dumm,102) it(4)
      call getarg(nb+5,dumm)
      read(dumm,102) it(5)
      call getarg(nb+6,dumm)
      read(dumm,104) irnt
 104  format(i6)
      call getarg(nb+7,dumm)
      read(dumm,106) samp
 106  format(f7.0)
c
c   get harmonics information
c
      read(lui,124) ipht
 124  format(a1)
      if(ipht.eq.'l') read(lui,126) rlong
 126  format(f8.0)
 13   read(lui,128,end=17) (idt(i,nin),i=1,6),amp(nin),phase(nin)
 128  format(6i2,2f9.0)
      if(ipht.eq.'l') phase(nin) = phase(nin) + idt(1,nin)*rlong
      nin=nin+1
c  stop reading on a flag or if arrays are full
      if(idt(1,nin-1).ge.0.and.nin.le.ncon) go to 13
 17   nin=nin-1
      if(idt(1,nin).lt.0) nin=nin-1
c
c  interpolate to larger set of harmonics
c
      ntt=nt
      call admint(amp,idt,phase,a,f,p,nin,ntt)
c  set up for first recursion, and normalize frequencies
      do i=1,ntt
        p(i) = dr*p(i)
        f(i) = samp*pi*f(i)/43200.d0
        wf(i) = f(i)
      enddo
 31   irhi = min(irli+599,irnt)
      np = irhi - irli + 1
c set up harmonic coefficients, compute tide, and write out
      do i=1,ntt
        hc(2*i-1) = a(i)*dcos(p(i))
        hc(2*i)  = -a(i)*dsin(p(i))
      enddo
      call recurs(x,np,hc,ntt,wf,scr)
      write(luo,130) (x(i),i=1,np)
 130  format(f13.5)
      if(irhi.eq.irnt) stop
      irli = irhi + 1
c  reset phases to the start of the new section
      do i=1,ntt
        p(i) = p(i) + np*f(i)
        p(i) = dmod(p(i),2.d0*pi)
      enddo
      go to 31
      end
      function mday(iy,m)
c
c  finds day number of day before start of month m, of year iy, in
c   Gregorian intercalation
c
      leap = 1 - (mod(iy,4)+3)/4
      if(mod(iy,100).eq.0.and.mod(iy,400).ne.0) leap=0
      mday = mod((367*(m-2-12*((m-14)/12)))/12+29,365) + leap*((9+m)/12)
      return
      end
