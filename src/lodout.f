c
c  Original date of RCS archived version is  Jul  5  2005
c
c
c     $Id: lodout.f,v 1.7 2013/03/11 02:21:04 agnew Exp agnew $
c
c     $Log: lodout.f,v $
c     Revision 1.7  2013/03/11 02:21:04  agnew
c     updated version number to 3.3.0.2
c
c     Revision 1.6  2012/06/25 17:14:29  agnew
c     updated version number to 3.3.0.1
c
c     Revision 1.5  2012/03/13 21:19:15  agnew
c      Modified polygon common block and output polygon names
c
c     Revision 1.4  2012/03/06 03:36:28  agnew
c     added lines for new grid options
c
c     Revision 1.3  2011/12/30 23:29:19  agnew
c     changed version number given in output, to 3.3.0
c
c     Revision 1.2  2011/12/30 23:03:36  agnew
c     removed unneeded statement numbers
c
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine lodout(cload,stnam,rlat,rlam,ht,modo)
c  puts the data into a single ASCII rescaled file for all
c  constituents, in order of frequency. Loads are written out in
c  amp and local or Greenwich phase, lags neg
c
      logical*1 ispoly,use
      character*80 stnam,mdfile
      character*80 grname
      character*80 polyf,polnam,polynm
      character*50 mdnam
      character*24 dates
      character*4 dsym
      character*1 modo,statgr,alpoly
      complex cload,toloc,phasor,cmp,dumm
      dimension cload(10),cmp(3)
      common/modpar/mdfile,dsym,icte(6),mdnam
      common/grinf/nring,rin(10),rout(10),rsz(10),statgr(10),grname
      common/polinf/xpoly(500,10),ypoly(500,10),np,npoly(10),
     1              polyf,polynm(10),use(10),ispoly
      common/polmor/polnam,alpoly
      common/coninf/dist,close,clat,clong,numnf,clnf,farnf
      data dr/.0174532935/
      data luo/6/
      write(luo,201) 'S',stnam, rlat, rlam, ht
 201  format(a1,2x,a40,3x,f10.4,f10.4,f10.0)
      write(luo,203) 'O',dsym,(icte(j),j=1,6),mdnam
 203  format(a1,1x,a4,5x,6i2,5x,a50)
      write(luo,205) 'G',grname
 205  format(a1,2x,a78)
      do i=1,nring
        if(statgr(i).eq.'F') write(luo,211) rin(i),rout(i),rsz(i)
 211    format('G    Rings from ',f6.2,' to ',f6.2,' spacing ',
     1   f4.2,' - detailed grid (ocean), seawater')
        if(statgr(i).eq.'L') write(luo,213) rin(i),rout(i),rsz(i)
 213    format('G    Rings from ',f6.2,' to ',f6.2,' spacing ',
     1   f4.2,' - detailed grid (land), freshwater')
        if(statgr(i).eq.'C') write(luo,215) rin(i),rout(i),rsz(i)
 215    format('G    Rings from ',f6.2,' to ',f6.2,' spacing ',
     1   f4.2,' - model grid, seawater')
        if(statgr(i).eq.'G') write(luo,217) rin(i),rout(i),rsz(i)
 217    format('G    Rings from ',f6.2,' to ',f6.2,' spacing ',
     1   f4.2,' - model grid, freshwater')
      enddo
      if(alpoly.ne.'N') then
         write(luo,221) polnam
 221     format('P',1x,a78)
         do i=1,np
            if(use(i).eqv..true.) write(6,223) polynm(i)
 223        format('P   Included polygon: ',a50)
            if(use(i).eqv..false.) write(6,225) polynm(i)
 225        format('P   Excluded polygon: ',a50)
         enddo
      endif
c following three lines depend on a specific routine, which returns time
c information in string   dates
      call fdate(dates)
      write(luo,231) dates
 231  format('C   Version 3.3.0.2 of load program, run at ',a24)
      write(luo,233) close,clat,clong
 233  format('C   closest nonzero load was ',f6.2,' degrees away, at',
     1 2f8.2)
      if(numnf.ne.0) write(luo,235) numnf,clnf,farnf
 235  format('C   ',i6,' zero loads found where ocean present, range',
     1 f7.2,'-',f7.2,' deg')
      if(modo.eq.'g') write(luo,241)
 241  format('L g          Phases are Greenwich, lags negative')
      if(modo.eq.'l') write(luo,243)
 243  format('L l          Phases are local, lags negative')
      write(luo,249)
 249  format('X')
      isp = icte(1)
      toloc = cmplx(cos(dr*isp*rlam),sin(dr*isp*rlam))
      if(modo.eq.'g') toloc = cmplx(1.,0.)
c
c      gravity
c
      cmp(1) = 1.e8*toloc*cload(1)
      dumm = conjg(cmp(1))
      cmp(1) = phasor(dumm)
      write(luo,251) cmp(1)
 251  format('g',9x,2f10.4)
c
c      potential height
c
      cmp(1) = 1.e3*toloc*cload(10)
      dumm = conjg(cmp(1))
      cmp(1) = phasor(dumm)
      write(luo,253) cmp(1)
 253  format('p',9x,2f10.4)
c
c      displacement (output as E N Z)
c
c  change sign of N and E (see CJ II-143, 5 Jan 82)
      cmp(1) = -1.e3*toloc*cload(4)
      cmp(2) = -1.e3*toloc*cload(3)
      cmp(3) =  1.e3*toloc*cload(2)
      do i = 1,3
        dumm = conjg(cmp(i))
        cmp(i) = phasor(dumm)
      enddo
      write(luo,255) cmp
 255  format('d',9x,6f10.4)
c
c      tilt (output as E N)
c
      cmp(1) = 1.e9*toloc*cload(6)
      cmp(2) = 1.e9*toloc*cload(5)
      do i = 1,2
        dumm = conjg(cmp(i))
        cmp(i) = phasor(dumm)
      enddo
      write(luo,257) cmp(1),cmp(2)
 257  format('t',9x,4f10.4)
c
c      strain (output as E, N, EN shear)
c
      cmp(1) = 1.e9*toloc*cload(8)
      cmp(2) = 1.e9*toloc*cload(7)
      cmp(3) = 1.e9*toloc*cload(9)
      do i = 1,3
        dumm = conjg(cmp(i))
        cmp(i) = phasor(dumm)
      enddo
      write(luo,259) cmp
 259  format('s',9x,6f10.4)
      return
      end
      function phasor(dumm)
c
c  returns, as a complex number, the amplitude and phase of the complex
c   quantity input
c
      complex phasor,dumm,cz
      data cz/(0.,0.)/
      if(dumm.eq.cz) phasor=cz
      if(dumm.ne.cz) phasor=
     1 cmplx(cabs(dumm),57.2958*atan2(aimag(dumm),real(dumm)))
      return
      end
