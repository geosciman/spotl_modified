c
c  Original date of RCS archived version is  Nov 24  2007
c
c
c     $Id: oclook.f,v 1.7 2013/03/11 02:01:28 agnew Exp agnew $
c
c     $Log: oclook.f,v $
c     Revision 1.7  2013/03/11 02:01:28  agnew
c     Fixed divide-by-zero if amp (so possibly rho) is zero
c
c     Revision 1.6  2012/11/20 16:29:32  agnew
c     fixed bug with coarse flag
c
c     Revision 1.5  2012/06/25 01:30:11  agnew
c     corrected use of density output from ocmodl by adding
c     density common block
c
c     Revision 1.4  2012/03/13 21:17:50  agnew
c      Modified polygon common block
c
c     Revision 1.3  2011/12/23 00:54:45  agnew
c     fixed index error for g a p r i option, modernized loops, cleaned up code
c     for e option
c
c     Revision 1.2  2011/11/27 04:31:33  agnew
c     larger dimension for polygon information
c
c     Revision 1.1  2011/09/23 22:39:02  agnew
c     Initial revision
c
c
c-------------------------------------
c   program oclook
c  looks at model values: either
c     prints out values for a particular place
c     outputs grid within an area (for plotting)
c     outputs values within an area
c
c
      logical*1 ispoly,use
      character*1 modo,fingrd
      character*4 dsym
      character*80 mdfile,dumm
      character*80 polyf,polynm
      character*50 mdnam
      character*40 stnam
      complex amp,toloc,phasor,oamp,cz
      dimension i1(2),i2(2)
c  common block holds polygon file information
      common/polinf/xpoly(500,10),ypoly(500,10),np,npoly(10),
     1              polyf,polynm(10),use(10),ispoly
      common /modpar/mdfile,dsym,icte(6),mdnam
      common /modlim/tlat,blat,wlong,elong,latc,longc
      common/ldens/rho
c
      phasor(toloc) =
     1 cmplx(cabs(toloc),57.2958*atan2(aimag(toloc),real(toloc)))
      data ispoly/.false./
      data cz/(0.,0.)/,oamp/(0.,0.)/
      data fingrd/'C'/
      data stnam/'  Interpolated ocean-tide value         '/
      if(iargc().lt.2.or.iargc().gt.7) then
	 write(6,100)
 100     format('usage: oclook modelfile e',/,
     1    ' or    oclook modelfile lat long d',/,
     2    ' or    oclook modelfile lat long o',/,
     3    ' or    oclook modelfile slat nlat wlong elong g',/,
     4    ' or    oclook modelfile slat nlat wlong elong a ',/,
     5    ' or    oclook modelfile slat nlat wlong elong p ',/,
     6    ' or    oclook modelfile slat nlat wlong elong r ',/,
     7    ' or    oclook modelfile slat nlat wlong elong i ',/,
     8    ' (in all cases but the first, a polyfile may be put ',/,
     9    '      on the end)')
	 stop
      endif
      call getarg(1,mdfile)
      if(iargc().eq.2) then
         call getarg(2,modo)
      endif
      if(iargc().eq.4.or.iargc().eq.5) then
         call getarg(2,dumm)
         read(dumm,102) rlato
 102     format(f10.0)
         call getarg(3,dumm)
         read(dumm,102) rlong
         call getarg(4,modo)
         if(iargc().eq.5) then
            ispoly = .true.
            call getarg(5,polyf)
         endif
      endif
      if(iargc().ge.6) then
         call getarg(2,dumm)
         read(dumm,102) slat
         call getarg(3,dumm)
         read(dumm,102) rnlat
         call getarg(4,dumm)
         read(dumm,102) ftlong
         call getarg(5,dumm)
         read(dumm,102) rilong
         call getarg(6,modo)
         if(iargc().eq.7) then
            ispoly = .true.
            call getarg(7,polyf)
         endif
      endif
c
c  done with getting arguments
c
      if(modo.eq.'d') then
	 fingrd = 'C'
         call ocmodl(rlato,rlong,fingrd,amp)
         if(amp.ne.cz) amp=amp/rho
         call celfnd(rlong,rlato,x,y,i,j,ind)
	 write(6,104) i,j,ind,x,y
 104     format('Cell indices are ',i5,i5,' (',i7,')',
     1       ' with fractional locations ',2f6.3)
	 if(amp.ne.cz) write(6,106) amp,phasor(amp)
	 if(amp.eq.cz) write(6,107)
 106     format('Cell amps (R&I,[amp&ph]) are ',2f10.4,/,' [',
     1          2f10.4,' ]')
 107     format('No ocean in this cell')
	 fingrd = 'F'
         call ocmodl(rlato,rlong,fingrd,amp)
         if(amp.ne.cz) amp=amp/rho
	 if(amp.ne.cz) write(6,108) amp,phasor(amp)
	 if(amp.eq.cz) write(6,109)
 108     format('Fine-grid (interpolated) amps (R&I,[amp&ph]) are ',
     1      2f10.4,/,' [',2f10.4,' ]')
 109     format('Interpolated tide is zero')
	 if(ispoly) then
           call chkgon(rlong,rlato,iok)
	   if(iok.eq.0) write(6,110)
 110       format('Point would be excluded by polygon')
         endif
      elseif(modo.eq.'o') then
         fingrd = 'F'
         call ocmodl(rlato,rlong,fingrd,amp)
         if(amp.ne.cz) amp=amp/rho
	 if(ispoly) then
	    call chkgon(rlong,rlato,iok)
	    if(iok.eq.0) amp = cz
         endif
	 if(amp.eq.cz) then
            write(6,119)
 119        format('Tide is zero or excluded by polygon')
	 endif
	 if(amp.ne.cz) then
            write(6,121) 'S',stnam, rlato, rlong, ht
 121        format(a1,2x,a40,3x,f10.4,f10.4,f10.0)
            write(6,123) 'O',dsym,(icte(j),j=1,6),mdnam
 123        format(a1,1x,a4,5x,6i2,5x,a50)
            write(6,125)
 125        format('L g          Phases are Greenwich, lags negative')
            write(6,127)
 127        format('X')
            amp = conjg(amp)
            amp = phasor(amp)
            write(6,129) amp
 129        format('o',9x,2f10.4)
	 endif
c
c  for options g a p r i, this section outputs either the grid cells
c    or the tidal values in the cells, within the region specified on
c    the command line
c  for option e, this section outputs the edges of the boundary between
c    zero and nonzero cells
c
      else
         call ocmodl(0.,0.,fingrd,amp)
         if(amp.ne.cz) amp=amp/rho
         if(modo.eq.'e') then
            slat=blat
            rnlat=tlat
            ftlong=wlong
            rilong=elong  
         endif
	 if(rilong.lt.ftlong) rilong=rilong+360
         call fndgde(ftlong,rilong,slat,rnlat,il,ih,jl,jh)
         if(il.eq.0.and.ih.eq.0) then
            write(6,132)
 132        format('Area is outside model')
	    stop
         endif
c following code takes care of case where longitudes specified span
c edge of model
	 if(ih.lt.il) then
	    ngrs = 2
	    i1(1) = 1
	    i2(1) = ih
	    i1(2) = il
	    i2(2) = longc
	 elseif(ih.eq.il) then
	    ngrs = 1
	    i1(1) = 1
	    i2(1) = longc
	 else
	    ngrs = 1
	    i1(1) = il
	    i2(1) = ih
	 endif
	 dlat = (tlat-blat)/latc
	 dlong = (elong-wlong)/longc
c
c now scan through cells by columns, then rows; for option e output
c   any EW-running edges between zero and nonzero, or at top and
c   bottom
c
	 do k=1,ngrs
	   do i=i1(k),i2(k)
	     do j=jl,jh
	       clat = tlat - (j-1)*dlat - dlat/2
	       clong = wlong + (i-1)*dlong + dlong/2
	       if(clong.gt.180) clong = clong - 360
               call ocmodl(clat,clong,fingrd,amp)
               if(amp.ne.cz) amp=amp/rho
	       if(ispoly) then
                 call chkgon(clong,clat,iok)
	         if(iok.eq.0) amp=cz
               endif
c write out EW edges between zero and nonzero cells, including the case
c   where the zero cell is one outside the grid
	       if(modo.eq.'e') then
                  if(amp.eq.cz.and.oamp.ne.cz) then
                     write(6,142) clong-dlong/2,clat+dlat/2
                     write(6,142) clong+dlong/2,clat+dlat/2
                     write(6,*) '370 100'
                  endif
                  if(amp.ne.cz.and.(j.eq.jl.or.oamp.eq.cz)) then
                     write(6,142) clong-dlong/2,clat+dlat/2
                     write(6,142) clong+dlong/2,clat+dlat/2
                     write(6,*) '370 100'
                  endif
                  if(amp.ne.cz.and.j.eq.jh) then
                     write(6,142) clong-dlong/2,clat+dlat/2
                     write(6,142) clong+dlong/2,clat+dlat/2
                     write(6,*) '370 100'
                  endif
               endif
c
c all other options, for writing out result in, or coordinates of,
c  individual grid cells
c
	       if(amp.ne.cz) then
                  toloc = phasor(amp)
	          if(modo.eq.'r') write(6,140) clong,clat,real(amp)
	          if(modo.eq.'i') write(6,140) clong,clat,aimag(amp)
	          if(modo.eq.'a') write(6,140) clong,clat,real(toloc)
	          if(modo.eq.'p') write(6,140) clong,clat,aimag(toloc)
 140              format(2f10.4,g14.6)
	          if(modo.eq.'g') then
                     write(6,142) clong-dlong/2,clat-dlat/2
                     write(6,142) clong-dlong/2,clat+dlat/2
                     write(6,142) clong+dlong/2,clat+dlat/2
                     write(6,142) clong+dlong/2,clat-dlat/2
                     write(6,142) clong-dlong/2,clat-dlat/2
                     write(6,*) '370 100'
 142                 format(2f10.4)
                  endif
	       endif
	       oamp = amp
             enddo
           enddo
         enddo
c
c for the e option, do a second scan, by rows, to find all the NS-running
c   edges
	 if(modo.eq.'e') then
            oamp=cz
	    do k=1,ngrs
	      do j=jl,jh
	        do i=i1(k),i2(k)
	           clat = tlat - (j-1)*dlat - dlat/2
	           clong = wlong + (i-1)*dlong + dlong/2
	           if(clong.gt.180) clong = clong - 360
                   call ocmodl(clat,clong,fingrd,amp)
                   if(amp.eq.cz.and.oamp.ne.cz) then
                      write(6,142) clong-dlong/2,clat+dlat/2
                      write(6,142) clong-dlong/2,clat-dlat/2
                      write(6,*) '370 100'
                   endif
                   if(amp.ne.cz.and.(i.eq.i1(k).or.oamp.eq.cz)) then
                      write(6,142) clong-dlong/2,clat+dlat/2
                      write(6,142) clong-dlong/2,clat-dlat/2
                      write(6,*) '370 100'
                   endif
                   if(amp.ne.cz.and.i.eq.i2(k)) then
                      write(6,142) clong+dlong/2,clat+dlat/2
                      write(6,142) clong+dlong/2,clat-dlat/2
                      write(6,*) '370 100'
                   endif
	           oamp = amp
                enddo
              enddo
            enddo
	 endif
      endif
      stop
      end
