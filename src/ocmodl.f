c
c  Original date of RCS archived version is  Jul  2  2006
c
c
c     $Id: ocmodl.f,v 1.6 2012/06/25 01:29:40 agnew Exp agnew $
c
c     $Log: ocmodl.f,v $
c     Revision 1.6  2012/06/25 01:29:40  agnew
c     added density common block, corrected comments
c
c     Revision 1.5  2012/03/13 21:17:32  agnew
c      Modified polygon common block
c
c     Revision 1.4  2012/03/06 03:38:09  agnew
c     added new grid options, density, and call to seadens
c
c     Revision 1.3  2011/11/27 04:31:57  agnew
c     larger dimension for polygon information
c
c     Revision 1.2  2011/11/18 16:52:49  agnew
c     corrected error in interpolation (incorrect indexing of surrounding cells
c     in one case)
c
c     Revision 1.1  2011/09/23 22:39:02  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine ocmodl(rlato,rlong,fingrd,amp)
c
c  For given latitude and east longitude (in degrees and ge 0) returns
c  the complex amplitude of the tide (in kg/m**2, with greenwich phase).
c  this amplitude is of course zero on land.
c    If fingrd is "F", then the land-sea grid used is a detailed one;
c  otherwise it is the one used by the ocean model.
c
c    On the first call, this calls a separate routine (ocstart) which
c  gets auxiliary information (passed to other routines in common) and
c  returns the two arrays of real and imaginary amplitude.  Later calls
c  simply go to the appropriate point of the arrays. 
c
c    The first call also reads in, if requested, information about
c  polygonal areas (up to 10) within which the point either must or
c  must not fall.
c
c   $$$$calls ocstart,lndsea,ibit
c
      save ir1,im1,newf
      character*1 fingrd
      integer*2 ir1,im1
      complex amp,bg,cz,bgs
      integer*4 ind
      logical*1 ispoly,use
      character*80 polyf,polynm
c
c  dimensioning of arrays below must be large enough to handle the largest
c  models
c
      dimension ir1(5000000),im1(5000000)
      dimension bg(4)
      common /modlim/tlat,blat,wlong,elong,latc,longc
c  common block holds polygon file information
      common/polinf/xpoly(500,10),ypoly(500,10),np,npoly(10),
     1              polyf,polynm(10),use(10),ispoly
      common/coninf/dist,close,clat,clong,numnf,clnf,farnf
      common/ldens/rhol
      data cz/(0.,0.)/,newf/1/
      if(newf.eq.1) then
c new files - load arrays (polygon information loaded into common)
         newf = 0
	 call ocstart(ir1,im1)
	 if(ispoly) call plstart
      endif
c ind is cell number; x and y are positions within cell, relative to center,
c in E and N, in fractions of a cell dimension
      call celfnd(rlong,rlato,x,y,i,j,ind)
      if(ind.gt.latc*longc.or.ind.le.0) then
         amp = cz
         return
      endif
c if cell is inside a polygon we don't want, or outside one we do, reject
c it without going further
      if(ispoly) then
	 call chkgon(rlong,rlato,iok)
	 if(iok.eq.0) then
            amp = cz
            return
         endif
      endif
      if(fingrd.eq.'C'.or.fingrd.eq.'G') then
c
c  large distance - return the model value (0 if on land for tidal models)
c
         amp = .001*cmplx(float(ir1(ind)),float(im1(ind)))
         if(fingrd.eq.'C') call seadens(rlato,rlong,rho)
         if(fingrd.eq.'G') rho = 1.e3
         amp=rho*amp
         rhol=rho
         return
      endif
      if(fingrd.eq.'F'.or.fingrd.eq.'L') then
c
c  close distance, or land load - use the
c  5-minute map file to decide if there is land or not. If there is
c  not (for F) or is (for L), bilinearly interpolate from the model,
c  filling out empty cells to make this possible (though if all
c  the cells needed are zero, return this).
c
         call lndsea(rlato,rlong,lnd)
         if((lnd.eq.1.and.fingrd.eq.'F').or.
     1      (lnd.eq.2.and.fingrd.eq.'L')) then
c  detailed grid says land (for F) or sea (for L) - return zero
            amp = cz
            return
         endif
c for L option, do not try to interpolate cells, density is 1.e3
         if(fingrd.eq.'L')then
           amp = .001*cmplx(float(ir1(ind)),float(im1(ind)))
           amp = 1.e3*amp
           return
         endif
c  now left with option F, and detailed map says ocean
c  load an interpolating array of four elements with values from adjacent
c  cells, depending on which quadrant of the cell we are in.
c  In the process note how many nonzero values we have; save one
         ix = 1
	 if(x.le.0) ix = -ix
         iy = -longc
	 if(y.le.0) iy = -iy
	 x=abs(x)
	 y=abs(y)
	 nnz=0
	 do i=1,4
	   if(i.eq.1) ibind=ind
	   if(i.eq.2) ibind=ind+ix
	   if(i.eq.3) ibind=ind+ix+iy
	   if(i.eq.4) ibind=ind+iy
           if(ibind.le.0.or.ibind.gt.latc*longc)
     1      bg(i) = cz
           if(ibind.gt.0.and.ibind.le.latc*longc)
     1      bg(i) = 0.001*cmplx(float(ir1(ibind)),float(im1(ibind)))
	   if(bg(i).ne.cz) nnz = nnz+1
	   if(bg(i).ne.cz) bgs = bg(i)
         enddo
c
c  now we have a grid for interpolation: the values are arranged as
c
c           4      3
c           1      2
c
c  where (1) is the value we are closest to.
c
c  we now fill out the zero values (if any) from the nonzero ones
c  (Except that if there are 0 or 4 nonzero values we needn't bother)
         if(nnz.eq.1) then
c  the one saved value is the nonzero one: fill out the whole array with this
            do i=1,4
              bg(i) = bgs
            enddo
	 endif
         if(nnz.eq.2) then
c  two values: either in a row, a column, or diagonal:
            if(bg(1).ne.cz.and.bg(2).ne.cz) then
               bg(4) = bg(1)
               bg(3) = bg(2)
            elseif(bg(2).ne.cz.and.bg(3).ne.cz) then
               bg(4) = bg(3)
               bg(1) = bg(2)
            elseif(bg(3).ne.cz.and.bg(4).ne.cz) then
               bg(1) = bg(4)
               bg(2) = bg(3)
            elseif(bg(4).ne.cz.and.bg(1).ne.cz) then
               bg(3) = bg(4)
               bg(2) = bg(1)
c     diagonals
            elseif(bg(4).ne.cz.and.bg(2).ne.cz) then
               bg(3) = 0.5*(bg(4)+bg(2))
               bg(1) = bg(3)
            elseif(bg(1).ne.cz.and.bg(3).ne.cz) then
               bg(4) = 0.5*(bg(1)+bg(3))
               bg(2) = bg(3)
            endif
	 endif
	 if(nnz.eq.3) then
c one missing value
	    if(bg(1).eq.cz) bg(1) = bg(4)+bg(2)-bg(3)
	    if(bg(2).eq.cz) bg(2) = bg(1)+bg(3)-bg(4)
	    if(bg(3).eq.cz) bg(3) = bg(4)+bg(2)-bg(1)
	    if(bg(4).eq.cz) bg(4) = bg(1)+bg(3)-bg(2)
	 endif
c  bilinear interpolation
	 amp = bg(1)*(1-x)*(1-y) + bg(2)*x*(1-y) +
     1         bg(4)*(1-x)*y     + bg(3)*x*y
c  if all the values are zero, we will return zero - but we pass some
c  information about this through common block coninf; the first time this
c  happens, set clnf and farnf to the current distance from the station,
c  and afterwards only farnf (since we always increase the distance).
	 if(nnz.eq.0) then
	    if(numnf.eq.0) then
	       clnf = dist
	       farnf = dist
	    else
	       farnf = dist
	    endif
	    numnf=numnf+1
	 endif
      endif
c
c  scale the amplitude by the sea-water density
c
      rho=0.0
      call seadens(rlato,rlong,rho)
      amp=rho*amp
      rhol=rho
      return
      end
