c
c  Original date of RCS archived version is  Nov 24  2007
c
c
c     $Id: fndgde.f,v 1.1 2011/09/23 22:39:01 agnew Exp agnew $
c
c     $Log: fndgde.f,v $
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine fndgde(ftlon,rilon,slat,rnlat,il,ih,jl,jh)
c
c  returns the cell indices for the corners of a grid
c   ftlon,rilon are the left and right longitudes wanted,
c   slat and rnlat the south and north limits
c  returns the cell indices that cover these, allowing for
c  the case when the limits given are outside the model; will not
c  give the full limits of indices when the longitudes cross the
c  edge of the model.
c
      common/modlim/tlat,blat,wlong,elong,latc,longc
c  put longitudes into model range
      dlongl = amod(ftlon-wlong,360.)
      if(dlongl.lt.0) dlongl = dlongl + 360
      dlongr = amod(rilon-wlong,360.)
      if(dlongr.lt.0) dlongr = dlongr + 360
c allow for limits completely outside the model
      if(slat.gt.tlat.or.rnlat.lt.blat.or.
     1   dlongl.gt.elong-wlong.or.dlongl.gt.elong-wlong) then
        il=0
	ih=0
	jl=0
	jh=0
      endif
      x = (longc*dlongl)/(elong-wlong)
      y = (latc*(tlat-rnlat))/(tlat-blat)
      il= int(x) + 1
      jl= int(y) + 1 
      if(rnlat.gt.tlat) jl = 1
      x = (longc*dlongr)/(elong-wlong)
      y = (latc*(tlat-slat))/(tlat-blat)
      ih= int(x) + 1
      jh= int(y) + 1
      if(slat.lt.blat) jh = latc
c if NS edge of model is in grid, ih will be le il (though both are > 0)
      return
      end
