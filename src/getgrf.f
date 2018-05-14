c
c  Original date of RCS archived version is  May 21  1998
c
c
c     $Id: getgrf.f,v 1.1 2011/09/23 22:39:01 agnew Exp agnew $
c
c     $Log: getgrf.f,v $
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine getgrf(grfile,num,ntot,ngr,fingrd)
      character*80 grfile
      character*80 grname
      character*1 fingrd,statgr
      save newf,llu
      common/green/beg,end,spc,grfn(200,7)
      common/grinf/nring,rin(10),rout(10),rsz(10),statgr(10),grname
      data newf/0/
      if(newf.eq.0) then
	 newf = 1
         llu =  3
         open(unit=llu,file=grfile,status='old',access='sequential',
     $    form="formatted")
	 read(llu,100) grname
 100     format(a)
      endif
      read(llu,102) ngreen,num,ntot,ngr,beg,end,spc,fingrd
 102  format(i1,i3,2i4,3f10.4,5x,a)
      read(llu,104) ((grfn(ii,j),j=1,7),ii=1,ngr)
 104  format(7e13.6)
      if(fingrd.eq.'f') fingrd = 'F'
      if(fingrd.eq.'c') fingrd = 'C'
      rin(num) = beg
      rout(num) = end
      rsz(num) = spc
      statgr(num) = fingrd
      nring=ntot
      return
      end
