c
c  Original date of RCS archived version is  Feb 23  2006
c
c
c     $Id: ldbxdr.f,v 1.1 2011/09/23 22:39:01 agnew Exp agnew $
c
c     $Log: ldbxdr.f,v $
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine ldbxdr(azimuth,distance,azstp,disstp)
      dimension xx(2,2),yy(2,2)
c  outputs the coordinates (longitude-latitude) of a load cell centered
c  at given azimuth and distance, with dimensions azstp and disstp. The
c  station location is passed to invspt in a common block
      data lu/6/
      do 7 i=1,2
      do 5 j=1,2
      idelsg = 2*i - 3
      ialpsg = 2*j - 3
      cornaz = azimuth+ialpsg*(azstp/2.)
      rad = distance+idelsg*(disstp/2.)
      call invspt(cornaz,rad,b,rlong)
      yy(i,j) = 90-b
      xx(i,j) = rlong
  5   continue
  7   continue
      write(lu,132) xx(1,1),yy(1,1)
      write(lu,132) xx(1,2),yy(1,2)
      write(lu,132) xx(2,2),yy(2,2)
      write(lu,132) xx(2,1),yy(2,1)
      write(lu,132) xx(1,1),yy(1,1)
      write(lu,132) 100.,0.
 132  format(2f13.6)
      return
      end
