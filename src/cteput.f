c
c  Original date of RCS archived version is  Mar 13  1995
c
c
c     $Id: cteput.f,v 1.1 2011/09/23 22:39:01 agnew Exp agnew $
c
c     $Log: cteput.f,v $
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine cteput(icc,amp,n,m)
c  Returns the amplitude, degree, and order of a tidal
c  constituent whose Doodson number (in C-T form) is input in array icc.
c  The constituents available are stored in the arrays idd (doodson
c  number) and tamp (Cartwright-Edden amplitude); these are the largest
c  constituents of the expansion (24 semidiurnl, 30 diurnal, 9 long-period,
c  4 terdirunal).
c
      dimension icc(6)
c  arrays containing information about all stored constituents
      dimension idd(6,67),tamp(67)
      data nt/67/
      data tamp/
     1  .63192, .29400, .12099, .07966, .02383,-.02358, .02298, .01932,
     2 -.01786, .01720, .01601, .00467,-.00466,-.00451, .00447, .00447,
     3  .00259,-.00246,-.00218, .00197, .00195, .00192,-.00190, .00180,
     4  .36878,-.26221,-.12203,-.05020, .05001,-.04945, .02062, .02062,
     5  .01129,-.00954,-.00947,-.00802, .00741,-.00730, .00723,-.00714,
     6 -.00664, .00525, .00414, .00409, .00395, .00394, .00343, .00342,
     7  .00293, .00289, .00216, .00194,-.00194,-.00180,-.06663,-.03518,
     8 -.02762,-.01276,-.00673,-.00583,-.00529,-.00288,-.00258, .00765,
     9  .00210, .00100,-.00043/
      data idd/
     1 2, 0, 0, 0, 0, 0,  2, 2,-2, 0, 0, 0,  2,-1, 0, 1, 0, 0,  2, 2, 0,
     2 0, 0, 0,  2, 2, 0, 0, 1, 0,  2, 0, 0, 0,-1, 0,  2,-1, 2,-1, 0, 0,
     3 2,-2, 2, 0, 0, 0,  2, 1, 0,-1, 0, 0,  2, 2,-3, 0, 0, 1,  2,-2, 0,
     4 2, 0, 0,  2,-3, 2, 1, 0, 0,  2, 1,-2, 1, 0, 0,  2,-1, 0, 1,-1, 0,
     5 2, 3, 0,-1, 0, 0,  2, 1, 0, 1, 0, 0,  2, 2, 0, 0, 2, 0,  2, 2,-1,
     6 0, 0,-1,  2, 0,-1, 0, 0, 1,  2, 1, 0, 1, 1, 0,  2, 3, 0,-1, 1, 0,
     7 2, 0, 1, 0, 0,-1,  2, 0,-2, 2, 0, 0,  2,-3, 0, 3, 0, 0,  1, 1, 0,
     8 0, 0, 0,  1,-1, 0, 0, 0, 0,  1, 1,-2, 0, 0, 0,  1,-2, 0, 1, 0, 0,
     9 1, 1, 0, 0, 1, 0,  1,-1, 0, 0,-1, 0,  1, 2, 0,-1, 0, 0,  1, 0, 0,
     1 1, 0, 0,  1, 3, 0, 0, 0, 0,  1,-2, 2,-1, 0, 0,  1,-2, 0, 1,-1, 0,
     2 1,-3, 2, 0, 0, 0,  1, 0, 0,-1, 0, 0,  1, 1, 0, 0,-1, 0,  1, 3, 0,
     3 0, 1, 0,  1, 1,-3, 0, 0, 1,  1,-3, 0, 2, 0, 0,  1, 1, 2, 0, 0, 0,
     4 1, 0, 0, 1, 1, 0,  1, 2, 0,-1, 1, 0,  1, 2,-2, 1, 0, 0,  1, 0, 2,
     5-1, 0, 0,  1,-1, 2, 0, 0, 0,  1, 3,-2, 0, 0, 0,  1, 1, 1, 0, 0,-1,
     6 1, 1,-1, 0, 0, 1,  1, 4, 0,-1, 0, 0,  1, 0,-2, 1, 0, 0,  1,-4, 2,
     7 1, 0, 0,  1,-2, 2,-1,-1, 0,  0, 2, 0, 0, 0, 0,  0, 1, 0,-1, 0, 0,
     8 0, 2, 0, 0, 1, 0,  0, 3, 0,-1, 0, 0,  0, 1,-2, 1, 0, 0,  0, 2,-2,
     9 0, 0, 0,  0, 3, 0,-1, 1, 0,  0, 2, 0,-2, 0, 0,  0, 2, 0, 0, 2, 0,
     1 3,0,0,0,0,0,  3,-1,0,1,0,0,  3,2,0,0,0,0,  3,0,0,0,-1,0/
c  see if Doodson numbers match
      do 15 kk=1,nt
      ii = 0
      do 14 i=1,6
 14   ii = ii + iabs(idd(i,kk)-icc(i))
      if(ii.eq.0) go to 16
 15   continue
         write(6,100) icc
 100     format('Doodson number ',6i2,' not available in list')
         stop
c  have a match
 16   amp = tamp(kk)
c  all order 0 1 2 tides in this list are degree 2
      if(icc(1).le.2) n = 2
      if(icc(1).eq.3) n = 3
      m = icc(1)
      return
      end
