c
c  Original date of RCS archived version is  May 13  1989
c
c
c     $Id: inpoly.f,v 1.1 2011/09/23 22:39:01 agnew Exp agnew $
c
c     $Log: inpoly.f,v $
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
      function inpoly(x,y,xv,yv,nv)
c
c   Returns 1 if point at (x,y) is inside polygon whose nv vertices
c   are at xv(1),yv(1);...,xv(nv),yv(nv)
c   Returns 0 if point is outside
c   Returns 2 if point is on edge or vertex
c
c     Rewritten by D. Agnew from the version by Godkin and Pulli,
c     in BSSA, Vol 74, pp 1847-1848 (1984)
c
c    $$$$$calls ncross (in this file)
c
      dimension xv(nv),yv(nv)
      inpoly = 0
      do 5 i=1,nv-1
      isc = ncross(xv(i)-x,yv(i)-y,xv(i+1)-x,yv(i+1)-y)
      if(isc.eq.4) then
c  on edge - know the answer
         inpoly = 2
         return
      endif
      inpoly = inpoly + isc
 5    continue
c  check final segment
      isc = ncross(xv(nv)-x,yv(nv)-y,xv(1)-x,yv(1)-y)
      if(isc.eq.4) then
         inpoly = 2
         return
      endif
      inpoly = inpoly + isc
      inpoly = inpoly/2
c  convert to all positive (a departure from the original)
      inpoly = iabs(inpoly)
      return
      end
      function ncross(x1,y1,x2,y2)
c
c   finds whether the segment from point 1 to point 2 crosses
c   the negative x-axis or goes through the origin (this is
c   the signed crossing number)
c     return value       nature of crossing
c         4               segment goes through the origin
c         2               segment crosses from below
c         1               segment ends on -x axis from below
c                          or starts on it and goes up
c         0               no crossing
c        -1               segment ends on -x axis from above
c                          or starts on it and goes down
c        -2               segment crosses from above
c
      if(y1*y2.gt.0) then
c all above (or below) axis
         ncross = 0
         return
      endif
      c12 = x1*y2
      c21 = x2*y1
      if(c12.eq.c21.and.x1*x2.le.0.) then
c through origin
         ncross = 4
         return
      endif
      if((y1.eq.0.and.x1.gt.0).or.(y2.eq.0.and.x2.gt.0)
     1 .or.((y1.lt.0).and.(c12.gt.c21)).or.((y1.gt.0).and.(c12.lt.c21))
     2 .or.(y1.eq.0.and.y2.eq.0.and.x1.lt.0.and.x2.lt.0)) then
c touches +x axis; crosses +x axis; lies entirely on -x axis
         ncross = 0
         return
      endif
      if(y1.ne.0.and.y2.ne.0) then
c  cross axis
         if(y1.lt.0) ncross = 2
         if(y1.gt.0) ncross = -2
         return
      endif
c one end touches -x axis - goes which way?
      if(y1.eq.0) then
         if(y2.lt.0) ncross = -1
         if(y2.gt.0) ncross = 1
      else
c y2=0 - ends on x-axis
         if(y1.lt.0) ncross = 1
         if(y1.gt.0) ncross = -1
      endif
      return
      end
