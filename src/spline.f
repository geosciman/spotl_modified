c
c  Original date of RCS archived version is  Feb 18  2005
c
c
c     $Id: spline.f,v 1.1 2011/09/23 22:39:02 agnew Exp agnew $
c
c     $Log: spline.f,v $
c     Revision 1.1  2011/09/23 22:39:02  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine spline(nn,x,u,s,a)
c$$$$ calls no other routines
c  finds array  s  for spline interpolator  eval.
c  nn  number of data points supplied (may be negative, see below)
c  x  array containing x-coordinates where function is sampled.  xx(1),xx(2),...
c     must be a strictly increasing sequence.
c  u  array containing sample values that are to be interpolated.
c  s  output array of 2nd derivative at sample points.
c  a  working space array of dimension at least  nn.
c  if the user wishes to force the derivatives at the ends of the series to     
c  assume specified values, he should put du(1)/dx and du(n)/dx in s1,s2
c  and call the routine with nn=-number of terms in series.  normally a parabola
c  is fitted through the 1st and last 3 points to find the slopes.
c  if less than 4 points are supplied, straight lines are fitted.
      dimension x(*),u(*),s(*),a(*)
c
      q(u1,x1,u2,x2)=(u1/x1**2-u2/x2**2)/(1.0/x1-1.0/x2)
c
      n=iabs(nn)
      if (n.le.3) go to 5000
      q1=q(u(2)-u(1),x(2)-x(1),u(3)-u(1),x(3)-x(1))
      qn=q(u(n-1)-u(n),x(n-1)-x(n),u(n-2)-u(n),x(n-2)-x(n))
      if (nn.gt.0) go to 1000
      q1=s(1)
      qn=s(2)
 1000 s(1)=6.0*((u(2)-u(1))/(x(2)-x(1)) - q1)
      n1= n - 1
      do 2000 i=2,n1
      s(i)= (u(i-1)/(x(i)-x(i-1)) - u(i)*(1.0/(x(i)-x(i-1))+
     + 1.0/(x(i+1)-x(i))) + u(i+1)/(x(i+1)-x(i)))*6.0
 2000 continue
      s(n)=6.0*(qn + (u(n1)-u(n))/(x(n)-x(n1)))
      a(1)=2.0*(x(2)-x(1))
      a(2)=1.5*(x(2)-x(1)) + 2.0*(x(3)-x(2))
      s(2)=s(2) - 0.5*s(1)
      do 3000 i=3,n1
      c=(x(i)-x(i-1))/a(i-1)
      a(i)=2.0*(x(i+1)-x(i-1)) - c*(x(i)-x(i-1))
      s(i)=s(i) - c*s(i-1)
 3000 continue
      c=(x(n)-x(n1))/a(n1)
      a(n)=(2.0-c)*(x(n)-x(n1))
      s(n)=s(n) - c*s(n1)
c  back substitiute
      s(n)= s(n)/a(n)
      do 4000 j=1,n1
      i=n-j
      s(i) =(s(i) - (x(i+1)-x(i))*s(i+1))/a(i)
 4000 continue
      return
c  series too short for cubic spline - use straight lines.
 5000 do 5500 i=1,n
 5500 s(i)=0.0
      return
      end
      function eval(y,nn,x,u,s)
c$$$$ calls no other routines
c  performs cubic spline interpolation of a function sampled at unequally
c  spaced intervals.  the routine spline  should be called to set up the array s
c  y  the coordinate at which function value is desired.
c  nn  number of samples of original function.
c  x  array containing sample coordinates. the sequence x(1),x(2).....x(nn)     
c     must be strictly increasing.
c  u  array containing samples of function at the coordinates x(1),x(2)...      
c  s  array containing the 2nd derivatives at the sample points.  found by the  
c     routine  spline, which must be called once before beginning interpolation.
c  if  y  falls outside range(x(1),x(nn))  the value at the nearest endpoint    
c  of the series is used.
      dimension x(1),u(1),s(1)
      common /startx/ istart
      data istart/1/
c
      nn=iabs(nn)
      if (y.lt.x(1))  go to 3000
      if (y.gt.x(nn)) go to 3100
c  locate interval (x(k1),x(k))  containing y
      if (y-x(istart)) 1200,1000,1000
c  scan up the x array
 1000 do 1100 k=istart,nn
      if (x(k).gt.y) go to 1150
 1100 continue
 1150 k1=k-1
      go to 1500
c  scan downwards in x array
 1200 do 1300 k=1,istart
      k1=istart-k
      if (x(k1).le.y) go to 1350
 1300 continue
 1350 k=k1+1
 1500 istart=k1
c  evaluate interpolate
      dy=x(k)-y
      dy1=y-x(k1)
      dk=x(k)-x(k1)
      deli=1./(6.0*dk)
      ff1=s(k1)*dy*dy*dy
      ff2=s(k)*dy1*dy1*dy1
      f1=(ff1+ff2)*deli
      f2=dy1*((u(k)/dk)-(s(k)*dk)/6.)
      f3=dy*((u(k1)/dk)-(s(k1)*dk)/6.0)
      eval=f1+f2+f3
      return
c  out of range.  substitute end values.
 3000 eval=u(1)
      return
 3100 eval=u(nn)
      return
      end


