c
c  Original date of RCS archived version is  Sep  7  1988
c
c
c     $Id: recurs.f,v 1.1 2011/09/23 22:39:02 agnew Exp agnew $
c
c     $Log: recurs.f,v $
c     Revision 1.1  2011/09/23 22:39:02  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine recurs(x,n,hc,nf,om,scr)
      double precision sc,scr
c  does sin and cos recursion to fill in data x, of length n, for
c  nf sines and cosines with frequenciies om (normalized so the
c  nyquist is pi). hc contains alternating cosine and sine coefficients
c  scr is a scratch array of length 3*nf (n.b. - it is double precision)
c
      dimension x(*),hc(*),scr(*),om(*)
c  set up for start of recursion by computing harmonic values
c  at start point and just before it
      do 3 i = 1,nf
      scr(3*i-2) = hc(2*i-1)
      scr(3*i-1) = hc(2*i-1)*cos(om(i)) -hc(2*i)*sin(om(i))
 3    scr(3*i) = 2.*dcos(dble(om(i)))
c  do recursion over data
      do 7 i = 1,n
      x(i) = 0.
c  do recursive computation for each harmonic
      do 5 j  = 1,nf
      x(i) = x(i) + scr(3*j-2)
      sc = scr(3*j-2)
      scr(3*j-2) = scr(3*j)*sc-scr(3*j-1)
 5    scr(3*j-1) = sc
 7    continue
      return
      end

