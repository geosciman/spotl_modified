c
c  Original date of RCS archived version is  Sep  7  1988
c
c
c     $Id: shells.f,v 1.1 2011/09/23 22:39:02 agnew Exp agnew $
c
c     $Log: shells.f,v $
c     Revision 1.1  2011/09/23 22:39:02  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine shells(x,k,n)
c  sorts an array x, of length n, sorting upward, and returns
c  an array k which may be used to key another array to the
c  sorted pattern (i.e., if we had an array f to which x 
c  corresponded before sorting, then after calling shells,
c  f(k(1)) will be the element of f corresponding to the
c  smallest x, f(k(2)) the next smallest, and so on.
c   revised 29-dec-82 so k is sorted in turn, the values of
c  k that point to identical values of x being put in increasing
c  order
c
c    calls no other routines
c
      dimension x(n),k(n)
      igap = n
      do 1 i = 1,n
 1    k(i) = i
 5    if(igap.le.1) go to 25
      igap = igap/2
      imax = n - igap
 10   iex = 0
      do 20 i = 1,imax
      ipl = i + igap
      if(x(i).le.x(ipl)) go to 20
      sv = x(i)
      ik = k(i)
      x(i) = x(ipl)
      k(i) = k(ipl)
      x(ipl) = sv
      k(ipl) = ik
      iex = iex + 1
 20   continue
      if(iex.gt.0) go to 10
      go to 5
c
c  now sort k's (for identical values of x, if any)
c
 25   j = 1
 30   if(j.ge.n) return
      if(x(j).eq.x(j+1)) go to 33
      j = j + 1
      go to 30
c  have at least two x's with the same value - see how long this is true
 33   l = j
 35   if(x(l).ne.x(l+1)) go to 38
      l = l + 1
      if(l.lt.n) go to 35
c  j and l are the indices within which x(i) does not change - sort k
 38   igap = l - j + 1
 40   if(igap.le.1) j = l + 1
      if(igap.le.1) go to 30
      igap = igap/2
      imax = l-j+1 - igap
 45   iex = 0
      do 50 i=1,imax
      ipl = i + igap + j - 1
      if(k(i+j-1).le.k(ipl)) go to 50
      ik = k(i+j-1)
      k(i+j-1) = k(ipl)
      k(ipl) = ik
      iex = iex + 1
 50   continue
      if(iex.gt.0) go to 45
      go to 40
      end

