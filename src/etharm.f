c
c  Original date of RCS archived version is  May 22  1998
c
c
c     $Id: etharm.f,v 1.1 2011/09/23 22:39:01 agnew Exp agnew $
c
c     $Log: etharm.f,v $
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine etharm(cta,n,m,pot,grav,disp,tlt,str)
c  For a Cartwright-Tayler constituent of amp cta, order n, species m,
c  returns the potential height, gravity, tilt (ns and ew) and strain
c  (ns, ew and shear) for a point whose position is specified in
c  common block obs (previously filled using subroutine sph).
c  Except for the potential, the amplitudes are given with sign adjusted
c  so the local potential would be positive.
c  Note that tlt(2) and str(3) are in quadrature to
c  the other two. tlt(1) is north; tlt(2) is east; str(1) is north;
c  str(2) is east; str(3) is se (theta-lambda). disp(1) is vertical,
c  disp(2) is north, disp(3) is east
c
c      calls ifact,p (both in this source file)
c
      dimension h(2),sk(2),sl(2),tlt(2),str(3),disp(3)
      common/obs/ct,st,clg,slg,dv,r,g
      data pi/3.14159265/,gm/9.7982776/,ae/6378140./
c  Love numbers (should be kept the same as in ertid for consistency).
      data h/.6114,.2891/,sk/.304,.09421/,sl/.0832,.0145/
c  ynm normalizing factor
      j = n-1
      rn = ((2*n+1)*ifact(n-m))/(4*pi*ifact(n+m))
      rn = sqrt(rn)
      if(mod(m,2).ne.0) rn = -rn
      potr =  cta*(gm/g)*(r**n)*rn*p(n,m,ct)
c potential height (relative to crust) in mm
      pot =  1000*potr*(1+sk(j)-h(j))
      grav = cta*(gm/ae)*n*(r**(n-1))*rn*p(n,m,ct)
c gravity in microgal
      grav = 1.e8*grav*(1.+ (2./n)*h(j) - ((n+1.)/n)*sk(j))
      sfac = (gm/(g*ae*(st**2)))*(r**(n-1))*rn
c strain in nanostrain
      str(1) = 1.e9*cta*sfac*(
     1  (h(j)*(st**2) + sl(j)*((n*ct)**2-n))*p(n,m,ct)
     2            -2.*sl(j)*(n-1)*(n+m)*ct*p(n-1,m,ct)
     3               + sl(j)*(n+m)*(n+m-1)*p(n-2,m,ct))
      str(2) = 1.e9*cta*sfac*(
     1  (h(j)*(st**2) + sl(j)*(n*(ct**2)-m**2))*p(n,m,ct)
     2                        -sl(j)*(n+m)*ct*p(n-1,m,ct))
      str(3) = 1.e9*cta*sfac*m*sl(j)*(
     1  (n-1)*ct*p(n,m,ct) - (n+m)*p(n-1,m,ct))
c tilt in nanorad
      tlt(1) = -1.e9*(1.+sk(j)-h(j))*cta*st*sfac*(n*ct*p(n,m,ct)
     1                           - (n+m)*p(n-1,m,ct))
      tlt(2) =  1.e9*(1.+sk(j)-h(j))*cta*st*sfac*m*p(n,m,ct)
      a = r*ae
c displacement in mm
      disp(1) = 1.e3*h(j)*potr
      disp(2) = (a*sl(j)/(1.+sk(j)-h(j)))*(tlt(1)/1.e6)
      disp(3) = (a*sl(j)/(1.+sk(j)-h(j)))*(tlt(2)/1.e6)
      return
      end
      function p(n,m,c)
c  computes some low-order legendre polynomials
      s = sqrt(1.-c*c)
      if(n.lt.m) p = 0
      if(n.eq.0.and.m.eq.0) p = 1.
      if(n.eq.1.and.m.eq.0) p = c
      if(n.eq.1.and.m.eq.1) p = s
      if(n.eq.2.and.m.eq.0) p = .5*(3*c*c - 1.)
      if(n.eq.2.and.m.eq.1) p = 3.*s*c
      if(n.eq.2.and.m.eq.2) p = 3.*s*s
      if(n.eq.3.and.m.eq.0) p =  .5*c*(5.*c*c - 3.)
      if(n.eq.3.and.m.eq.1) p = 1.5*s*(5.*c*c - 1.)
      if(n.eq.3.and.m.eq.2) p = 15.*s*s*c
      if(n.eq.3.and.m.eq.3) p = 15.*(s**3)
      return
      end
      function ifact(n)
      ifact = 1
      if(n.le.1) return
      do 1 i=2,n
 1    ifact = i*ifact
      return
      end
