c
c  Original date of RCS archived version is  Jul  3  2005
c
c
c     $Id: ertid.f,v 1.1 2011/09/23 22:39:01 agnew Exp agnew $
c
c     $Log: ertid.f,v $
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
c     program erthtd
c
c     BOMM subroutine written by Jon Berger  November 1969
c     astronomy revised by Judah Levine (after J C Harrison) Sept 1973
c     updated by Karen Young  March 1977, for PDP 11
c     emended by Duncan Agnew June 1978 to speed execution
c     solar astronomy redone and rigid-earth potential added by Duncan
c     Agnew June 1979
c     tide generating part put in separate subroutine and computation
c     of Munk-Cartwright coefficients added by Duncan Agnew Feb 1981,
c     July 1982
c
c      This version rewritten for export, using F77 calls for I/O,
c     by Duncan Agnew Mar 1987
c
c     tides are calculated from harmonics 2 through 3 for
c     the lunar terms, 2 for the solar terms.
c     love numbers h(n), k(n), and l(n) are set by data statements.
c
c     gravity tide is in microgals, plus for up acceleration
c     rigid-earth potential is in meters (v/g), plus for up
c     tilt tide is in nanoradians (nr)......  1 nr = 2.063e-4 arcsec
c     strain tide is in parts in e9, plus for extension
c
c     azimuths for strain and tilt are reckoned
c     positive clockwise from north (east from north)
c
      character*1 ioptn
      character*40 filnam
      double precision tt,ts,te,d,zhr,yhr,t
      real iyr,jyr,iday,jday,k,l,lamda
      dimension serout(60,7),ntw(3),tilt(2),strain(3),a(14),iget(7),
     1 ispc(4),filnam(7)
      common/bpos/dsz,dcs,dsl,dcl,ssz,scz,ssl,scl,dpar,sdist
      common/love/h(3),k(3),l(3)
      common/azimut/azt(3,2),azs(3,3)
      common/tdiff/etmut
      common/sunny/moon
      data h/.6114,.2891,.175/,k/.304,.09421,.043/,
     1 l/.0832,.0145,.0103/
      data irl/0/,iflag/0/,ntotl/0/,iget/7*0/,ispc/4*0/,ntw/3*0/,
     1 moon/0/
c
c   data statements for input and output unit numbers (on terminal I/O)
c    and blocklength (in reals) of the data arrays
c
      data inun/5/,ioun/6/,nptpb/60/
      write(ioun,701)
701   format(' Start time- - -  (year,day,hour)   ',t45,$)
      read(inun,*) iyr,iday,zhr
      iyr=iyr-1900
      write(ioun,702)
702   format(' End time- - -  (year,day,hour)  ',t45,$)
      read(inun,*) jyr,jday,yhr
      jyr=jyr-1900
      write(ioun,703)
703   format(' Time interval (hours)-      ',$)
      read(inun,*) d
c   find times in hours from 0 hr, 1 jan 1900
      ts=zhr+24.d0*(iday-1)+8760.d0*iyr+24.d0*int((iyr-1)/4)
      te=yhr+24.d0*(jday-1)+8760.d0*jyr+24.d0*int((jyr-1)/4)
 5    write(ioun,7035)
 7035 format(' Type t for theoretical tides, m for Munk-'
     1 ,'Cartwright coefficients ',$)
      read(inun,*) ioptn
      if(ioptn.ne.'m'.and.ioptn.ne.'t') go to 5
      if(ioptn.eq.'t') go to 69
 6    write(ioun,7036)
 7036 format(' Type 2 for order 2, 3 for 3, 4 for both ',$)
      read(inun,*) iord
      if(iord.lt.2.or.iord.gt.4) go to 6 
      nspc = min0(iord+1,4)
      write(ioun,7037)
 7037 format(' Type 1 to get long-period tides ',$)
      read(inun,*) ispc(1)
      write(ioun,7038)
 7038 format(' Type 1 to get diurnal tides ',$)
      read(inun,*) ispc(2)
      write(ioun,7039)
 7039 format(' Type 1 to get semidiurnal tides ',$)
      read(inun,*) ispc(3)
      if(nspc.eq.4) write(ioun,7040)
 7040 format(' Type 1 to get terdiurnal tides ',$)
      if(nspc.eq.4) read(inun,*) ispc(4)
      do 20 i=1,nspc
      if(ispc(i).le.0) go to 20
      if((iord.eq.2.or.iord.eq.4).and.i.le.3) iget(i)=1
      if(iord.ge.3) iget(i+3) = 1
 20   continue
      ntotl = 0
      do 21 i=1,7
 21   if(iget(i).eq.1) ntotl=ntotl+1
      go to 741
 69   write(ioun,704)
704   format(' North latitude-      ',$)
      read(inun,*) theta
      write(ioun,705)
705   format(' East longitude-      ',$)
      read(inun,*) lamda
      write(ioun,706)
 706  format(' Number of gravity (0 to 1), tilt (0-2), and'
     $,' strain (0-3) tides wanted')
      do 710 i=1,3
      read(inun,*) ntw(i)
  710 ntotl = ntotl + iabs(ntw(i))
      if(ntw(2).eq.0) go to 732
      write(ioun,707)
 707  format(' Tiltmeter azimuth(s) - positive clockwise from'
     $ ,' north ',$)
      do 731 i = 1,ntw(2)
      read(inun,*) azt(1,i)
      azt(2,i) = cos(azt(1,i)/57.2958)
 731  azt(3,i) = sin(azt(1,i)/57.2958)
 732  if(ntw(3).eq.0) go to 741
      write(ioun,708)
 708  format(' Strainmeter azimuth(s) - positive clockwise from'
     $ ,' north ',$)
      do 733 i = 1,ntw(3)
      read(inun,*) azs(1,i)
      azs(2,i) = cos(azs(1,i)/57.2958)
 733  azs(3,i) = sin(azs(1,i)/57.2958)
 741  write(ioun,709) ntotl
 709  format(' Enter output filenames (total of ',i1,')')
c
c  open direct-access files, nptpb real*4's per record
c
      do 7099 i=1,ntotl
      read(inun,7098) filnam(i)
 7098 format(a40)
      iunit = i
c skip units 5 and 6
      if(i.ge.5) iunit = i+2
      open(unit=iunit,file=filnam(i))
CC     $ form='unformatted',recl=4*nptpb)
 7099 continue
      write(ioun,903) ts,te,d
 903  format(' Time runs from',f12.4,' to',f12.4,
     $' with interval',f12.4)
c
c     done asking questions - begin execution
c
      i = 1
      tt = ts
      call sph(theta,lamda,0.)
      etmut = 41.184 + iyr - 70
 10   if(tt.gt.te) go to 749
      t = (tt+12.d0 + (etmut/3600.d0))/876600.d0
c   t is ephemeris time in julian centuries from 12 hr 0 jan 1900
      call ephem(t)
      if(ioptn.eq.'m') go to 80
c
c     calculate normailzed total tides
c
      call elastd(grav,tilt,strain,ntw)
      if(ntw(1).eq.0) go to 71
c   gravity
      serout(i,1) = 1.e8*grav
      if(ntw(1).eq.-1) serout(i,1) = grav
 71   if(ntw(2).eq.0) go to 73
c   tilt
      do 72 kk = 1,ntw(2)
      idd = iabs(ntw(1)) + kk
 72   serout(i,idd) = 1.e9*tilt(kk)
 73   if(ntw(3).eq.0) go to 3
c   strain
      do 74 kk = 1,ntw(3)
      idd = iabs(ntw(1)) + ntw(2) + kk
 74   serout(i,idd) = 1.e9*strain(kk)
      go to 3
 80   call mcco(a,a(8))
      idd = 0
      do 81 kk =1,7
      if(iget(kk).eq.0) go to 81
      idd = idd + 1
      serout(i,idd) = a(kk)
      serout(i+1,idd) = a(7+kk)
 81   continue
c  double increment for mc coefficients
      i = i + 1
 3    i = i + 1
      tt = tt + d
      if(i.eq.nptpb+1) go to 750
      go to 10
 749  iflag = 1
      if(i.eq.1) go to 770
c
c      write data out
c
750   do 751 kk = 1,ntotl
      npts = nptpb
      if(iflag.eq.1) npts = i - 1
      iunit = kk
      if(kk.ge.5) iunit = kk+2
CCC      write(iunit,rec=irl) (serout(jj,kk),jj=1,nptpb)
      write(iunit,7501) (serout(jj,kk),jj=1,npts)
 7501 format(f13.5)
 751  continue
      i = 1
      irl = irl + npts
      if(iflag.eq.0) go to 10
770   write(ioun,771) irl
771   format(//' each output series is ',i8,' terms long'/)
      stop
      end
      subroutine sph(grlat,elong,ht)  
c   for a point at geographical north latitude grlat, east longitude elong      
c   (in degrees), and height ht (in meters), finds geocentric position   
c   and local g using formulae for a spheroid
c  
      common/obs/cthet,sthet,clong,slong,dvert,radn,g      
      data gn,ae,f,rm,dr/9.798277692,6378140.,.00335281,.00344978,
     $ .01745329252/    
      clong = cos(elong*dr)    
      slong = sin(elong*dr)    
c   latitude difference 
      dvert = f*(1.+.5*f)*sin(2.*grlat*dr) - .5*f*f*sin(4.*grlat*dr)     
      gcclat = (3.1415926535898/2.) - (grlat*dr - dvert)   
      cthet = cos(gcclat)      
      sthet = sin(gcclat)      
c   geocentric radius   
      radn = 1. - f*(cthet**2)*(1. + 1.5*f*(sthet**2))     
c   formulae for g are from jeffreys, 4.022 and 4.023      
      g = gn*(1. + f - 1.5*rm + f*(f-(27./14.)*rm) + (2.5*rm - f -f*(f-  
     $ (39./14.)*rm))*(cthet**2) - (f/2.)*(7.*f-15.*rm)*((cthet*sthet)   
     $**2))      
c   free air correction 
      g = g - g*(2.*ht*(1.+f+rm - 2.*f*(cthet**2))/ae)     
      return     
      end 
      subroutine ephem(t)
      double precision t
c   t is ephemeris time in julian centuries from 12 hr 0 jan 1900
c  (in ordinary civil reckoning this is greenwich noon on 31 december
c  1899). if the difference between ephemeris and unversal time is
c  not put in (see common tdiff below), t should be in universal time
c  which is (nearly) the time ordinarily kept by clocks.
c   computes positions of the sun and moon at time t, returning results
c  in common block bpos. the solar ephemeris uses the mean sun.
c    Derived from J. Levine's revision (after J. C. Harrison)
c  of an earthtide program by J. Berger and W. E. farrell, with small
c  alterations by D. C. Agnew, partly after M. Wimbush. present
c  subroutine version by d. c. agnew.
      double precision ts,hr,pi20,t2,ll,xxx,hs,ps,es,psig,chmp,shmp,
     $ ls,sz,cz,psis,hm,pm,nm,bl,bls,bf,bd,tlatm,tlongm,plx,sinmla,
     $ cosmla,sinmln,cosmln,at1,at2,ram,w,sinw,cosw,sls
c  common block bpos contains, in order:
c    sine and cosine of colatitude of sublunar point
c    sine and cosine of east longitude of sublunar point
c    sine and cosine of colatitude of subsolar point
c    sine and cosine of east longitude of subsolar point
c    the lunar sine parallax in arc seconds
c    the solar distance in astronomical units
      common/bpos/dsz,dcz,dsl,dcl,ssz,scz,ssl,scl,dpar,sdist
c
c  common block containing the difference ephemeris minus
c  universal time, in seconds. if this is not known it should
c  be set to zero, and the argument to the program should be
c  universal rather than ephemeris time.
c
      common/tdiff/etmutc
c  common block containing the instruction on which ephemeris to compute
c  moon =   0  - both sun and moon
c           1  - moon only
c           2  - sun only
      common/sunny/moon
      fcos(xxx) = cos(sngl(xxx))
      fsin(xxx) = sin(sngl(xxx))
      data pi20/62.8318530717958d0/
c   compute universal time in hours
      ts = 876600.d0*t -12.d0 - (etmutc/3600.d0)
      hr = dmod(ts,24.d0)
c   compute obliquity of the ecliptic
      w = .409319747d0 - .0002271107d0*t
      cosw = dcos(w)
      sinw = dsin(w)
      t2 = t*t
      if(moon.eq.1) go to 7
c
c     compute solar constants for given t
c
      hs = 4.881627982482d0 + 628.3319508731d0*t + 0.523598775578d-5*t2
      hs = dmod(dmod(hs,pi20)+pi20,pi20)
      ps = 4.908229466993d0+ 0.03000526416690d0*t + 0.790246300201d-5*t2
      es = 0.01675104d0 - 0.00004180d0*t - 0.000000126d0*t2
      psig = 0.2617993877971d0*(hr-12.) + hs
      chmp = dcos(hs-ps)
      shmp = dsin(hs-ps)
      ls = hs + shmp*es*(2.+2.5*es*chmp)
      sls = dsin(ls)
      cz = sinw*sls
      sz = dsqrt(1.-cz**2)
      psis = datan2(cosw*sls,dcos(ls))
      rbarr = 1. + es*(chmp + es*(chmp-shmp)*(chmp+shmp))
      ll = psis - psig
      scz = cz
      ssz = sz
      ssl = sin(ll)
      scl = cos(ll)
      sdist = 1./rbarr
c
c     compute lunar constants for given t
c
 7    if(moon.eq.2) return
      hm=4.7199666d0+8399.7091449d0*t-.0000198d0*t2
      pm=5.83515154d0+71.01804120839d0*t-.180205d-3*t2
      nm=4.523601515d0-33.75714624d0*t+.3626406335d-4*t2
c.......  bl bls bf bd are the fundamental arguments of browns theory
      bl=hm-pm
      bls=hs-ps
      bf=hm-nm
      bd=hm-hs
c.........lunar lat long and parallax from brown.  latter two from
c......improved lunar ephemeris, latitude from ras paper of 1908...
      tlongm=hm+.10976*fsin(bl)-.02224*fsin(bl-2.*bd)+0.01149*fsin(2.*bd
     1)+0.00373*fsin(2.*bl)-.00324*fsin(bls)-.00200*fsin(2.*bf)-0.00103*
     2fsin(2.*bl-2.*bd)-.00100*fsin(bl+bls-2.*bd)+0.00093*fsin(bl+2.*bd)
     3-.00080*fsin(bls-2.*bd)+.00072*fsin(bl-bls)-.00061*fsin(bd)-.00053
     4*fsin(bl+bls)
      tlatm=.08950*fsin(bf)+.00490*fsin(bl+bf)-.00485*fsin(bf-bl)-.00303
     1 *fsin(bf-2.*bd)+.00097*fsin(2.*bd+bf-bl)-.00081*fsin(bl+bf-2.*bd)
     2+.00057*fsin(bf+2.*bd)
      plx=(3422.45+186.54*fcos(bl)+34.31*fcos(bl-2.*bd)+28.23*fcos(2.*bd
     1)+10.17*fcos(2.*bl)+3.09*fcos(bl+2.*bd)+1.92*fcos(bls-2.*bd)+1.44*
     2fcos(bl+bls-2.*bd)+1.15*fcos(bl-bls)-0.98*fcos(bd)-0.95*fcos(bl+
     3bls)-0.71*fcos(bl-2.*bf)+0.62*fcos(3.*bl)+0.60*fcos(bl-4.*bd))
      sinmla=dsin(tlatm)
      cosmla=dcos(tlatm)
      sinmln=dsin(tlongm)
      cosmln=dcos(tlongm)
c...convert from celestial lat and long according to explan suppl of
c......na and le page 26
      cz=cosmla*sinmln*sinw+sinmla*cosw
      sz=dsqrt(1.-cz**2)
      at1=cosmla*sinmln*cosw-sinmla*sinw
      at2=cosmla*cosmln
      ram=datan2(at1,at2)
      ll=ram-psig
      dcz = cz
      dsz = sz
      dsl = sin(ll)
      dcl = cos(ll)
      dpar = plx
      return
      end
      subroutine mcco(a,b)     
c   computes a and b tidal potential coefficients as in Munk and Cartwright     
c    1966, in meters for a sphere with the earths equatorial radius.     
c   coefficients are returned in a and b, in order of increasing l (l = 2 or 3) 
c    and for each l in order of increasing m.  b(1) and b(4) are zero.   
c  
      dimension a(7),b(7),rk(4),c(7),cl(4),sl(4),pa(2),rnorm(7)   
c   common with solar and lunar positions    
      common/bpos/coor(8),par,sdist   
c   tidal coefficients are from 1976 iau constants  
      data rk/.3583728,.0059463,.16457875,.00000702/,cl(1)/1./,sl(1)/0./ 
      data ifl/0/,fpi/12.566370614/,c(1)/1./ 
c   on first call compute normalizing coefficiants  
      if(ifl.ne.0) go to 3     
      ifl = 1    
      do 1 i = 2,7      
 1    c(i) = (i-1.)*c(i-1)     
      do 2 l = 2,3      
      mu = l + 1 
      do 2 m = 1,mu     
      i = 3*(l/3) + m   
      t = 1. + (m+1)/3  
 2    rnorm(i) =((-1)**(m+1))*t*sqrt((fpi*c(l-m+2))/((2.*l+1)*c(l+m)))   
c   find normalized parallax   
 3    pa(1) = par/3422.448     
      pa(2) = 1./sdist  
      do 4 i = 1,7      
      a(i) = 0.  
 4    b(i) = 0.  
      do 6 n = 1,2      
c   n = 1 for lunar part, 2 for solar 
      do 5 i = 1,3      
      cl(i+1) = coor(4*n)*cl(i) - coor(4*n-1)*sl(i) 
 5    sl(i+1) = coor(4*n)*sl(i) + coor(4*n-1)*cl(i) 
      cz = coor(4*n-2)  
      sz = coor(4*n-3)  
      c(1) = (3.*cz*cz - 1.)/2.
      c(2) = 3.*cz*sz   
      c(3) = 3.*sz*sz   
      c(4) = (5./3.)*cz*(c(1) - .4)   
      c(5) = 1.5*sz*(5.*cz*cz - 1.)   
      c(6) = 5.*cz*c(3) 
      c(7) = 5.*sz*c(3) 
      do 6 l = 2,3      
      mu = l + 1 
      do 6 m = 1,mu     
      k = 3*(l/3) + m   
      c(k) = rnorm(k)*rk(2*n+l-3)*(pa(n)**(l+1))*c(k)      
      a(k) = a(k) + c(k)*cl(m) 
 6    b(k) = b(k) + c(k)*sl(m) 
      return     
      end 
      subroutine elastd(grav,tilt,strain,ntw)
      save gdc,tnsdc,etdc,eldc,potdc,re,del,dimin
      real k,l
      dimension tilt(2),strain(3),ntw(3),amrat(2),rbor(2),pa(2),p(3),
     1 pp(3),ppp(3),e(3),rkr(3),tltcor(2),del(3),dimin(3)
c  computes the earth tides on an elastic earth, given the solar and lunar
c  positions and the place of observation.  degree 2 is used for the solar
c  tides, 2 through 4 for the lunar tides. the results are returned in
c  grav, tilt, and strain. ntw(1) gives the number of gravity tides
c  ntw(2) the number of tilt tides, ntw(3) the number of strain tides.
c  as the dimensioning shows, at most one gravity tide, two tilts, and
c  three strains may be computed.  if ntw(1) = -1, the program will
c  put the equilibrium potential height for a rigid spherical earth in
c  grav.  the units are m/s**2, radians, extension, and meters.
c  the sign convention is positive for upward potential
c  height, a decrease in g, a tilt in the azimuth given, and
c  extensional strain.
c  that part of the tidal signal which comes from the permanent
c  deformation is subtracted out, using the coefficient of .31455 m
c  found by Cartwright and Edden for the dc tide.
c
c   based (very closely) on J. Berger earth tide program.
c   converted to subroutine by D. Agnew Nov 1979.
c   computation of the potential and subtraction of dc tides added
c   by D. Agnew Jan 1980.
c
c  common block with observatory information
      common/obs/cth,sth,clg,slg,dif,radn,gl
c  common block with love numbers
      common/love/h(3),k(3),l(3)
c  common block with strainmeter and tiltmeter azimuths
      common/azimut/azt(3,2),azs(3,3)
c  common block with lunar and solar colat and long, lunar sine parallax,
c  and solar distance
      common/bpos/coor(8),par,sdist
c  data for mean parallaxes, a times m over m(earth), equatorial radius.
      data rbor/1.6592496e-2,4.2635233e-5/,amrat/78451.25,
     1  2.1235762e12/,a/6.37814e6/
      data g1/9.79828/,g2/9.82022/,ppp(1)/3./,iflag/0/
c  on first call compute factors for gravity and tilt, and dc tides
c  at the given latitude.
      if(iflag.eq.1) go to 3
      iflag=1
      do 1 i = 1,3
      del(i) = 1. + (2./(i+1.))*h(i) - ((i+2.)/(i+1.))*k(i)
 1    dimin(i) = 1. + k(i) - h(i)
c  dc gravity tide is also known as the Honkasalo correction
c  **note that the love numbers for an elastic earth are used
c    in computing the dc tide as well
      gdc = -3.0481e-7*(3*cth**2 - 1.)*del(1)*radn
      tnsdc = -9.1445e-7*cth*sth*dimin(1)*radn/gl
      etdc = -1.555e-8*(h(1)*(3.*cth**2-1.) - 6.*l(1)*(2.*cth**2 -1.))
      eldc = -1.555e-8*(h(1)*(3.*cth**2-1.) - 6.*l(1)*cth**2)
      potdc = .0992064*(1.-3*cth**2)
      re = 1./(radn*a)
c  zero out arrays
 3    do 4 i = 1,2
      tilt(i) = 0.
      e(i) = 0.
 4    tltcor(i) = 0.
      e(3) = 0.
      grav = 0.
      gnth = 0.
c  compute normalized parallax
      pa(1) = par/3422.45
      pa(2) = 1./sdist
c  in outer loop, ii = 1 for moon, 2 for sun
      do 40 ii = 1,2
      id = 3
      if(ii.eq.2) id = 1
      ir = 4*(ii-1)
c  find cosine of zenith angle, potential constants, legendre polynomials
c  and their derivatives, and derivatives of the cosine of the zenith
c  angle.
      cll = clg*coor(ir+4) + slg*coor(ir+3)
      sll = slg*coor(ir+4) - clg*coor(ir+3)
      cz = coor(ir+2)
      sz = coor(ir+1)
      cu = cth*cz + sth*sz*cll
      xi = rbor(ii)*pa(ii)*radn
      cc = amrat(ii)*rbor(ii)*pa(ii)
      rkr(1) = cc*xi*xi
      rkr(2) = rkr(1)*xi
      rkr(3) = rkr(2)*xi
      p(1) = 0.5*(3*cu*cu - 1.)
      pp(1)=3*cu
      if(ii.eq.2) go to 5
      p(2) = .5*cu*(5.*cu*cu - 3.)
      p(3) = .25*(7.*cu*p(2) - 3.*p(1))
      pp(2) = 1.5*(5.*cu*cu - 1.)
      pp(3) = .25*(7.*p(2) + cu*pp(2)) - 3.*pp(1)
      ppp(2) = 15.*cu
      ppp(3) = 7.5*(7.*cu*cu - 1.)
 5    cut = -sth*cz + cth*sz*cll
      cutt = -cu
      cul = -sth*sz*sll
      cull = -sth*sz*cll
      cutl = -cth*sz*sll
      if(ntw(1).eq.0) go to 15
      do 10 j = 1,id
      if(ntw(1).eq.-1) grav = grav + rkr(j)*p(j)
 10   if(ntw(1).eq.1) grav = grav + del(j)*(j+1)*rkr(j)*p(j)*g1*re
      gnth = gnth - dimin(1)*rkr(1)*pp(1)*g1*cut*re
 15   if(ntw(2).eq.0) go to 30
      do 25 kk = 1,ntw(2)
      do 20 j=1,id
 20   tilt(kk)=tilt(kk)-(dimin(j)*rkr(j)*pp(j)*re)*(cut*azt(2,kk)
     $ - cul*azt(3,kk)/sth)
 25   tltcor(kk) = tltcor(kk) + del(1)*2.*rkr(1)*p(1)*re
 30   if(ntw(3).eq.0) go to 40
      do 35 j = 1,id
      e(1) = e(1) + rkr(j)*(l(j)*(pp(j)*cut*cth/sth + (ppp(j)*cul*cul
     1 + pp(j)*cull)/(sth*sth)) + h(j)*p(j))
      e(2) = e(2) + rkr(j)*(l(j)*(ppp(j)*cut*cut + pp(j)*cutt) +
     1  h(j)*p(j))
      e(3) = e(3) + 2*(rkr(j)*l(j)*(ppp(j)*cul*cut + pp(j)*cutl -
     1 cth*pp(j)*cul/sth))/sth
  35  continue
  40  continue
c  do ellipticity corrections, convert strains to strainmeter
c  azimuths, and remove dc terms.
      if(ntw(1).eq.1) grav = grav + gnth*dif - gdc
      if(ntw(1).eq.-1) grav = grav - potdc
      do 45 i =1,ntw(2)
 45   tilt(i) = tilt(i)*(g1/gl) -tltcor(i)*dif*azt(2,i) -tnsdc*azt(2,i)
      do 50 i = 1,ntw(3)
 50   strain(i) = (e(1)*azs(3,i)**2 - e(3)*azs(3,i)*azs(2,i)
     $+e(2)*azs(2,i)**2)*re*(g1/g2) -etdc*azs(2,i)**2 -eldc*azs(3,i)**2
      return
      end
