c
c  Original date of RCS archived version is  Oct 15  2006
c
c
c     $Id: tdfrph.f,v 1.2 2011/11/18 18:38:09 agnew Exp agnew $
c
c     $Log: tdfrph.f,v $
c     Revision 1.2  2011/11/18 18:38:09  agnew
c     added 2008 leap second
c
c     Revision 1.1  2011/09/23 22:39:02  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine tdfrph(idood,freq,phase)
      double precision freq,phase
c
c   Given the Doodson number of a tidal constituent (in idood), returns
c  the frequency and phase.  Phase is returned in degrees and frequency
c  in cycles/day.
c
c   Note that phases must be decreased by 90 degrees if the sum of the order
c  and the species number is odd (as for the 2nd degree diurnals, and 3rd
c  degree low-frequency and semidiurnals).
c   These phases may need further adjustment to allow for the spherical
c  harmonic normalization used; e.g. for that used for the potential by
c  Cartwright and Tayler, 180 degrees must be added for (species,order)
c  = (1,2), (1,3), or (3,3).
c
c  Common block date stores time information, in UT
c
c    calls toymd, leap, juldat, etutc (all for timekeeping)
c
      save itmsave,d,dd
      double precision f1,f2,f3,f4,f5
      double precision fd1,fd2,fd3,fd4,fd5
      double precision d(6),dd(6)
      double precision dayfr,djd,t
      dimension idood(6),itm2(6),itmsave(5)
      common/date/itm(5)
      data itmsave/5*0/
c
c  test to see if time has changed; if so, set the phases and frequencies
c   for each of the Doodson arguments
c
      intial=0
      do i=1,5
        if(itm(i).ne.itmsave(i)) initial=1
      enddo
      if(initial.eq.1) then
        do i=1,5
           itmsave(i) = itm(i)
        enddo
c
c convert times to Julian days (UT) then to Julian centuries from J2000.00
c   (ET)
c
        call toymd(itm,itm2)
        jd = juldat(itm2)
        dayfr=  itm(3)/24.d0 + itm(4)/1440.d0 + itm(5)/84600.d0
        year=itm(1)+(itm(2)+dayfr)/(365+leap(itm(1)))
        call etutc(year,delta)
        djd= jd - 0.5d0 + dayfr
        t = (djd - 2451545.0d0 + delta/86400.d0)/36525.d0
c
c IERS expressions for the Delauney arguments, in degrees
c
        f1 =     134.9634025100d0 +
     1    t*( 477198.8675605000d0 +
     1    t*(      0.0088553333d0 +
     1    t*(      0.0000143431d0 +
     1    t*(     -0.0000000680d0 ))))
        f2 =     357.5291091806d0 +
     1    t*(  35999.0502911389d0 +
     1    t*(     -0.0001536667d0 +
     1    t*(      0.0000000378d0 +
     1    t*(     -0.0000000032d0 ))))
        f3 =      93.2720906200d0 +
     1    t*( 483202.0174577222d0 +
     1    t*(     -0.0035420000d0 +
     1    t*(     -0.0000002881d0 +
     1    t*(      0.0000000012d0 ))))
        f4 =     297.8501954694d0 +
     1    t*( 445267.1114469445d0 +
     1    t*(     -0.0017696111d0 +
     1    t*(      0.0000018314d0 +
     1    t*(     -0.0000000088d0 ))))
        f5 =     125.0445550100d0 +
     1    t*(  -1934.1362619722d0 +
     1    t*(      0.0020756111d0 +
     1    t*(      0.0000021394d0 +
     1    t*(     -0.0000000165d0 ))))
c
c  convert to Doodson (Darwin) variables
c
        d(1) = 360.d0*dayfr - f4
        d(2) = f3 + f5
	d(3) = d(2) - f4
        d(4) = d(2) - f1
        d(5) = -f5
        d(6) = d(3) - f2 
c
c   find frequencies of Delauney variables (in cycles/day), and from these
c    the same for the Doodson arguments
c 
        fd1 =  0.0362916471 + 0.0000000013*t
        fd2 =  0.0027377786
        fd3 =  0.0367481951 - 0.0000000005*t
        fd4 =  0.0338631920 - 0.0000000003*t
        fd5 = -0.0001470938 + 0.0000000003*t
        dd(1) = 1.d0 - fd4
        dd(2) = fd3 + fd5
        dd(3) = dd(2) - fd4
        dd(4) = dd(2) - fd1
        dd(5) = -fd5
        dd(6) = dd(3) - fd2 
      endif
c
c   end of intialization (likely to be called only once)
c
c  compute phase and frequency of the given tidal constituent
c
      freq=0.d0
      phase=0.d0
      do i = 1,6
         freq =   freq + idood(i)*dd(i)
         phase = phase + idood(i)*d(i)
      enddo
c adjust phases to fall in the positive range 0 to 360
      phase = dmod(phase,360.d0)
      if(phase.lt.0.d0) phase = phase + 360.d0
      return
      end
      subroutine toymd(it1,it2)
      dimension it1(*),it2(*)
c
c  converts times in it1 given as year and day to year-month-day in it2
c
      idn(m) = mod((367*(m-2-12*((m-14)/12)))/12+29,365)
      mon(jj,m) = (12*(jj-29-m))/367 + 2 + (jj-200)/169
      it2(1)=it1(1)
      it2(2)=mon(it1(2),leap(it1(1)))
      it2(3) = it1(2) - idn(it2(2)) - leap(it2(1))*((9+it2(2))/12)
      return
      end
      function leap(iy)
c
c returns 1 if year is a leap year, by Clavian (Gregorian) rule
c
      leap = 1 - (mod(iy,4)+3)/4
      if(mod(iy,100).eq.0.and.mod(iy,400).ne.0) leap=0
      return
      end
      function juldat(it)
      dimension it(*)
c
c  Julian Date from Gregorian date, Algorithm from p. 604, Explanatory
c    Supplement Amer Ephemeris & Nautical Almanac (cf Comm CACM, 11, 657 (1968)
c    and 15, 918 (1972)) Valid for all positive values of Julian Date
c
      juldat=(1461*(it(1)+4800+(it(2)-14)/12))/4
     1     + (367*(it(2)-2-12*((it(2)-14)/12)))/12
     2     - (3*((it(1)+4900+(it(2)-14)/12)/100))/4+it(3)-32075
      return
      end
      subroutine etutc(year,delta)
c
c   input year (and fraction thereof), 1700-2006
c   output delta = et - utc in seconds
c   tables from p 90, Explanatory Supplement to the A.E., Dr. R A Broucke,
c   JPL, and the leap.sec table in GAMIT
c   utc (and ut) is time usually used (eg in time signals)
c
c
      dimension d(142),tx(39),ty(39),st(24),si(24)
      data d/5.15,4.64,5.36,3.49,3.27,2.45,4.03,1.76,3.3,1.,2.42,.94,2.3
     11,2.27,-.22,.03,-.05,-.06,-.57,.03,-.47,.98,-.86,2.45,.22,.37,2.79
     2,1.2,3.52,1.17,2.67,3.06,2.66,2.97,3.28,3.31,3.33,3.23,3.6,3.52,
     34.27,2.68,2.75,2.67,1.94,1.39,1.66,.88,.33,-.17,-1.88,-3.43,-4.05,
     4-5.77,-7.06,-7.36,-7.67,-7.64,-7.93,-7.82,-8.35,-7.91,-8.03,-9.14,
     5-8.18,-7.88,-7.62,-7.17,-8.14,-7.59,-7.17,-7.94,-8.23,-7.88,-7.68,
     6-6.94,-6.89,-7.11,-5.87,-5.04,-3.9,-2.87,-.58,.71,1.8,3.08,4.63,5.
     786,7.21,8.58,10.5,12.1,12.49,14.41,15.59,15.81,17.52,19.01,18.39,1
     89.55,20.36,21.01,21.81,21.76,22.35,22.68,22.94,22.93,22.69,22.94,2
     93.2,23.31,23.63,23.47,23.68,23.62,23.53,23.59,23.99,23.8,24.2,
     124.99,24.97,25.72,26.21,26.37,26.89,27.68,28.13,28.94,29.42,29.66,
     130.29,30.96,31.09,31.59,32.06,31.82,32.69,33.05,33.16,33.59/
      data tx/61.5,
     2 62.    ,62.5     ,63.      ,63.5     ,64.      ,64.5     ,65.   ,
     3 65.5   ,66.      ,66.5     ,67.      ,67.5     ,68.      ,68.25 ,
     4 68.5   ,68.75    ,69.      ,69.25    ,69.5     ,69.75    ,70.   ,
     5 70.25  ,70.5     ,70.75    ,71.      ,71.085   ,71.162   ,71.247,
     6 71.329 ,71.414   ,71.496   ,71.581   ,71.666   ,71.748   ,71.833,
     7 71.915 ,71.999   ,72.0/
      data ty/33.59,
     2 34.032 ,34.235   ,34.441   ,34.644   ,34.95    ,35.286   ,35.725,
     3 36.16  ,36.498   ,36.968   ,37.444   ,37.913   ,38.39    ,38.526,
     4 38.76  ,39.      ,39.238   ,39.472   ,39.707   ,39.946   ,40.185,
     5 40.42  ,40.654   ,40.892   ,41.131   ,41.211   ,41.284   ,41.364,
     6 41.442 ,41.522   ,41.600   ,41.680   ,41.761   ,41.838   ,41.919,
     7 41.996 ,42.184   ,42.184/
c   update st,si, and nstep to allow for added leap seconds
      data st/1972.5,1973.0,1974.0,1975.0,1976.0,1977.0,1978.0,1979.0,
     1        1980.0,1981.5,1982.5,1983.5,1985.5,1988.0,1990.0,
     2        1991.0,1992.5,1993.5,1994.5,1996.0,1997.5,1999.0,
     3        2006.0,2008.0/
      data si/24*1./
      data nstep/24/
      if(year.gt.1972.) go to 6
      if(year.gt.1961.5) go to 2
      if(year.gt.1820.5) go to 1
c  for oldest epochs, use approximation
      if(year.ge.1785.) delta = 6.
      if(year.lt.1785.) delta = (year-1750.)/5.
      if(year.lt.1700.) delta = 0.
      return
c   for 1820.5 to 1961.5, data is spaced at yearly intervals
 1    n = year - 1819.5
      frac = year - (1819.5 + n)
      delta = (d(n+1) - d(n))*frac + d(n)
      return
c   for 1961.5 to 1972.0, interpolate between unequispaced data
 2    do 3 i = 1,38
      if(year-1900..eq.tx(i)) go to 4
      if(year-1900..lt.tx(i)) go to 5
 3    continue
 4    delta = ty(i)
      return
 5    delta=ty(i-1)+(ty(i)-ty(i-1))*((year-1900.-tx(i-1))/(tx(i)-tx(i-1)
     1))
      return
c   after 1972 et-utc has only step offsets. st is the array of step times,
c   si the step sizes (an added second is +1.)
 6    delta = 42.184
      do 7 i = 1,nstep
      if(year.ge.st(i)) delta = delta + si(i)
      if(year.lt.st(i)) return
 7    continue
      return
      end
