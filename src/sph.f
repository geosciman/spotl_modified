c
c  Original date of RCS archived version is  Mar  7  1996
c
c
c     $Id: sph.f,v 1.1 2011/09/23 22:39:02 agnew Exp agnew $
c
c     $Log: sph.f,v $
c     Revision 1.1  2011/09/23 22:39:02  agnew
c     Initial revision
c
c
c-------------------------------------
      subroutine sph(grlat,elong,ht)  
c   for a point at geographical north latitude grlat, east longitude elong      
c   (in degrees), and height ht (in meters), finds geocentric position   
c   and local g using formulae for a spheroid
c  
      common/obs/cthet,sthet,clong,slong,dvert,radn,g      
      data gn,ae,f,rm,dr/9.798277692,6378137.,.00335281,.00344978,
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
