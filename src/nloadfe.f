c
c  Original date of RCS archived version is  Jul  3  2005
c
c
c     $Id: nloadf.f,v 1.2 2012/03/06 03:32:43 agnew Exp agnew $
c
c     $Log: nloadf.f,v $
c     Revision 1.2  2012/03/06 03:32:43  agnew
c     removed density (now done in ocmodl)
c
c     Revision 1.1  2011/09/23 22:39:02  agnew
c     Initial revision
c
c
c-------------------------------------
c  program nloadf
c
c   computes the load tide caused by the oceans, using a station-centered
c  grid and integrated green functions.
c    calls:  ignewt (for the newtonian green functions)
c            ocmodl (for the ocean load - calls ocstart to read model in)
c            getcl  (to get command-line parameters)
c            getgrf (to get Green functions)
c            lodout (to write results to standard output)
c
      character*80 stnam,mdfile,grfile
      character*50 mdnam
      character*1 modo,fingrd
      complex disp,grav,tilt,pothi,str,amp,cmp,cload
      dimension disp(3),tilt(2),str(3),cload(10)
c  matrix of green functions holds: vert disp, horiz disp, elastic gravity,
c                                   elastic tilt, e theta theta, e lam lam
      common /green/beg,end,spc,grfn(200,7)
      common /modpar/mdfile,dsym,icte(6),mdnam
      common /stloc/ct,st,rlam,ht
      common/coninf/dist,close,clat,clong,numnf,clnf,farnf
      data dr/.0174532935/,rd/57.29577951/,tpi/6.283185307/,minaz/150/,
     1 an/6.371e6/,cz/(0.,0.)/
      data close/0.0/,numnf/0/,nclose/0/
      an = 1./an
c  getcl读取nloadf函数命令行数据，分析是否有足够的参数变量，否则程序自动退出，
c  并显示错误。
      call getcl(stnam,rlat,rlam,ht,mdfile,grfile,modo,distm)
      ct = cos(dr*(90.-rlat))
      st = sin(dr*(90.-rlat))
      grav = cz
      pothi = cz
      tilt(1) = cz
      tilt(2) = cz
      do 1 i=1,3
      disp(i) = cz
 1    str(i) = cz
      distbr = endn
c  convolutions start here - read in Green functions, which come in
c  ntot sections; num is the current one, fingrd is whether or not to
c  use the detailed land-sea grid
c  读取格林函数文件，其中num是第j个采样间隔，ntot是M距离范围的划分，总共划分的份数
c  nf是Nj每个采样间隔中的采样数目，fingrd为水密度的选择，尤其是水陆相交处
c  这些参数还有待于进一步理解
 3    call getgrf(grfile,num,ntot,nf,fingrd)
      stp = dr*spc
c  if this is the first convolution interval see if (1) there is ocean at
c  the station location and (2) if it is below sea level. If both are true,
c  compute the gravity load (only) for a disk from distance 0 to the start
c  of the first integrated Green function
      if(num.eq.1) then
         call ocmodl(rlat,rlam,fingrd,amp)
c        write(*,*) amp, beg
	 if(amp.ne.cz.and.ht.lt.0) then
	    del=dr*beg/2
            call ignewt(del,2*del,g,t,pot)
c            elastic part getting by setting g=0.
            g = 0
            grav = grav + g*amp*tpi
	    close = 0.
	    clat = rlat
	    clong = rlam
	    nclose = 1
	 endif
      endif
c
c  loop through distance begins here
c
      do 13 ii=1,nf
      dist = beg + (ii-1.)*spc
      del = dr*dist
      cd = cos(del)
      sd = sin(del)
      call ignewt(del,stp,g,t,pot)
      g = 0
      g = g + grfn(ii,3)
      t = t + grfn(ii,4)
      pot = pot + grfn(ii,7)
      u = grfn(ii,1)
      v = grfn(ii,2)
      ett = grfn(ii,5)
      ell = grfn(ii,6)
c find the number of azimuth steps (scaled with sin(del), but kept above
c some minimum)
      naz = 360.*sd
      if(dist.lt.180.) naz = max0(naz,minaz)
      azstp = tpi/naz
      azstpd = 360./naz
      azimuth = azstpd/2
      caz = cos(azstp/2)
      saz = sin(azstp/2)
      saztp = sin(azstp)
      caztp = cos(azstp)
      stpfac = 2*saz
c
c  loop through azimuth begins here
c
      do 11 jj=1,naz
c do the spherical trigonometry to find the cell lat and long
      cb = cd*ct + sd*st*caz
      if(abs(cb).gt.1) cb = cb/abs(cb)
      sb = sqrt(1.-cb*cb)
      rlato = 90 - rd*acos(cb) 
      rlong = 0.
c  if near the north pole leave longitude at zero- otherwise find it
      if(sb.gt.1.e-3) then
         sb = 1./sb
         sg = sd*saz*sb
         cg = (st*cd - sd*ct*caz)*sb
         rlong = rlam + rd*atan2(sg,cg)
      endif
      call ocmodl(rlato,rlong,fingrd,amp)
      if(amp.ne.(0.,0.).and.modo.ne.'m') then
c if this is the first nonzero amplitude, save the distance and location
	 if(nclose.eq.0) then
	    close = dist
	    clat = rlato
	    clong = rlong
	    nclose = 1
         endif
         cmp = amp*stpfac
c  compute the loads.  disp is vert,ns,ew; tilt is ns,ew;
c                      strain is ns,ew,ne shear
         disp(1) = disp(1) + u*amp*azstp
         disp(2) = disp(2) + v*caz*cmp
         disp(3) = disp(3) + v*saz*cmp
         grav = grav + g*amp*azstp
	 pothi = pothi + pot*amp*azstp
         tilt(1) = tilt(1) + t*caz*cmp
         tilt(2) = tilt(2) + t*saz*cmp
         aa = .5*(caz*caz - saz*saz)
         bb = .5*azstp + aa*stpfac
         cc = .5*azstp - aa*stpfac
         str(1) = str(1) + (ett*bb + ell*cc)*amp
         str(2) = str(2) + (ell*bb + ett*cc)*amp
         str(3) = str(3) + (ett-ell)*saz*caz*saztp*caztp*amp
      elseif(amp.ne.(0.,0.).and.modo.eq.'m'.and.dist.le.distm) then
c  otherwise, output locations of cell corners, for later plotting
         call ldbxdr(azimuth,dist,azstpd,spc)
      endif
c  use recursion to step sin and cos of the azimuth
      azimuth = azimuth + azstpd
      xx = saz*caztp + caz*saztp
      caz = caz*caztp - saz*saztp
      saz = xx
 11   continue
 13   continue
      if(num.lt.ntot) go to 3
      if(modo.eq.'m') stop
c
c  done with convolution - add to disk file, and convert to amp and phase
c
      cload(1) = grav
      cload(10) = pothi
      cload(5) = tilt(1)
      cload(6) = tilt(2)
      do 15 i=1,3
      cload(i+1) = disp(i)
      cload(i+6) = str(i)
 15   continue
      call lodout(cload,stnam,rlat,rlam,ht,modo)
      stop
      end
