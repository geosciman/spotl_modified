      subroutine lndsea(rlat,rlong,lnd)
c   for given north latitude and east longitude (in degrees) returns
c  the flag lnd = 1 if land, 2 if water (cf paul revere)
c    ***calls ibit
c
      save index,mapful
      double precision clat,clong,rlong2
      integer*2 index
      integer*4 map,mapful
      dimension index(259200),map(32),mapful(525536)
c     data nlat/360/,nlong/720/,nfine/32/
      data ifrst/0/
      if(ifrst.eq.0) then
	 ifrst=1
         llu = 2
         open(unit=llu,file='lndsea.ind',status='old',access=
     $       'sequential',form="unformatted")
         read(llu) index
         close(llu)
	 nbl=1
         open(unit=llu,file='lndsea.bit',status='old',access=
     $       'sequential',form="unformatted")
 3       read(llu,end=7) map
	 do 5 i=1,32
 5       mapful(32*(nbl-1)+i) = map(i)
	 nbl=nbl+1
	 go to 3
 7       continue
	 close(llu)
      endif
c special branches for the polar regions, to avoid possible indexing trouble
      if(rlat.ge.85.) then
	 lnd=2
	 return
      endif
      if(rlat.le.-87.) then
	 lnd=1
	 return
      endif
      rlong2 = rlong
      if(rlong.lt.0) rlong2 = rlong + 360.
      if(rlong.gt.360) rlong2 = rlong - 360.
      clat = 90-rlat
      ilat = int(2*clat) + 1
      ilong = int(2*rlong2) + 1
      ind = 720*(ilat-1) + ilong
      if(index(ind).lt.0) then
	 lnd=-index(ind)
	 return
      endif
c are in a cell that is mixed land and sea: need to read the bit
c compute the lat and long of the SW corner
      clat = 90 - 0.5*ilat
      clong = 0.5*(ilong-1)
      idlat=min(int(64*(rlat-clat)) + 1,32)
      idlong=min(int(64*(rlong2-clong)) + 1,32)
c idlat and idlong are the indices to the bit-mapped version (32x32 bits)
c of the cell; the min function is to take care of the case when rlat
c or rlong are exactly 0.5 degrees from clat and clong.
c Now convert these to word and bit
      nwd = idlong + 32*(index(ind)-1)
      nbit = idlat
      iv = ibit(mapful(nwd),nbit)
      if(iv.eq.1) lnd = 1
      if(iv.eq.0) lnd = 2
      return
      end
      function ibit(id,n)
c
c  function to return the n-th bit (counting to the right) of a long
c  integer id
c
c    ***calls system routine and - this performs a bitwise and on
c       its two arguments
c
      save ibp,ifl
      integer*4 id,ibp(32),it
      data ifl/0/
      if(ifl.eq.1) go to 3
c  on first call, set up an array of integers, each with exactly one bit
c  set at successively lower values. ibp(1) should have only its highest
c  bit set; that this is true should be verified for each machine
      ifl = 1
      ibp(32) = 1
      do 1 i=2,32
 1    ibp(33-i) = 2*ibp(34-i)
 3    it = ispand(id,ibp(n))
      ibit = 0
      if(it.ne.0) ibit = 1
      return
      end
