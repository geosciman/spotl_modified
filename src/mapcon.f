c
c  Original date of RCS archived version is  Aug  3  2001
c
c
c     $Id: mapcon.f,v 1.2 2012/06/25 04:31:35 agnew Exp agnew $
c
c     $Log: mapcon.f,v $
c     Revision 1.2  2012/06/25 04:31:35  agnew
c     removed unused variables ind and indfi (only in data statement)
c
c     Revision 1.1  2011/09/23 22:39:01  agnew
c     Initial revision
c
c
c-------------------------------------
c  program mapcon; converts between binary and ASCII versions of the map files
c
      character*80 mdfil1,mdfil2
      character*1 cdir
      integer*2 index
      integer*4 map
c index is 360x720 integers, map 32x32 bits (64/degree). The following
c dimension and data statements would need to be changed for a different
c cell size and resolution.
      dimension index(259200),map(32)
      if(iargc().lt.1.or.iargc().gt.3.or.iargc().eq.2) then
	 write(6,100)
 100     format('usage: mapcon r file1 file2 > asciifile',/,
     1     'to convert binary files file1 and file2 to ASCII',/,
     2     '    (file1 is the bitmap, file2 the index)',/,
     2     'or     cat asciifile | mapcon f',/,
     3     'to convert ASCII to the binary files lndsea.bit and',
     4     ' lndsea.ind (names are hardwired)')
	 stop
      endif
      call getarg(1,cdir)
      if(cdir.ne.'r'.and.cdir.ne.'f') then
         write(6,100)
         stop
      endif
      if(cdir.eq.'r') then
         if(iargc().ne.3) then
	    write(6,100)
            stop
         endif
         call getarg(2,mdfil1)
         call getarg(3,mdfil2)
         llu = 2
	 luo = 6
         open(unit=llu,file=mdfil2,status='old',access=
     $       'sequential',form="unformatted")
         read(llu) index
         close(llu)
	 write(luo,120) index
 120     format(12i6)
         open(unit=llu,file=mdfil1,status='old',access=
     $       'sequential',form="unformatted")
 3       read(llu,end=7) map
	 write(luo,122) map
 122     format(6i12)
	 go to 3
 7       continue
	 close(llu)
      endif
      if(cdir.eq.'f') then
         if(iargc().ne.1) then
	    write(6,100)
            stop
         endif
	 ncells=0
         llu = 5
	 luo = 2
         open(unit=luo,file='lndsea.ind',status='new',access=
     $       'sequential',form="unformatted")
         read(llu,120) index
	 write(luo) index
         close(luo)
         open(unit=luo,file='lndsea.bit',status='new',access=
     $       'sequential',form="unformatted")
 13      read(llu,122,end=17) map
	 ncells = ncells+1
	 write(luo) map
	 go to 13
 17      continue
	 close(luo)
	 write(6,124) ncells,32*ncells
 124  format(i6,' mixed cells written, dimension mapful to ',i9)
      endif
      stop
      end
