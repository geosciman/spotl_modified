c  program polymake
c
c  reads (from standard input) a set of polygon names preceded by + or -
c   and writes out a polygon file for use by nloadf
c
c    $Id: polymake.f,v 1.3 2011/11/28 01:09:17 agnew Exp $
c
c    $Log: polymake.f,v $
c    Revision 1.3  2011/11/28 01:09:17  agnew
c    modified common block to match new polydata
c
c    Revision 1.2  2011/11/19 01:32:48  agnew
c    dimension changed for number of coordinates
c
c    Revision 1.1  2011/11/18 16:57:11  agnew
c    Initial revision
c
c
c
      character*80 line
      character*1 code
      character*30 inname
      character*30 names
      dimension code(50),inname(50)
c  common block with data (in subroutine polydata)
      common/polys/coords(1516),ind(23),names(22),namlen(22)
      data num/22/
      nin=1
 3    read(5,100,end=7) line
 100  format(a80)
      code(nin)=line(1:1)
      inname(nin)=line(3:32)
c check for comments or bad code symbol
      if(code(nin).eq.'#') go to 3
      if(code(nin).ne.'+'.and.code(nin).ne.'-') then
        write(6,101) code(nin)
 101    format('Code symbol "',a1,'" is not one of + - or #')
        go to 3
      endif
c
c loop through names to check them
c
      ifnd=0
      do j=1,num
        m=namlen(j)
        if(inname(nin)(1:m).eq.names(j)(1:m)) then
          ifnd=1
          nin=nin+1
        endif
      enddo
      if(ifnd.eq.0) write(0,102) inname(nin)
 102  format('Polygon ',a30,' not in database')
      go to 3
 7    nin=nin-1
c
c  write the start of the file
c
      write(6,104)
 104  format('Polygon file produced by polymake')
      write(6,106) nin
 106  format(i3)
c
c loop through names and for each one write out the coordinates
c

      do i=1,nin
        write(6,108) inname(i)
 108    format(a30)
        do j=1,num
          m=namlen(j)
          if(inname(i)(1:m).eq.names(j)(1:m)) then
            write(6,110) ind(j+1)-ind(j)
 110        format(i3)
            write(6,112) code(i)
 112        format(a1)
            do k=2*ind(j)-1,2*ind(j+1)-3,2
              write(6,114) coords(k),coords(k+1)
 114          format(f9.3,f8.3)
            enddo
          endif
        enddo
      enddo
      stop
      end
