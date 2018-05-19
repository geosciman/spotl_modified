# machine independent makefile
# 

all: ocst.a mapst.a nloadf oclook loadcomb harprp hartid ertid mapcheck modcon mapcon polymake

BIN = ../bin/

# uncomment the following lines for g77 linux/mac compilation
#  for Macs, get gcc by installing Xcode from apple; g77 from
#  http://hpc.sourceforge.net/
#
#FTN = g77
#FFLAGS = -O2  -Wunused -Wno-globals -fno-f2c -fno-automatic -fno-backslash -fno-globals
#CC = gcc
#CFLAGS = -O -c -archi386
# on some systems line -archi386 may need to be removed

# uncomment the following lines for gfortran compilation (on some Macs)
FTN = gfortran
FFLAGS = -O3 -Wuninitialized -fno-f2c -fno-automatic -fno-range-check -fno-backslash
CC = /usr/bin/gcc
CFLAGS = -c

# uncomment the following lines for Intel ifort compilation
#FTN = ifort
#FFLAGS = -O2 -f77rtl -fixed -no-cxxlib -no-automatic -assume nobscc
#CC = gcc
#CFLAGS = -c

# uncomment the following lines for hp-optimized compilation
#FTN = fort77
#FFLAGS = +O3 +U77 -K
#CC = cc
#CFLAGS = +O3 -c

# uncomment the following lines for sun-optimized compilation
#FTN = f77
#FFLAGS = -O3
#CC = cc
#CFLAGS = -O3 -c

ertid: ertid.f
	$(FTN) $(FFLAGS) ertid.f  -o $(BIN)ertid
	rm -f *.o

loadcomb: loadcomb.f etharm.f sph.f cteput.f
	$(FTN) $(FFLAGS) loadcomb.f etharm.f sph.f cteput.f -o $(BIN)loadcomb
	rm -f *.o

harprp: harprp.f
	$(FTN) $(FFLAGS) harprp.f -o $(BIN)harprp
	rm -f *.o

hartid: hartid.f admint.f recurs.f spline.f shells.f tdfrph.f
	$(FTN) $(FFLAGS) hartid.f admint.f recurs.f spline.f shells.f tdfrph.f -o $(BIN)hartid
	rm -f *.o

mapcheck: mapcheck.f $(MAPLIB)
	$(FTN)  $(FFLAGS) mapcheck.f $(MAPLIB) -o $(BIN)mapcheck
	rm -f *.o

mapcon: mapcon.f
	$(FTN) $(FFLAGS) mapcon.f  -o $(BIN)mapcon
	rm -f *.o

modcon: modcon.f
	$(FTN) $(FFLAGS) modcon.f  -o $(BIN)modcon
	rm -f *.o

nloadf: nloadf.f getcl.f ignewt.f getgrf.f lodout.f ldbxdr.f invspt.f seadens.f $(OCLIB) $(MAPLIB)
	$(FTN) $(FFLAGS) nloadf.f getcl.f ignewt.f getgrf.f lodout.f ldbxdr.f seadens.f invspt.f $(OCLIB) $(MAPLIB) -o $(BIN)nloadf
	rm -f *.o

oclook: oclook.f fndgde.f seadens.f $(OCLIB) $(MAPLIB)
	$(FTN) $(FFLAGS) oclook.f fndgde.f seadens.f $(OCLIB) $(MAPLIB) -o $(BIN)oclook
	rm -f *.o

polymake: polymake.f polydata.f
	$(FTN) $(FFLAGS) polymake.f polydata.f -o $(BIN)polymake 
	rm -f *.o

OCLIB = ocst.a
MAPLIB = mapst.a


$(OCLIB): \
	$(OCLIB)(celfnd.o)	$(OCLIB)(ocstart.o) $(OCLIB)(ocmodl.o)
	ranlib $(OCLIB)

$(MAPLIB): \
	$(MAPLIB)(polygon.o)	$(MAPLIB)(inpoly.o) $(MAPLIB)(lndsea3.o) \
	$(MAPLIB)(ispand.o)
	ranlib $(MAPLIB)

.f.a:
	$(FTN) -c $(FFLAGS) $<
	ar rv $@ $*.o
	rm -f $*.o

.c.a:
	$(CC) $(CFLAGS) $<
	ar rv $@ *.o
	rm -f $*.o
