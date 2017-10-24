CC=gcc
FC=gfortran

OBJECTS=xdrfile.o xdrfile_trr.o

PROGS=gau2gro gro2gau ddip_nonpolar gromacs_link.sh
BINDIR=/path/to/bin


all: gau2gro gro2gau ddip_nonpolar


gau2gro:xdrfile.c xdrfile_trr.c gau2gro.c
	$(CC) $^ -o $@

gro2gau:gro2gau.f90
	$(FC) $^ -o $@

ddip_nonpolar:ddip_nonpolar.f90
	$(FC) $^ -o $@

install: $(PROGS) $(BINDIR)
	cp $(PROGS) $(BINDIR) -v

$(BINDIR):
	@echo ""
	@echo "To install the programs, please, specify the BINDIR as:"
	@echo ""
	@echo " make install BINDIR=/path/to/bin"
	@echo ""
	@echo "E.g."
	@echo " make install BINDIR=~/bin/"
	@echo ""
	@false

clean:
	rm gro2gau gau2gro ddip_nonpolar
