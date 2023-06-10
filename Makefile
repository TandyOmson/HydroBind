#include $(AMBERHOME)/config.h


run_min_nmode : run_xmin run_xnmode

run_xmin:  xmin
	./xmin complex.prmtop complex.rst7 > xmin.out
	/bin/rm -f foo.rst7


run_xnmode:  xnmode
	./nmode complex.prmtop complex.rst7 > xnmode.out
	/bin/rm -f vecs

#---------------------------------------------------------------------------
xmin:  xmin.o
	$(CC) -o xmin xmin.o -L$(AMBERHOME)/lib  $(FLIBS) $(LM)

xnmode:  xnmode.o
	$(CC) -o xnmode xnmode.o -L$(AMBERHOME)/lib  $(FLIBS) $(LM)

#---------------------------------------------------------------------------
.c.o:
	$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -I$(INCDIR) -o $@ $<

tmd.MPI.o: tmd.c
	mpicc -DMPI -c $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -I$(INCDIR) -o $@ $<

clean:
	/bin/rm -f *.o nmode
