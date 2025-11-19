# DO NOT CHANGE ANYTHING HERE!!! ##################################
# (unless you know *exactly* what you're doing...) 

# Seb RV, august 5 2024. I sure hope I know what Im doing. I removed tsplibconvert because
# we don't need it anymore.
# Making it always parallel

# objects and headers for tsp_sa serial 
TOBJ =  edge_wt.o move.o tsp_sa.o savestate.o initialize.o\
        ../lam/distributions.o ../lam/error.o  ../lam/lsa.o ../lam/random.o

# these 2 lines are for parallel tsp_sa-mpi 
TPOBJ = edge_wt.o move-mpi.o tsp_sa-mpi.o  savestate-mpi.o initialize.o \
				../lam/distributions.o  ../lam/lsa-mpi.o ../lam/error.o ../lam/random.o

#calc_ave_error_bar
TEBOBJ = calc_ave_error_bar.o

#curve_fit
TCFOBJ = curve_fit.o

#printscore
TPSOBJ =  printscore.o initialize.o edge_wt.o \
					../lam/error.o

SOURCES = `ls *.c`

# Below here be build-targets...

all: $(TSPEXECS)

calc_ave_error_bar: $(TEBOBJ)
	$(CC) -o calc_ave_error_bar $(CFLAGS) $(TEBOBJ) $(LIBS) 

curve_fit: $(TCFOBJ)
	$(CC) -o curve_fit $(CFLAGS) $(TCFOBJ) $(LIBS) 

initialize: intialize.c
	$(MPICC) -c $(CFLAGS) $(VFLAG) $(MPIFLAGS) initialize.c
# parallel stuff

tsp_sa.mpi: $(TPOBJ)
	$(MPICC) -o tsp_sa.mpi $(MPIFLAGS) $(TPOBJ) $(LIBS)

move-mpi.o: move.c
	$(MPICC) -c -o move-mpi.o $(MPIFLAGS) $(CFLAGS) move.c

tsp_sa-mpi.o: tsp_sa.c $(TPOBJ)
	$(MPICC) -c -o tsp_sa-mpi.o $(MPIFLAGS) $(CFLAGS) $(VFLAGS) tsp_sa.c


savestate-mpi.o: savestate.c
	$(MPICC) -c -o savestate-mpi.o $(MPIFLAGS) $(CFLAGS) savestate.c

printscore.o: printscore.c
	$(CC) -c $(CFLAGS) $(VFLAGS)  printscore.c

printscore: $(TPSOBJ)
	$(CC) -o printscore $(CFLAGS) $(TPSOBJ) $(LIBS)

tsp_sa:tsp_sa.mpi

# ... and here be the cleanup and make deps targets

clean:
	rm -f *.o core*

Makefile: ${FRC}
	rm -f $@
	cp basic.mk $@
	echo "#Automatically generated dependencies list#" >> $@
	${CC} $(CFLAGS) -M ${SOURCES} >> $@
	chmod -w $@


