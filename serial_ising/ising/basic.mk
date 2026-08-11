# DO NOT CHANGE ANYTHING HERE!!! ##################################
# (unless you know *exactly* what you're doing...) 

# Seb RV, august 5 2024. I sure hope I know what Im doing. I removed isinglibconvert because
# we don't need it anymore.
# Making it always parallel

# objects and headers for ising_sa serial 
TOBJ =  move.o ising_sa.o savestate.o initialize.o brute_force.o\
        ../lam/distributions.o ../lam/error.o  ../lam/lsa.o ../lam/random.o

# these 2 lines are for parallel ising_sa 
#TPOBJ =  move.o ising_sa.o savestate.o initialize.o \
#				../lam/distributions.o  ../lam/lsa.o ../lam/error.o ../lam/random.o \
#				../lam/assign_groups.o

#calc_ave_error_bar
TEBOBJ = calc_ave_error_bar.o

#curve_fit
TCFOBJ = curve_fit.o

#printscore
TPSOBJ =  printscore.o ../lam/error.o

SOURCES = `ls *.c`

# Below here be build-targets...

all: $(ISINGEXECS)

calc_ave_error_bar: $(TEBOBJ)
	$(CC) -o calc_ave_error_bar $(CFLAGS) $(TEBOBJ) $(LIBS) 

curve_fit: $(TCFOBJ)
	$(CC) -o curve_fit $(CFLAGS) $(TCFOBJ) $(LIBS) 

initialize.o: initialize.c
	$(CC) -c $(CFLAGS) $(VFLAGS) initialize.c
# parallel stuff

ising_sa.o: ising_sa.c
	$(CC) -c -o ising_sa.o $(CFLAGS) $(VFLAGS) ising_sa.c

move.o: move.c
	$(CC) -c -o move.o $(CFLAGS) move.c

brute_force.o: brute_force.c
	$(CC) -c -o brute_force.o $(CFLAGS) brute_force.c

ising_sa: ising_sa.o $(TOBJ)
	$(CC) -o ising_sa $(CFLAGS) $(VFLAGS) $(TOBJ) \
		$(LIBS)

savestate.o: savestate.c
	$(CC) -c -o savestate.o $(CFLAGS) savestate.c

printscore.o: printscore.c
	$(CC) -c $(CFLAGS) $(VFLAGS)  printscore.c

printscore: $(TPSOBJ)
	$(CC) -o printscore $(CFLAGS) $(TPSOBJ) $(LIBS)

ising_sa:ising_sa.o

# ... and here be the cleanup and make deps targets

clean:
	rm -f *.o core*

Makefile: ${FRC}
	rm -f $@
	cp basic.mk $@
	echo "#Automatically generated dependencies list#" >> $@
	${CC} $(CFLAGS) -MM ${SOURCES} >> $@
	chmod -w $@


