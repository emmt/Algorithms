F2CFLAGS = -Aas
F2CLIBS = -lf2c
CC = gcc
CPPFLAGS = -I.
CFLAGS = -Wall -O2

F2CSOURCES = addcon.c adjtol.c conres.c delcon.c eqcons.c fgcalc.c getd.c \
             getfes.c getmin.c initzu.c ktvec.c lsrch.c main.c minflc.c \
             minfun.c newcon.c satact.c sdegen.c sdirn.c stepbd.c zbfgs.c

F2COBJECTS = $(F2CSOURCES:.c=.o)

# The following filter command can be used to convert scientific exponential
# notation for floating-point number from C to FORTRAN.  Just use `cat` if
# no conversion is needed.
#
# FILTER = sed 's/\([0-9]\)[eE]\([-+]\?[0-9]\)/\1D\2/g'
FILTER = cat

default:
	echo $(F2COBJECTS)

%.c: ../fortran/%.f
	f2c $(F2CFLAGS) $<

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

tolmin-f2c: $(F2COBJECTS)
	$(CC) $^ -o $@ $(F2CLIBS) -lm

tolmin.o: tolmin.c

tolmin: tolmin.o
	$(CC) $^ -o $@ -lm

output.txt: tolmin
	./$< | $(FILTER) >$@

output-f2c.txt: tolmin-f2c
	./$< | $(FILTER) >$@

test: output.txt
	diff output-orig.txt $<

clean:
	$(RM) *~ *.o

veryclean: clean
	$(RM) $(F2CSOURCES)