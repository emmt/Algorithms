.PHONY: clean test-orig test-f2c

#---- Customization begin here.

# Compilers (only the C compiler is needed to build the library):
CC = gcc
FC = gfortran

# Compiler and linker flags:
#STRICT = -ansi -Wpedantic
CFLAGS = ${STRICT} -Wall -O2 -I. -DNDEBUG
FFLAGS = ${CFLAGS}
LDFLAGS =

# Command and flags to build a static libarry:
AR = ar
ARFLAGS = rv

#---- Customization end here.

srcdir=.
TMPDIR = /tmp
NAME = bobyqa
LIBNAME = lib${NAME}.a

OBJS = ${NAME}.o

SINGLE = -DSINGLE_PRECISION
DOUBLE = -USINGLE_PRECISION

default: library

library: ${LIBNAME}

test: test-code

clean:
	rm -f *~ *.o ${LIBNAME} *.diff
	rm -f ${NAME}-f2c.c
	rm -f ${NAME}-orig ${NAME}-f2c ${NAME}-test

very-clean: clean
	rm -f *.out

${LIBNAME}: ${OBJS}
	${AR} ${ARFLAGS} $@ $^

${NAME}.o: ${srcdir}/${NAME}.c ${srcdir}/${NAME}.h
	${CC} ${CFLAGS} ${DOUBLE} -UTESTING -c $< -o $@


# Rule to produce output of tests (the exponential notation uses a 'D' in
# FORTRAN, an 'E' in C, so we have to filter this).

%.out: %
	./$< | sed 's/\([0-9]\)D\([-+]\?[0-9]\)/\1E\2/g' > $@


# ### Original FORTRAN code.

${NAME}-orig: ${srcdir}/${NAME}.f
	${FC} ${FFLAGS} $< -o $@ -lm


### Conversion with F2C.

${NAME}-f2c.c: ${srcdir}/${NAME}.f
	rm -f ${TMPDIR}/${NAME}.c
	f2c -Aas -d${TMPDIR} $<
	mv ${TMPDIR}/${NAME}.c $@

${NAME}-f2c: ${srcdir}/${NAME}-f2c.c
	${CC} ${CFLAGS} $< -o $@ -lf2c -lm

test-f2c: ${NAME}-f2c.out ${NAME}-orig.out
	@if diff $^ >${NAME}-f2c.diff; then \
	    echo "test of F2C version successful"; \
	else \
	    echo "test of F2C version failed (see file ${NAME}-f2c.diff)"; \
	fi


### Conversion and optimizations by hand.

${NAME}-test.o: ${srcdir}/${NAME}.c ${srcdir}/${NAME}.h
	${CC} ${CFLAGS} ${DOUBLE} -DTESTING -c $< -o $@

${NAME}-test: ${NAME}-test.o
	${CC} ${LDFLAGS} $< -o $@ -lm

test-code: ${NAME}-test.out ${NAME}-orig.out
	@if diff $^ >${NAME}-test.diff; then \
	    echo "test of C version successful"; \
	else \
	    echo "test of C version failed (see file ${NAME}-test.diff)"; \
	fi

