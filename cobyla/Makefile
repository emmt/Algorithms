.PHONY: clean test-orig test-f2c

srcdir=.
TMPDIR = /tmp
NAME = cobyla
NAME0 = ${NAME}-orig
NAME1 = ${NAME}-f2c
NAME2 = ${NAME}-fwrap
NAME3 = ${NAME}-revcom

CC = gcc
FC = gfortran

CFLAGS = -Wall -O2 -I. -DNDEBUG
FFLAGS = ${CFLAGS}

SINGLE_FLAGS = ${CFLAGS} -DSINGLE_PRECISION
DOUBLE_FLAGS = ${CFLAGS} -USINGLE_PRECISION

TESTING_SINGLE_FLAGS = ${SINGLE_FLAGS} -DTESTING
TESTING_DOUBLE_FLAGS = ${DOUBLE_FLAGS} -DTESTING

REVCOM = -DREVCOM

default: test

clean:
	rm -f *~ *.o *.diff
	rm -f ${NAME}-double ${NAME}-single
	rm -f ${NAME0}-double ${NAME0}-single
	rm -f ${NAME1}-double ${NAME1}-single
	rm -f ${NAME3}-double ${NAME3}-single


# Rule to produce output of tests.
%.out: %
	./$< > $@

# ### Original FORTRAN code.
#
# ${NAME0}-double: ${srcdir}/${NAME}.f
#  	${FC} ${FFLAGS} $< -o $@ -lm
#
# ${NAME0}-single: ${srcdir}/orig/${NAME}.f
#  	${FC} ${FFLAGS} $< -o $@ -lm
#

### Conversion with F2C.

${NAME1}-double.c: ${srcdir}/${NAME}.f
	rm -f ${TMPDIR}/${NAME}.c
	f2c -Aas -d${TMPDIR} $<
	cp -a ${TMPDIR}/${NAME}.c $@

${NAME1}-single.c: ${srcdir}/orig/${NAME}.f
	rm -f ${TMPDIR}/${NAME}.c
	f2c -Aas -d${TMPDIR} $<
	cp -a ${TMPDIR}/${NAME}.c $@

${NAME1}-double: ${srcdir}/${NAME1}-double.c
	${CC} ${CFLAGS} $< -o $@ -lf2c -lm

${NAME1}-single: ${srcdir}/${NAME1}-single.c
	${CC} ${CFLAGS} $< -o $@ -lf2c -lm

${NAME1}-double.out: ${NAME1}-double
	./$< > $@

${NAME1}-single.out: ${NAME1}-single
	./$< > $@

### Conversion and optimizations by hand.

${NAME}-double: ${srcdir}/${NAME}.c ${srcdir}/${NAME}.h
	${CC} ${TESTING_DOUBLE_FLAGS} $< -o $@ -lm

${NAME}-single: ${srcdir}/${NAME}.c ${srcdir}/${NAME}.h
	${CC} ${TESTING_SINGLE_FLAGS} $< -o $@ -lm

test: test-double test-single

test-double: ${NAME}-double.out ${NAME0}-double.out
	@if diff ${NAME0}-double.out ${NAME}-double.out >${NAME}-double.diff; then \
	    echo "test-double successful"; \
	else \
	    echo "test-double failed (see file ${NAME}-double.diff)"; \
	fi

test-single: ${NAME}-single.out ${NAME0}-single.out
	@if diff ${NAME0}-single.out ${NAME}-single.out >${NAME}-single.diff; then \
	    echo "test-single successful"; \
	else \
	    echo "test-single failed (see file ${NAME}-single.diff)"; \
	fi

cobyla.o: ${srcdir}/cobyla.c ${srcdir}/cobyla.h
	${CC} ${CFLAGS} -c $< -o $@

### Reverse communication version.

${NAME3}-double: ${srcdir}/${NAME3}.c
	${CC} ${TESTING_DOUBLE_FLAGS} ${REVCOM} $< -o $@ -lm

${NAME3}-single: ${srcdir}/${NAME3}.c
	${CC} ${TESTING_SINGLE_FLAGS} ${REVCOM} $< -o $@ -lm
