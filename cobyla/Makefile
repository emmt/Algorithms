.PHONY: clean test-orig test-f2c

srcdir=.
TMPDIR = /tmp
NAME = cobyla

CC = gcc
FC = gfortran

CFLAGS = -Wall -O2 -I. -DNDEBUG
FFLAGS = ${CFLAGS}

SINGLE = -DSINGLE_PRECISION
DOUBLE = -USINGLE_PRECISION

REVCOM = -DREVCOM

default: test

clean:
	rm -f *~ *.o *.diff
	rm -f ${NAME}-double ${NAME}-single
	rm -f ${NAME}-orig-double ${NAME}-orig-single
	rm -f ${NAME}-f2c-double ${NAME}-f2c-single
	rm -f ${NAME}-fwrap-double ${NAME}-fwrap-single
	rm -f ${NAME}-revcom-double ${NAME}-revcom-single


# Rule to produce output of tests.
%.out: %
	./$< > $@

# ### Original FORTRAN code.
#
# ${NAME}-orig-double: ${srcdir}/${NAME}.f
#  	${FC} ${FFLAGS} $< -o $@ -lm
#
# ${NAME}-orig-single: ${srcdir}/orig/${NAME}.f
#  	${FC} ${FFLAGS} $< -o $@ -lm
#

### Conversion with F2C.

${NAME}-f2c-double.c: ${srcdir}/${NAME}.f
	rm -f ${TMPDIR}/${NAME}.c
	f2c -Aas -d${TMPDIR} $<
	cp -a ${TMPDIR}/${NAME}.c $@

${NAME}-f2c-single.c: ${srcdir}/orig/${NAME}.f
	rm -f ${TMPDIR}/${NAME}.c
	f2c -Aas -d${TMPDIR} $<
	cp -a ${TMPDIR}/${NAME}.c $@

${NAME}-f2c-double: ${srcdir}/${NAME}-f2c-double.c
	${CC} ${CFLAGS} $< -o $@ -lf2c -lm

${NAME}-f2c-single: ${srcdir}/${NAME}-f2c-single.c
	${CC} ${CFLAGS} $< -o $@ -lf2c -lm

${NAME}-f2c-double.out: ${NAME}-f2c-double
	./$< > $@

${NAME}-f2c-single.out: ${NAME}-f2c-single
	./$< > $@

### Conversion and optimizations by hand.

${NAME}-double: ${srcdir}/${NAME}.c ${srcdir}/${NAME}.h
	${CC} ${CFLAGS} ${DOUBLE} -DTESTING $< -o $@ -lm

${NAME}-single: ${srcdir}/${NAME}.c ${srcdir}/${NAME}.h
	${CC} ${CFLAGS} ${SINGLE} -DTESTING $< -o $@ -lm

test: test-double test-single

test-double: ${NAME}-double.out ${NAME}-orig-double.out
	@if diff ${NAME}-orig-double.out ${NAME}-double.out \
	    >${NAME}-double.diff; then \
	    echo "test-double successful"; \
	else \
	    echo "test-double failed (see file ${NAME}-double.diff)"; \
	fi

test-single: ${NAME}-single.out ${NAME}-orig-single.out
	@if diff ${NAME}-orig-single.out ${NAME}-single.out \
	    >${NAME}-single.diff; then \
	    echo "test-single successful"; \
	else \
	    echo "test-single failed (see file ${NAME}-single.diff)"; \
	fi

cobyla.o: ${srcdir}/cobyla.c ${srcdir}/cobyla.h
	${CC} ${CFLAGS} -c $< -o $@

### Reverse communication version.

${NAME}-revcom-double: ${srcdir}/${NAME}.c ${srcdir}/${NAME}.h
	${CC} ${CFLAGS} ${DOUBLE} -DTESTING_REVCOM $< -o $@ -lm

${NAME}-revcom-single: ${srcdir}/${NAME}.c ${srcdir}/${NAME}.h
	${CC} ${CFLAGS} ${SINGLE} -DTESTING_REVCOM $< -o $@ -lm

test-revcom: test-revcom-double test-revcom-single

test-revcom-double: ${NAME}-revcom-double.out ${NAME}-orig-double.out
	@if diff ${NAME}-orig-double.out ${NAME}-revcom-double.out \
	    >${NAME}-revcom-double.diff; then \
	    echo "test-double successful"; \
	else \
	    echo "test-double failed (see file ${NAME}-revcom-double.diff)"; \
	fi

test-revcom-single: ${NAME}-revcom-single.out ${NAME}-orig-single.out
	@if diff ${NAME}-orig-single.out ${NAME}-revcom-single.out \
	    >${NAME}-revcom-single.diff; then \
	    echo "test-single successful"; \
	else \
	    echo "test-single failed (see file ${NAME}-revcom-single.diff)"; \
	fi

### FORTRAN wrapper.

${NAME}-fwrap-double: ${srcdir}/${NAME}.c ${srcdir}/${NAME}.h
	${CC} ${CFLAGS} ${DOUBLE} -DTESTING_FWRAP $< -o $@ -lm

${NAME}-fwrap-single: ${srcdir}/${NAME}.c ${srcdir}/${NAME}.h
	${CC} ${CFLAGS} ${SINGLE} -DTESTING_FWRAP $< -o $@ -lm

test-fwrap: test-fwrap-double test-fwrap-single

test-fwrap-double: ${NAME}-fwrap-double.out ${NAME}-orig-double.out
	@if diff ${NAME}-orig-double.out ${NAME}-fwrap-double.out \
	    >${NAME}-fwrap-double.diff; then \
	    echo "test-double successful"; \
	else \
	    echo "test-double failed (see file ${NAME}-fwrap-double.diff)"; \
	fi

test-fwrap-single: ${NAME}-fwrap-single.out ${NAME}-orig-single.out
	@if diff ${NAME}-orig-single.out ${NAME}-fwrap-single.out \
	    >${NAME}-fwrap-single.diff; then \
	    echo "test-single successful"; \
	else \
	    echo "test-single failed (see file ${NAME}-fwrap-single.diff)"; \
	fi
