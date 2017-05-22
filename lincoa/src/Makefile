.PHONY: clean test test-code test-f2c

#---- Customizations begin here -----------------------------------------------

# Compilers (only the C compiler is needed to build the library):
CC = gcc
FC = gfortran

# Compiler and linker flags:
#STRICT = -ansi -Wpedantic
CPPFLAGS = -I. -DNDEBUG
CFLAGS = -Wall -O2
FFLAGS = ${CFLAGS}
LDFLAGS =

# Command and flags to build a static library:
AR = ar
ARFLAGS = rv

# Set the following to true/false depending whether your terminal does support
# ANSI colors.
COLOR = true

#---- Customizations end here -------------------------------------------------

srcdir=.
TMPDIR = /tmp
NAME = lincoa
LIBNAME = lib${NAME}.a

default:

test: test-code

${NAME}-f2c: ${NAME}-f2c.o ${NAME}.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@ -lf2c -lm

${NAME}-f2c.o: ${NAME}-f2c.c ${NAME}.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

${NAME}.o: ${NAME}.c ${NAME}.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

${NAME}-test: ${NAME}-test.o ${NAME}.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@ -lm

${NAME}-test.o: ${NAME}.c ${NAME}.h
	$(CC) -DTESTING=1 $(CPPFLAGS) $(CFLAGS) -c $< -o $@

${NAME}-test.out: ${NAME}-test
	./$< | sed 's/\([0-9]\)[eE]\([-+]\?[0-9]\)/\1D\2/g' >$@

test-code: ${NAME}-test.out ${srcdir}/../fortran/output.txt
	@if diff $^ >${NAME}-test.diff; then \
	    msg="test of C version successful"; \
	    if ${COLOR}; then \
	        /bin/echo -e "\e[1;32m$${msg}\e[0;39;49m"; \
	    else \
	        echo "$${msg}"; \
	    fi; \
	else \
	    msg="test of C version failed (see file ${NAME}-test.diff)"; \
	    if ${COLOR}; then \
	        /bin/echo -e "\e[1;31m$${msg}\e[0;39;49m"; \
	    else \
	        echo "$${msg}"; \
	    fi; \
	fi

${NAME}-f2c.out: ${NAME}-f2c
	./$< | sed 's/\([0-9]\)[eE]\([-+]\?[0-9]\)/\1D\2/g' >$@

clean:
	rm -f *~ *.o ${NAME}-f2c ${NAME}
