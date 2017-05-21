
CC = gcc
CPPFLAGS = -I.
CFLAGS = -Wall -O2

default:

lincoa-f2c: lincoa-f2c.o lincoa.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@ -lf2c -lm

lincoa-f2c.o: lincoa-f2c.c lincoa.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

lincoa: lincoa.o lincoa.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@ -lm

lincoa.o: lincoa.c lincoa.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

lincoa-test: lincoa-test.c lincoa.o lincoa.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $< lincoa.o -o $@ -lm

lincoa.out: lincoa
	./$< | sed 's/\([0-9]\)[eE]\([-+]\?[0-9]\)/\1D\2/g' >$@

lincoa-f2c.out: lincoa-f2c
	./$< | sed 's/\([0-9]\)[eE]\([-+]\?[0-9]\)/\1D\2/g' >$@

test: output.txt
	diff "$<" "../fortran/output.txt"

clean:
	rm -f *~ *.o lincoa-f2c lincoa
