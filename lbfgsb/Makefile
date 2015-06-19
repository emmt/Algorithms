# where are the sources? (automatically filled in by configure script)
srcdir=.

# these values filled in by    yorick -batch make.i
Y_MAKEDIR=
Y_EXE=
Y_EXE_PKGS=
Y_EXE_HOME=
Y_EXE_SITE=
Y_HOME_PKG=

# ----------------------------------------------------- optimization flags

# options for make command line, e.g.-   make COPT=-g TGT=exe
COPT=$(COPT_DEFAULT)
TGT=$(DEFAULT_TGT)

FFLAGS=$(COPT)
FC=gfortran

# ------------------------------------------------ macros for this package

PKG_NAME=ylbgsb
PKG_I=${srcdir}/lbfgsb.i

OBJS=lbfgsb-wrap.o lbfgsb.o linpack.o blas.o timer.o

# change to give the executable a name other than yorick
PKG_EXENAME=yorick

# PKG_DEPLIBS=-Lsomedir -lsomelib   for dependencies of this package
PKG_DEPLIBS=-lgfortran
# set compiler (or rarely loader) flags specific to this package
PKG_CFLAGS=
PKG_LDFLAGS=

# list of additional package names you want in PKG_EXENAME
# (typically $(Y_EXE_PKGS) should be first here)
EXTRA_PKGS=$(Y_EXE_PKGS)

# list of additional files for clean
PKG_CLEAN=driver1.for driver2.for driver3.for \
          lbfgsb.for linpack.for blas.for timer.for

# autoload file for this package, if any
PKG_I_START=${srcdir}/lbfgsb-autoload.i
# non-pkg.i include files for this package, if any
PKG_I_EXTRA=

# -------------------------------- standard targets and rules (in Makepkg)

%.o: ${srcdir}/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

%.o: %.for
	$(FC) $(FFLAGS) $(PLUG_PIC) -o $@ -c $<

# We prefix all function/subroutine names with 'lb' to avoid conflicts
# (note that the suffix must be changed to avoid circular dependency).
%.for: ${srcdir}/lbfgsb-3.0/%.f
	sed -r <$< >$@ \
	    -e 's/(call|subroutine) +d(axpy|copy|scal|pofa|trsl)/\1 lb\2/g' \
	    -e 's/(^|[^a-z0-9_])d(dot|nrm2)([^a-z0-9_]|$$)/\1lb\2\3/g' \
	    -e 's/(call|subroutine) +(setu|main|cmpr|errc|hpso|lnsr|prn[123]|)lb/\1 lb\2/g' \
	    -e 's/(call|subroutine) +active/\1 lbactv/g' \
	    -e 's/(call|subroutine) +bmv/\1 lbbmv/g' \
	    -e 's/(call|subroutine) +cauchy/\1 lbcchy/g' \
	    -e 's/(call|subroutine) +formk/\1 lbfrmk/g' \
	    -e 's/(call|subroutine) +formt/\1 lbfrmt/g' \
	    -e 's/(call|subroutine) +freev/\1 lbfrev/g' \
	    -e 's/(call|subroutine) +matupd/\1 lbmupd/g' \
	    -e 's/(call|subroutine) +projgr/\1 lbprjg/g' \
	    -e 's/(call|subroutine) +subsm/\1 lbsubs/g' \
	    -e 's/(call|subroutine) +dcsrch/\1 lbsrch/g' \
	    -e 's/(call|subroutine) +dcstep/\1 lbcstp/g' \
	    -e 's/(call|subroutine) +timer/\1 lbtimr/g'

# set macros Makepkg uses in target and dependency names
# DLL_TARGETS, LIB_TARGETS, EXE_TARGETS
# are any additional targets (defined below) prerequisite to
# the plugin library, archive library, and executable, respectively
PKG_I_DEPS=$(PKG_I)
Y_DISTMAKE=distmake

include $(Y_MAKEDIR)/Make.cfg
include $(Y_MAKEDIR)/Makepkg
include $(Y_MAKEDIR)/Make$(TGT)

# override macros Makepkg sets for rules and other macros
# see comments in Y_HOME/Makepkg for a list of possibilities

# if this package built with mpy: 1. be sure mpy appears in EXTRA_PKGS,
# 2. set TGT=exe, and 3. uncomment following two lines
# Y_MAIN_O=$(Y_LIBEXE)/mpymain.o
# include $(Y_MAKEDIR)/Makempy

# configure script for this package may produce make macros:
# include output-makefile-from-package-configure

# reduce chance of yorick-1.5 corrupting this Makefile
MAKE_TEMPLATE = protect-against-1.5

# ------------------------------------- targets and rules for this package

# simple example:
#myfunc.o: myapi.h
# more complex example (also consider using PKG_CFLAGS above):
#myfunc.o: myapi.h myfunc.c
#	$(CC) $(CPPFLAGS) $(CFLAGS) -DMY_SWITCH -o $@ -c myfunc.c

# -------------------------------------------------------- end of Makefile
