# these values filled in by yorick -batch make.i
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

# ------------------------------------------------ macros for this package

PKG_NAME=ycobyla
PKG_I=cobyla.i

OBJS=ycobyla.o cobyla.o

# change to give the executable a name other than yorick
PKG_EXENAME=yorick

# PKG_DEPLIBS=-Lsomedir -lsomelib   for dependencies of this package
PKG_DEPLIBS=
# set compiler (or rarely loader) flags specific to this package
PKG_CFLAGS= -I.. -DINTEGER=long -USINGLE_PRECISION
PKG_LDFLAGS=

# list of additional package names you want in PKG_EXENAME
# (typically $(Y_EXE_PKGS) should be first here)
EXTRA_PKGS=$(Y_EXE_PKGS)

# list of additional files for clean
PKG_CLEAN=

# autoload file for this package, if any
PKG_I_START=
# non-pkg.i include files for this package, if any
PKG_I_EXTRA=

# -------------------------------- standard targets and rules (in Makepkg)

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

ycobyla.o: ycobyla.c ../cobyla.h
cobyla.o: ../cobyla.c ../cobyla.h

# -------------------------------------------------------- end of Makefile
