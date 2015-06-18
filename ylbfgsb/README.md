C wrapper and Yorick plugin for L-BFGS-B
========================================

This package provides wrapper code to use L-BFGS-B algorithm in your C or
[Yorick](http://yorick.github.com/) code.

Directory [lbfgsb-3.0](./lbfgsb-3.0) contains original code of L-BFGS-B algorithm
(version 3.0) by R.  H. Byrd, P. Lu, J.  Nocedal and C. Zhu (see ref. [1]
and [2]).  This code can be downloaded at
http://users.iems.northwestern.edu/~nocedal/lbfgsb.html


Usage in C
----------

In addition to the FORTRAN code in directory [lbfgsb-3.0](./lbfgsb-3.0),
files [`lbfgsb.h`](./lbfgsb.h) and [`lbfgsb-wrap.c`](./lbfgsb-wrap.c) are
all you need.  You may have to edit the basic types defined at the
beginning of [`lbfgsb.h`](./lbfgsb.h) to match the settings of your FORTRAN
compiler.

Note that the build process renames all FORTRAN functions and subroutines
to be prefixed by `lb` to avoid conflicts.


Installation of the Yorick plugin
---------------------------------

In short, building and installing the Yorick plug-in can be as quick as:
```
    cd BUILD_DIR
    SRC_DIR/configure
    make
    make install
```
where `BUILD_DIR` is the build directory (at your convenience) and `SRC_DIR` is
the source directory of the plug-in code.  The build and source directories
can be the same in which case, call `./configure` to configure for building.

To use the plug-in, start the Yorick interpreter and type:
```
    #include "lbfgsb.i"
        help, lbfgsb_setup
```
More detailled explanations are given below.

0. You must have [Yorick](http://yorick.github.com/) installed on your
   machine.

1. Unpack the plug-in code somewhere.

2. Configure for compilation.  The are two possibilities:

   * For an **in-place build**, go to the source directory of the plug-in
     code and run the configuration script:
     ```
         cd SRC_DIR
         ./configure
     ```
     To see the configuration options, call:
     ```
         ./configure --help
     ```

   * To compile in a **different build directory**, say `BUILD_DIR`, create
     the build directory, go to the build directory, and run the
     configuration script:
     ```
         mkdir -p BUILD_DIR
         cd BUILD_DIR
         SRC_DIR/configure
     ```
     where `SRC_DIR` is the path to the source directory of the plug-in code.
     To see the configuration options, call:
     ```
         SRC_DIR/configure --help
     ```

3. Compile the code:
   ```
       make
   ```

4. Install the plug-in in Yorick directories:
   ```
       make install
   ```


References
----------
The L-BFGS-B algorithm is described in:

1. R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, *A limited memory algorithm
   for bound constrained optimization*, SIAM J. Scientific Computing 16
   (1995), no. 5, pp. 1190--1208.

2. C.  Zhu, R.H. Byrd, P. Lu, J.  Nocedal, *L-BFGS-B: a limited memory
   FORTRAN code for solving bound constrained optimization problems*,
   Tech. Report, NAM-11, EECS Department, Northwestern University, 1994.



License
-------

* The FORTRAN code of L-BFGS-B (version 3.0) has been released under the "New
  BSD License".

* The Yorick and C code have been released under the MIT "expat" License.
