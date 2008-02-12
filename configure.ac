#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

# make sure $ACLOCAL_FLAGS are used during a rebuild.
AC_SUBST([ACLOCAL_AMFLAGS], ["-I $ac_macro_dir \${ACLOCAL_FLAGS}"])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_MKDIR_P
AC_LANG(C++)

# Checks for libraries.
AX_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_CHECK_HEADER(Utils/Clonable.h, ,)
AC_CHECK_HEADER(NumCalc/Matrix.h, ,)
AC_CHECK_HEADER(Seq/Sequence.h,  ,)

# Checks for library functions.
AX_EXT_HAVE_LIB(/lib /usr/lib /usr/local/lib, bpp-numcalc,main, ,[AC_MSG_ERROR([Library bpp-numcalc not found.])])
AX_EXT_HAVE_LIB(/lib /usr/lib /usr/local/lib, bpp-utils,main, ,[AC_MSG_ERROR([Library bpp-utils not found.])])
AX_EXT_HAVE_LIB(/lib /usr/lib /usr/local/lib, bpp-seq,main, ,[AC_MSG_ERROR([Library bpp-seq not found.])])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
