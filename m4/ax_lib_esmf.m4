# ===========================================================================
#      https://www.gnu.org/software/autoconf-archive/ax_lib_esmf.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_ESMF()
#
# DESCRIPTION
#
#   This macro checks the availability of the ESMF library.
#
#   The --with-esmf option is added, accepting one of the following values:
#
#     no    - do not check for the ESMF library 
#     yes   - do check for ESMF library in standard locations 
#     check - same as yes 
#     path  - absolute path for ESMF makefile fragment (esmf.mk).
#
#   If ESMF is successfully found and build support is provided for the
#   current language (_AC_LANG), this macro defines the following variables:
#
#   ESMF_VERSION
#   ESMF_XX
#   ESMF_XXFLAGS
#   ESMF_LDFLAGS
#   ESMF_LIBS
#
#   where XX is a short (uppercase) signature of _AC_LANG (e.g. CXX, FC).
#
#   If successful, the preprocessor macro HAVE_ESMF is defined and set to 1
#   and with_esmf is set to "yes". Otherwise, with_esmf="no".
#
#   Your configuration script can test $with_esmf to take any further
#   actions. 
#
#   To use the macro, one would add the following lines to "configure.ac"
#   before AC_OUTPUT:
#
#     dnl Check for ESMF support
#     AX_LIB_ESMF()
#
#   One could test $with_esmf for the outcome or display it as follows
#
#     echo "ESMF support:  $with_esmf"
#
# LICENSE
#
#   Copyright (c) 2019 University Corporation for Atmospheric Research,
#   Massachusetts Institute of Technology, Geophysical Fluid Dynamics Laboratory,
#   University of Michigan, National Centers for Environmental Prediction,
#   Los Alamos National Laboratory, Argonne National Laboratory,
#   NASA Goddard Space Flight Center.
#   All rights reserved.
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 1

AC_DEFUN([AX_LIB_ESMF], [

AC_REQUIRE([AC_PROG_GREP])
AC_REQUIRE([AC_PROG_SED])

dnl Assume ESMF is required if no argument
m4_define(_ax_lib_esmf_arg, [m4_normalize(m4_default([$1],[yes]))])

m4_case(_ax_lib_esmf_arg,
        [check],  [m4_define([_ax_lib_esmf_arg], [yes])],
        [yes|no], [],
        [])

dnl Add a default --with-esmf configuration option.
AC_ARG_WITH([esmf],
  AS_HELP_STRING(
    [--with-esmf=[yes/no/PATH]],
      [provide (PATH) or retrieve (yes) location of ESMF makefile fragment (default: _ax_lib_esmf_arg)]
  ),
  [],[]
)

if test "x$with_esmf" = x ; then
  with_esmf=_ax_lib_esmf_arg
fi

dnl Set absolute path of ESMF makefile fragment
AS_CASE([$with_esmf],
  [check|yes], [ax_lib_esmf_mk=${ESMFMKFILE}],
  [no],        [ax_lib_esmf_mk=""],
  [ax_lib_esmf_mk=$with_esmf]
)

dnl Set default values for defined variables
ESMF_FC=""
ESMF_FFLAGS=""
ESMF_FLIBS=""
ESMF_CXX=""
ESMF_CXXFLAGS=""
ESMF_LIBS=""
ESMF_LDFLAGS=""
ESMF_VERSION=""

dnl Set default return value
with_esmf=no

dnl Check if ESMF makefile fragment exists
if test "x$ax_lib_esmf_mk" != x ; then
  AC_MSG_CHECKING([for ESMF library])
  m4_warn([cross],
          [cannot check for file existence when cross compiling])dnl
  AC_CACHE_VAL([ax_cv_file_esmf_mk],
    [test "$cross_compiling" = yes &&
      AC_MSG_ERROR([cannot check for file existence when cross compiling])
      if test -r "$ax_lib_esmf_mk"; then
        AS_VAR_SET([ax_cv_file_esmf_mk], [yes])
      else
        AS_VAR_SET([ax_cv_file_esmf_mk], [no])
     fi])
  AS_VAR_IF([ax_cv_file_esmf_mk], [yes], 
            [AC_SUBST(ESMFMKFILE, [$ax_lib_esmf_mk])], 
            [ax_lib_esmf_mk=""
             AC_MSG_RESULT([no])])
fi

if test "x$ax_lib_esmf_mk" != x ; then

  dnl Check if ESMF supports current language
  AS_CASE("_AC_LANG",
    [C++],     [ax_lib_esmf_compiler=CXX],
    [Fortran], [ax_lib_esmf_compiler=F90],
    [AC_MSG_ERROR([ESMF does not support language: _AC_LANG])]
  )

  dnl Set output variables for current language
  ESMF_[]_AC_LANG_PREFIX[]=`$GREP ESMF_${ax_lib_esmf_compiler}COMPILER $ax_lib_esmf_mk | $SED 's/.*=//'`
  ESMF_[]_AC_LANG_PREFIX[]FLAGS=`$GREP ESMF_${ax_lib_esmf_compiler}COMPILEPATHS $ax_lib_esmf_mk | $SED 's/.*=//'`
  ESMF_LDFLAGS=`$GREP '^ESMF_${ax_lib_esmf_compiler}LINKPATHS' $ax_lib_esmf_mk | $SED 's/.*=//'`
  ESMF_LIBS=`$GREP ESMF_${ax_lib_esmf_compiler}ESMFLINKLIBS $ax_lib_esmf_mk | $SED 's/.*=//'`
  ESMF_VERSION=`$GREP ESMF_VERSION_STRING= $ax_lib_esmf_mk | $SED 's/.*=//'`

  AC_SUBST(ESMF_[]_AC_LANG_PREFIX[])
  AC_SUBST(ESMF_[]_AC_LANG_PREFIX[]FLAGS)
  AC_SUBST([ESMF_FLIBS])
  AC_SUBST([ESMF_LDFLAGS])
  AC_SUBST([ESMF_LIBS])
  AC_DEFINE([HAVE_ESMF], [1], [Defined if you have ESMF support])

  with_esmf=yes
  AC_MSG_RESULT([yes (version $[ESMF_VERSION])])

fi
]) # AX_LIB_ESMF
