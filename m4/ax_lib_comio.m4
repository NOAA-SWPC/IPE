# ===========================================================================
#      https://www.gnu.org/software/autoconf-archive/ax_lib_comio.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_COMIO()
#
# DESCRIPTION
#
#   This macro provides tests of the availability of the COMIO library.
#
#   The macro adds a --with-comio option accepting one of three values:
#
#     no   - do not check for the COMIO library.
#     yes  - do check for COMIO library in standard locations.
#     path - complete path to comio-config or COMIO installation prefix
#
#   It checks for COMIO compilation helper script (comio-config) and
#   calls the script to set the proper compiler options.
#
#   If COMIO is successfully found, this macro calls
#
#     AC_SUBST(COMIO_VERSION)
#     AC_SUBST(COMIO_FC)
#     AC_SUBST(COMIO_FFLAGS)
#     AC_SUBST(COMIO_FLIBS)
#     AC_SUBST(COMIO_LDFLAGS)
#     AC_DEFINE(HAVE_COMIO)
#
#   It also sets
#
#     with_comio="yes"
#
#   If COMIO is disabled or not found, this macros sets
#
#     with_comio="no"
#
#   Your configuration script can test $with_comio to take any further
#   actions. COMIO_F{FLAGS,LIBS} should be used when building Fortran
#   applications.
#
#   To use the macro, one would add the following lines to "configure.ac"
#   before AC_OUTPUT:
#
#     dnl Check for COMIO support
#     AX_LIB_COMIO()
#
#   One could test $with_comio for the outcome or display it as follows
#
#     echo "COMIO support:  $with_comio"
#
# LICENSE
#
#   Copyright (c) 2020 NOAA/ESRL/SWPC development team
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 1

AC_DEFUN([AX_LIB_COMIO], [

AC_REQUIRE([AC_FC_MODULE_EXTENSION])

dnl Check first argument is one of the recognized values.
dnl Fail eagerly if is incorrect as this simplifies case statements below.
if   test "m4_normalize(m4_default([$1],[]))" = ""        ; then
  dnl Assume 'yes' if no argument
  ax_lib_comio_enable="yes"
elif test "m4_normalize(m4_default([$1],[]))" = "yes"  ; then
  ax_lib_comio_enable="yes"
elif test "m4_normalize(m4_default([$1],[]))" = "no"; then
  ax_lib_comio_enable="no"
else
  AC_MSG_ERROR([
    Unrecognized value for AX[]_LIB_COMIO within configure.ac.
    If supplied, argument 1 must be either 'yes' or 'no'.
  ])
fi

dnl Add a default --with-comio configuration option.
AC_ARG_WITH([comio],
  AS_HELP_STRING(
    [--with-comio=[yes/no/PATH]],
    m4_case(m4_normalize([$1]),
            [full path to comio-config or base directory of COMIO installation])
  ),
  [if test "$withval" = "no"; then
     with_comio="no"
   elif test "$withval" = "yes"; then
     with_comio="yes"
   else
     with_comio="yes"
     if test -d "${withval}" ; then
       ax_lib_comio_prefix="${withval}"
       COMIO_CONFIG="${withval}/bin/comio-config"
     elif test -x "${withval}" ; then
       COMIO_CONFIG="${withval}"
     fi
   fi],
  [with_comio=$ax_lib_comio_enable]
)

dnl Set defaults to blank
COMIO_FC=""
COMIO_FFLAGS=""
COMIO_FLIBS=""
COMIO_LDFLAGS=""
COMIO_VERSION=""

dnl Look for COMIO compile helper script and use it to find compile options
if test "$with_comio" = "yes"; then
  if test -z "$COMIO_CONFIG"; then
    dnl Check for comio-config in the path
    AC_PATH_PROGS([COMIO_CONFIG], [comio-config], [])
  else
    dnl Check if provided comio-config works
    AC_PATH_PROGS([COMIO_CONFIG], [$COMIO_CONFIG], [])
  fi
  AC_MSG_CHECKING([for COMIO library])
  if test -z "$COMIO_CONFIG"; then
    AC_MSG_RESULT([no])
    AC_MSG_WARN([

Unable to locate COMIO compilation helper script 'comio-config'.
Please specify either --with-comio=<HELPER> as the full path to
COMIO helper script, or --with-comio=<LOCATION> as the full path
to COMIO base installation directory.
COMIO support is being disabled (equivalent to --with-comio=no).
])
    with_comio="no"
  else
    COMIO_VERSION=$(eval $COMIO_CONFIG --version)
    COMIO_FC=$(eval $COMIO_CONFIG --fc)
    COMIO_FFLAGS=$(eval $COMIO_CONFIG --fflags)
    for arg in $(eval $COMIO_CONFIG --flibs); do
      case "$arg" in
        -L*) COMIO_LDFLAGS="$COMIO_LDFLAGS $arg"
             ;;
          *) COMIO_FLIBS="$COMIO_FLIBS $arg"
             ;;
      esac
    done
    AC_MSG_RESULT([yes (version $[COMIO_VERSION])])
  fi

  dnl Check if the library works
  AC_MSG_CHECKING([whether COMIO library works])
  AC_LANG_PUSH([Fortran])
  ax_lib_comio_save_FCFLAGS="${FCFLAGS}"
  FCFLAGS="${FCFLAGS} $COMIO_FFLAGS"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[dnl
       use comio
    ])],
    [AC_MSG_RESULT(yes)],
    [dnl library cannot compile
     AC_MSG_RESULT(no)
     FCFLAGS=$ax_lib_comio_save_FCFLAGS
     AC_MSG_ERROR([Incompatible COMIO library])]
  )
  FCFLAGS=$ax_lib_comio_save_FCFLAGS
  AC_LANG_POP()

  AC_SUBST([COMIO_FC])
  AC_SUBST([COMIO_FFLAGS])
  AC_SUBST([COMIO_FLIBS])
  AC_SUBST([COMIO_LDFLAGS])
  AC_SUBST([COMIO_VERSION])
  AC_DEFINE([HAVE_COMIO], [1], [Defined if you have COMIO support])
fi
])
