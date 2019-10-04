# This file is part of Autoconf 2.69.                  -*- Autoconf -*-
# Fortran languages support.
# Copyright (C) 2001, 2003-2012 Free Software Foundation, Inc.

# This file is part of Autoconf.  This program is free
# software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Under Section 7 of GPL version 3, you are granted additional
# permissions described in the Autoconf Configure Script Exception,
# version 3.0, as published by the Free Software Foundation.
#
# You should have received a copy of the GNU General Public License
# and a copy of the Autoconf Configure Script Exception along with
# this program; see the files COPYINGv3 and COPYING.EXCEPTION
# respectively.  If not, see <http://www.gnu.org/licenses/>.

# Written by David MacKenzie, with help from
# Franc,ois Pinard, Karl Berry, Richard Pixley, Ian Lance Taylor,
# Roland McGrath, Noah Friedman, david d zuhn, and many others.

# serial 0


## ------------------------------- ##
## Few compilers' characteristics. ##
## ------------------------------- ##

# AC_FC_PP_SRCEXT(EXT, [ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
# --------------------------------------------------------------
# Like AC_FC_SRCEXT, set the source-code extension used in Fortran (FC) tests
# to EXT (which defaults to f).  Also, look for any necessary additional
# FCFLAGS needed to allow this extension for preprocessed Fortran, and store
# them in the output variable FCFLAGS_<EXT> (e.g. FCFLAGS_f90 for EXT=f90).
# If successful, call ACTION-IF-SUCCESS.  If unable to compile preprocessed
# source code with EXT, call ACTION-IF-FAILURE, which defaults to failing with
# an error message.
#
# Some compilers allow preprocessing with either a Fortran preprocessor or
# with the C preprocessor (cpp).  Prefer the Fortran preprocessor, to deal
# correctly with continuation lines, `//' (not a comment), and preserve white
# space (for fixed form).
#
# (The flags for the current source-code extension, if any, are stored in
# $ac_fcflags_srcext and used automatically in subsequent autoconf tests.)
#
# For ordinary extensions like f90, etcetera, the modified FCFLAGS
# are needed for IBM's xlf*.  Also, for Intel's ifort compiler, the
# $FCFLAGS_<EXT> variable *must* go immediately before the source file on the
# command line, unlike other $FCFLAGS.  Ugh.
#
# Known extensions that enable preprocessing by default, and flags to force it:
# GNU: .F .F90 .F95 .F03 .F08, -cpp for most others,
#      -x f77-cpp-input for .f77 .F77; -x f95-cpp-input for gfortran < 4.4
# SGI: .F .F90, -ftpp or -cpp for .f .f90, -E write preproc to stdout
#      -macro_expand enable macro expansion everywhere (with -ftpp)
#      -P preproc only, save in .i, no #line's
# SUN: .F .F95, -fpp for others; -xpp={fpp,cpp} for preprocessor selection
#      -F preprocess only (save in lowercase extension)
# IBM: .F .F77 .F90 .F95 .F03, -qsuffix=cpp=EXT for extension .EXT to invoke cpp
#      -WF,-qnofpp -WF,-qfpp=comment:linecont:nocomment:nolinecont
#      -WF,-qlanglvl=classic or not -qnoescape (trigraph problems)
#      -d no #line in output, -qnoobject for preprocessing only (output in .f)
#      -q{no,}ppsuborigarg substitute original macro args before expansion
# HP:  .F, +cpp={yes|no|default} use cpp, -cpp, +cpp_keep save in .i/.i90
# PGI: -Mpreprocess
# Absoft: .F .FOR .F90 .F95, -cpp for others
# Cray: .F .F90 .FTN, -e Z for others; -F enable macro expansion everywhere
# Intel: .F .F90, -fpp for others, but except for .f and .f90, -Tf may also be
#        needed right before the source file name
# PathScale: .F .F90 .F95, -ftpp or -cpp for .f .f90 .f95
#         -macro_expand for expansion everywhere, -P for no #line in output
# Lahey: .F .FOR .F90 .F95, -Cpp
# NAGWare: .F .F90 .F95, .ff .ff90 .ff95 (new), -fpp for others
# Compaq/Tru64: .F .F90, -cpp, -P keep .i file, -P keep .i file
# f2c: .F, -cpp
# g95: .F .FOR .F90 .F95 .F03, -cpp -no-cpp, -E for stdout
AC_DEFUN([AC_FC_PP_SRCEXT],
[AC_LANG_PUSH(Fortran)dnl
AC_CACHE_CHECK([for Fortran flag to compile preprocessed .$1 files],
		ac_cv_fc_pp_srcext_$1,
[ac_ext=$1
ac_fcflags_pp_srcext_save=$ac_fcflags_srcext
ac_fcflags_srcext=
ac_cv_fc_pp_srcext_$1=unknown
case $ac_ext in #(
  [[fF]]77) ac_try=f77-cpp-input;; #(
  *) ac_try=f95-cpp-input;;
esac
for ac_flag in none -ftpp -fpp -Tf "-fpp -Tf" -xpp=fpp -Mpreprocess "-e Z" \
               -cpp -xpp=cpp -qsuffix=cpp=$1 "-x $ac_try" +cpp -Cpp; do
  test "x$ac_flag" != xnone && ac_fcflags_srcext="$ac_flag"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#if 0
#include <ac_nonexistent.h>
      choke me
#endif]])],
    [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#if 1
#include <ac_nonexistent.h>
      choke me
#endif]])],
       [],
       [ac_cv_fc_pp_srcext_$1=$ac_flag; break])])
done
rm -f conftest.$ac_objext conftest.$1
ac_fcflags_srcext=$ac_fcflags_pp_srcext_save
])
if test "x$ac_cv_fc_pp_srcext_$1" = xunknown; then
  m4_default([$3],
             [AC_MSG_ERROR([Fortran could not compile preprocessed .$1 files])])
else
  ac_fc_srcext=$1
  if test "x$ac_cv_fc_pp_srcext_$1" = xnone; then
    ac_fcflags_srcext=""
    FCFLAGS_[]$1[]=""
  else
    ac_fcflags_srcext=$ac_cv_fc_pp_srcext_$1
    FCFLAGS_[]$1[]=$ac_cv_fc_pp_srcext_$1
  fi
  AC_SUBST(FCFLAGS_[]$1)
  $2
fi
AC_LANG_POP(Fortran)dnl
])# AC_FC_PP_SRCEXT


# AC_FC_LINE_LENGTH([LENGTH], [ACTION-IF-SUCCESS],
#		    [ACTION-IF-FAILURE = FAILURE])
# ------------------------------------------------
# Look for a compiler flag to make the Fortran (FC) compiler accept long lines
# in the current (free- or fixed-format) source code, and adds it to FCFLAGS.
# The optional LENGTH may be 80, 132 (default), or `unlimited' for longer
# lines.  Note that line lengths above 254 columns are not portable, and some
# compilers (hello ifort) do not accept more than 132 columns at least for
# fixed format.  Call ACTION-IF-SUCCESS (defaults to nothing) if successful
# (i.e. can compile code using new extension) and ACTION-IF-FAILURE (defaults
# to failing with an error message) if not.  (Defined via DEFUN_ONCE to
# prevent flag from being added to FCFLAGS multiple times.)
# You should call AC_FC_FREEFORM or AC_FC_FIXEDFORM to set the desired format
# prior to using this macro.
#
# The known flags are:
# -f{free,fixed}-line-length-N with N 72, 80, 132, or 0 or none for none.
# -ffree-line-length-none: GNU gfortran
# -ffree-line-length-huge: g95 (also -ffixed-line-length-N as above)
#       -qfixed=132 80 72: IBM compiler (xlf)
#                -Mextend: Cray
#            -132 -80 -72: Intel compiler (ifort)
#                          Needs to come before -extend_source because ifort
#                          accepts that as well with an optional parameter and
#                          doesn't fail but only warns about unknown arguments.
#          -extend_source: SGI compiler
#  -W, -WNN (132, 80, 72): Absoft Fortran
#     +es, +extend_source: HP Fortran (254 in either form, default is 72 fixed,
#			   132 free)
#            -w, (-)-wide: Lahey/Fujitsu Fortran (255 cols in fixed form)
#                      -e: Sun Fortran compiler (132 characters)
#                    -132: NAGWare
#         -72, -f, -Wf,-f: f2c (a weak form of "free-form" and long lines).
#                  /XLine: Open Watcom
AC_DEFUN_ONCE([AC_FC_LINE_LENGTH],
[AC_LANG_PUSH([Fortran])dnl
m4_case(m4_default([$1], [132]),
  [unlimited], [ac_fc_line_len_string=unlimited
	               ac_fc_line_len=0
                       ac_fc_line_length_test='
      subroutine longer_than_132(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,'\
'arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19)'],
  [132],            [ac_fc_line_len=132
		       ac_fc_line_length_test='
      subroutine longer_than_80(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,'\
'arg10)'],
  [80],             [ac_fc_line_len=80
		       ac_fc_line_length_test='
      subroutine longer_than_72(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)'],
  [m4_warning([Invalid length argument `$1'])])
: ${ac_fc_line_len_string=$ac_fc_line_len}
AC_CACHE_CHECK(
[for Fortran flag needed to accept $ac_fc_line_len_string column source lines],
	       [ac_cv_fc_line_length],
[ac_cv_fc_line_length=unknown
ac_fc_line_length_FCFLAGS_save=$FCFLAGS
for ac_flag in none \
	       -ffree-line-length-none -ffixed-line-length-none \
	       -ffree-line-length-huge \
	       -ffree-line-length-$ac_fc_line_len \
	       -ffixed-line-length-$ac_fc_line_len \
	       -qfixed=$ac_fc_line_len -Mextend \
	       -$ac_fc_line_len -extend_source \
	       -W$ac_fc_line_len -W +extend_source +es -wide --wide -w -e \
               -f -Wf,-f -xline
do
  test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_line_length_FCFLAGS_save $ac_flag"
  AC_COMPILE_IFELSE([[$ac_fc_line_length_test
      end subroutine]],
		    [ac_cv_fc_line_length=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ac_fc_line_length_FCFLAGS_save
])
if test "x$ac_cv_fc_line_length" = xunknown; then
  m4_default([$3],
	     [AC_MSG_ERROR([Fortran does not accept long source lines], 77)])
else
  if test "x$ac_cv_fc_line_length" != xnone; then
    FCFLAGS="$FCFLAGS $ac_cv_fc_line_length"
  fi
  $2
fi
AC_LANG_POP([Fortran])dnl
])# AC_FC_LINE_LENGTH


# AC_FC_MODULE_EXTENSION
# ----------------------
# Find the Fortran 90 module file extension.  The module extension is stored
# in the variable FC_MODEXT and empty if it cannot be determined.  The result
# or "unknown" is cached in the cache variable ac_cv_fc_module_ext.
AC_DEFUN([AC_FC_MODULE_EXTENSION],
[AC_CACHE_CHECK([Fortran 90 module extension], [ac_cv_fc_module_ext],
[AC_LANG_PUSH(Fortran)
mkdir conftest.dir
cd conftest.dir
ac_cv_fc_module_ext=unknown
AC_COMPILE_IFELSE([[
      module conftest_module
      contains
      subroutine conftest_routine
      write(*,'(a)') 'gotcha!'
      end subroutine
      end module]],
  [ac_cv_fc_module_ext=`ls | sed -n 's,conftest_module\.,,p'`
   if test x$ac_cv_fc_module_ext = x; then
dnl Some F90 compilers use upper case characters for the module file name.
     ac_cv_fc_module_ext=`ls | sed -n 's,CONFTEST_MODULE\.,,p'`
   fi])
cd ..
rm -rf conftest.dir
AC_LANG_POP(Fortran)
])
FC_MODEXT=$ac_cv_fc_module_ext
if test "$FC_MODEXT" = unknown; then
  FC_MODEXT=
fi
AC_SUBST([FC_MODEXT])dnl
])


# AC_FC_MODULE_FLAG([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ---------------------------------------------------------------------
# Find a flag to include Fortran 90 modules from another directory.
# If successful, run ACTION-IF-SUCCESS (defaults to nothing), otherwise
# run ACTION-IF-FAILURE (defaults to failing with an error message).
# The module flag is cached in the ac_cv_fc_module_flag variable.
# It may contain significant trailing whitespace.
#
# Known flags:
# gfortran: -Idir, -I dir (-M dir, -Mdir (deprecated), -Jdir for writing)
# g95: -I dir (-fmod=dir for writing)
# SUN: -Mdir, -M dir (-moddir=dir for writing;
#                     -Idir for includes is also searched)
# HP: -Idir, -I dir (+moddir=dir for writing)
# IBM: -Idir (-qmoddir=dir for writing)
# Intel: -Idir -I dir (-mod dir for writing)
# Absoft: -pdir
# Lahey: -mod dir
# Cray: -module dir, -p dir (-J dir for writing)
#       -e m is needed to enable writing .mod files at all
# Compaq: -Idir
# NAGWare: -I dir
# PathScale: -I dir  (but -module dir is looked at first)
# Portland: -module dir (first -module also names dir for writing)
# Fujitsu: -Am -Idir (-Mdir for writing is searched first, then '.', then -I)
#                    (-Am indicates how module information is saved)
AC_DEFUN([AC_FC_MODULE_FLAG],[
AC_CACHE_CHECK([Fortran 90 module inclusion flag], [ac_cv_fc_module_flag],
[AC_LANG_PUSH([Fortran])
ac_cv_fc_module_flag=unknown
mkdir conftest.dir
cd conftest.dir
AC_COMPILE_IFELSE([[
      module conftest_module
      contains
      subroutine conftest_routine
      write(*,'(a)') 'gotcha!'
      end subroutine
      end module]],
  [cd ..
   ac_fc_module_flag_FCFLAGS_save=$FCFLAGS
   # Flag ordering is significant for gfortran and Sun.
   for ac_flag in -M -I '-I ' '-M ' -p '-mod ' '-module ' '-Am -I'; do
     # Add the flag twice to prevent matching an output flag.
     FCFLAGS="$ac_fc_module_flag_FCFLAGS_save ${ac_flag}conftest.dir ${ac_flag}conftest.dir"
     AC_COMPILE_IFELSE([[
      program main
      use conftest_module
      call conftest_routine
      end program]],
       [ac_cv_fc_module_flag="$ac_flag"])
     if test "$ac_cv_fc_module_flag" != unknown; then
       break
     fi
   done
   FCFLAGS=$ac_fc_module_flag_FCFLAGS_save
])
rm -rf conftest.dir
AC_LANG_POP([Fortran])
])
if test "$ac_cv_fc_module_flag" != unknown; then
  FC_MODINC=$ac_cv_fc_module_flag
  $1
else
  FC_MODINC=
  m4_default([$2],
    [AC_MSG_ERROR([unable to find compiler flag for module search path])])
fi
AC_SUBST([FC_MODINC])
# Ensure trailing whitespace is preserved in a Makefile.
AC_SUBST([ac_empty], [""])
AC_CONFIG_COMMANDS_PRE([case $FC_MODINC in #(
  *\ ) FC_MODINC=$FC_MODINC'${ac_empty}' ;;
esac])dnl
])
