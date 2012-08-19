# SYNOPSIS
#
#   AX_EIGEN
#
# DESCRIPTION
#
#   Test for the EIGEN libraries of a particular version (or newer)
#
#   If no path to the installed eigen library is given the macro searchs
#   under /usr, /usr/local, /opt and /opt/local and evaluates the
#   $EIGEN_ROOT environment variable. 
#	Adapted from AX_BOOST_BASE
#
#   This macro calls:
#
#     AC_SUBST(EIGEN_CPPFLAGS) / AC_SUBST(EIGEN_LDFLAGS)
#
#   And sets:
#
#     HAVE_EIGEN
#
# LICENSE
#
#   Copyright (c) 2010 Cole Trapnell <cole@cs.umd.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_EIGEN],
[
AC_ARG_WITH([eigen],
AS_HELP_STRING([--with-eigen@<:@=DIR@:>@], [use EIGEN libraries (default is yes) - it is possible to specify the root directory for EIGEN (optional)]),
[
if test "$withval" = "no"; then
want_eigen="no"
elif test "$withval" = "yes"; then
want_eigen="yes"
ac_eigen_path=""
else
want_eigen="yes"
ac_eigen_path="$withval"
fi
],
[want_eigen="yes"])

if test "x$want_eigen" = "xyes"; then
AC_MSG_CHECKING(for eigenlib)
succeeded=no

dnl first we check the system location for eigen libraries
if test "$ac_eigen_path" != ""; then
EIGEN_CPPFLAGS="-I$ac_eigen_path/include"
else
for ac_eigen_path_tmp in /usr /usr/local /opt /opt/local ; do
if test -d "$ac_eigen_path_tmp/include/eigen" && test -r "$ac_eigen_path_tmp/include/eigen"; then
EIGEN_CPPFLAGS="-I$ac_eigen_path_tmp/include"
break;
fi
done
fi

CPPFLAGS_SAVED="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $EIGEN_CPPFLAGS"
export EIGEN_CPPFLAGS

AC_LANG_PUSH(C++)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
@%:@include <Eigen/Dense>
]], [[
]])],[
AC_MSG_RESULT(yes)
succeeded=yes
found_system=yes
],[
])
AC_LANG_POP([C++])

CPPFLAGS="$CPPFLAGS $EIGEN_CPPFLAGS"
export CPPFLAGS
LDFLAGS="$LDFLAGS $EIGEN_LDFLAGS"
export LDFLAGS
export EIGEN_CPPFLAGS

if test "$succeeded" == "yes" ; then
AC_SUBST(EIGEN_CPPFLAGS)
AC_DEFINE(HAVE_EIGEN,,[define if the EIGEN library is available])
fi

CPPFLAGS="$CPPFLAGS_SAVED"
LDFLAGS="$LDFLAGS_SAVED"
fi

])
