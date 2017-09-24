# SYNOPSIS
#
#   AX_BAM
#
# DESCRIPTION
#
#   Test for the BAM libraries (htslib or bamlib)
#
#   If no path to the installed bam library is given the macro searchs
#   under /usr, /usr/local, /opt and /opt/local and evaluates the
#   $BAM_ROOT environment variable.
#	Adapted from AX_BOOST_BASE
#
#   This macro calls:
#
#     AC_SUBST(BAM_CPPFLAGS) / AC_SUBST(BAM_LDFLAGS)
#
#   And sets:
#
#     HAVE_BAM, HAVE_HTSLIB
#
# LICENSE
#
#   Copyright (c) 2010 Cole Trapnell <cole@cs.umd.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_BAM],
[
AC_ARG_WITH([bam],
	AS_HELP_STRING([--with-bam@<:@=DIR@:>@], [use BAM libraries (default is yes) - it is possible to specify the root directory for BAM (optional)]),
	[
    if test "$withval" = "no"; then
		want_bam="no"
    elif test "$withval" = "yes"; then
        want_bam="yes"
        ac_bam_path=""
    else
	    want_bam="yes"
        ac_bam_path="$withval"
	fi
    ],
    [want_bam="yes"])


AC_ARG_WITH([bam-libdir],
        AS_HELP_STRING([--with-bam-libdir=LIB_DIR],
        [Force given directory for bam libraries. Note that this will overwrite library path detection, so use this parameter only if default library detection fails and you know exactly where your bam libraries are located.]),
        [
        if test -d $withval
        then
                ac_bam_lib_path="$withval"
        else
                AC_MSG_ERROR(--with-bam-libdir expected directory name)
        fi
        ],
        [ac_bam_lib_path=""]
)

if test "x$want_bam" = "xyes"; then
#	bam_lib_version_req=ifelse([$1], ,1.20.0,$1)
#	bam_lib_version_req_shorten=`expr $bam_lib_version_req : '\([[0-9]]*\.[[0-9]]*\)'`
#	bam_lib_version_req_major=`expr $bam_lib_version_req : '\([[0-9]]*\)'`
#	bam_lib_version_req_minor=`expr $bam_lib_version_req : '[[0-9]]*\.\([[0-9]]*\)'`
#	bam_lib_version_req_sub_minor=`expr $bam_lib_version_req : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
#	if test "x$bam_lib_version_req_sub_minor" = "x" ; then
#		bam_lib_version_req_sub_minor="0"
#    	fi
#	WANT_BAM_VERSION=`expr $bam_lib_version_req_major \* 100000 \+  $bam_lib_version_req_minor \* 100 \+ $bam_lib_version_req_sub_minor`
	AC_MSG_CHECKING(for htslib or bamlib)
	succeeded=no

	USE_HTSLIB=0

	dnl first we check the system location for bam libraries
	SEARCH_DIRS="/usr /usr/local /opt /opt/local"
	if test "$ac_bam_path" != ""; then
	    SEARCH_DIRS=$ac_bam_path
	fi
	for ac_bam_path_tmp in $SEARCH_DIRS ; do
		if test -d "$ac_bam_path_tmp/include/htslib" && test -r "$ac_bam_path_tmp/include/htslib"; then
			BAM_LDFLAGS="-L$ac_bam_path_tmp/lib"
			BAM_CPPFLAGS="-I$ac_bam_path_tmp/include"
			BAM_LIB_PARAM="-lhts"
			USE_HTSLIB=1
			AC_DEFINE(HAVE_HTSLIB,,[define if the BAM library is htslib])
			break;
		fi
		if test -d "$ac_bam_path_tmp/include/bam" && test -r "$ac_bam_path_tmp/include/bam"; then
			BAM_LDFLAGS="-L$ac_bam_path_tmp/lib"
			BAM_CPPFLAGS="-I$ac_bam_path_tmp/include"
			BAM_LIB_PARAM="-lbam"
			break;
		fi
	done

    dnl overwrite ld flags if we have required special directory with
    dnl --with-bam-libdir parameter
    if test "$ac_bam_lib_path" != ""; then
       BAM_LDFLAGS="-L$ac_bam_lib_path"
    fi

	CPPFLAGS_SAVED="$CPPFLAGS"
	CPPFLAGS="$CPPFLAGS $BAM_CPPFLAGS"
	export CPPFLAGS

	LDFLAGS_SAVED="$LDFLAGS"
	LDFLAGS="$LDFLAGS $BAM_LDFLAGS"
	export LDFLAGS

	if test "$USE_HTSLIB" == "0"; then

	   AC_LANG_PUSH(C++)
	        	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
					@%:@include <bam/bam.h>
					]], [[]])],[
       AC_MSG_RESULT(yes)
	   succeeded=yes
	   found_system=yes
       	],[
       	])
	   AC_LANG_POP([C++])

	else

	   AC_LANG_PUSH(C++)
	        	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
					@%:@include <htslib/hts.h>
					]], [[]])],[
       AC_MSG_RESULT(yes)
	   succeeded=yes
	   found_system=yes
       	],[
       	])
	   AC_LANG_POP([C++])

    fi

	if test "$succeeded" != "yes" ; then
		AC_MSG_ERROR([[We could not detect the bam libraries. Try installing htslib or libbam, or specifying a path to one with --with-bam=/some/path]])
	else
		BAM_LIB="$BAM_LIB_PARAM"
		AC_SUBST(BAM_CPPFLAGS)
		AC_SUBST(BAM_LDFLAGS)
		AC_SUBST(BAM_LIB)
		AC_DEFINE(HAVE_BAM,,[define if the BAM library is available])
	fi

        CPPFLAGS="$CPPFLAGS_SAVED"
       	LDFLAGS="$LDFLAGS_SAVED"
fi

])
