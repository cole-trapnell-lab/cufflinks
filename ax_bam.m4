# SYNOPSIS
#
#   AX_BAM
#
# DESCRIPTION
#
#   Test for the BAM libraries of a particular version (or newer)
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
#     HAVE_BAM
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
	AC_MSG_CHECKING(for bamlib)
	succeeded=no

	dnl first we check the system location for bam libraries
	if test "$ac_bam_path" != ""; then
		BAM_LDFLAGS="-L$ac_bam_path/lib"
		BAM_CPPFLAGS="-I$ac_bam_path/include"
	else
		for ac_bam_path_tmp in /usr /usr/local /opt /opt/local ; do
			if test -d "$ac_bam_path_tmp/include/bam" && test -r "$ac_bam_path_tmp/include/bam"; then
				BAM_LDFLAGS="-L$ac_bam_path_tmp/lib"
				BAM_CPPFLAGS="-I$ac_bam_path_tmp/include"
				break;
			fi
		done
	fi

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

	AC_LANG_PUSH(C++)
     	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
	@%:@include <bam/bam.h>
	]], [[
	]])],[
        AC_MSG_RESULT(yes)
	succeeded=yes
	found_system=yes
       	],[
       	])
	AC_LANG_POP([C++])

	dnl if we found no bam with system layout we search for bam libraries
	dnl built and installed without the --layout=system option or for a staged(not installed) version
	if test "x$succeeded" != "xyes"; then
		_version=0
		if test "$ac_bam_path" != ""; then
			if test -d "$ac_bam_path" && test -r "$ac_bam_path"; then
				for i in `ls -d $ac_bam_path/include/bam-* 2>/dev/null`; do
					_version_tmp=`echo $i | sed "s#$ac_bam_path##" | sed 's/\/include\/bam-//' | sed 's/_/./'`
					V_CHECK=`expr $_version_tmp \> $_version`
					if test "$V_CHECK" = "1" ; then
						_version=$_version_tmp
					fi
					VERSION_UNDERSCORE=`echo $_version | sed 's/\./_/'`
					BAM_CPPFLAGS="-I$ac_bam_path/include/bam-$VERSION_UNDERSCORE"
				done
			fi
		else
			for ac_bam_path in /usr /usr/local /opt /opt/local ; do
				if test -d "$ac_bam_path" && test -r "$ac_bam_path"; then
					for i in `ls -d $ac_bam_path/include/bam-* 2>/dev/null`; do
						_version_tmp=`echo $i | sed "s#$ac_bam_path##" | sed 's/\/include\/bam-//' | sed 's/_/./'`
						V_CHECK=`expr $_version_tmp \> $_version`
						if test "$V_CHECK" = "1" ; then
							_version=$_version_tmp
	               					best_path=$ac_bam_path
						fi
					done
				fi
			done

			VERSION_UNDERSCORE=`echo $_version | sed 's/\./_/'`
			BAM_CPPFLAGS="-I$best_path/include/bam-$VERSION_UNDERSCORE"
            if test "$ac_bam_lib_path" = ""
            then
               BAM_LDFLAGS="-L$best_path/lib"
            fi

	    		if test "x$BAM_ROOT" != "x"; then
				if test -d "$BAM_ROOT" && test -r "$BAM_ROOT" && test -d "$BAM_ROOT/stage/lib" && test -r "$BAM_ROOT/stage/lib"; then
					version_dir=`expr //$BAM_ROOT : '.*/\(.*\)'`
					stage_version=`echo $version_dir | sed 's/bam_//' | sed 's/_/./g'`
			        	stage_version_shorten=`expr $stage_version : '\([[0-9]]*\.[[0-9]]*\)'`
					V_CHECK=`expr $stage_version_shorten \>\= $_version`
                    if test "$V_CHECK" = "1" -a "$ac_bam_lib_path" = "" ; then
						AC_MSG_NOTICE(We will use a staged bam library from $BAM_ROOT)
						BAM_CPPFLAGS="-I$BAM_ROOT"
						BAM_LDFLAGS="-L$BAM_ROOT/stage/lib"
					fi
				fi
	    		fi
		fi

		CPPFLAGS="$CPPFLAGS $BAM_CPPFLAGS"
		export CPPFLAGS
		LDFLAGS="$LDFLAGS $BAM_LDFLAGS"
		export LDFLAGS

		AC_LANG_PUSH(C++)
	     	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
		@%:@include <bam/version.hpp>
		]], [[
		]])],[
        	AC_MSG_RESULT(yes)
		succeeded=yes
		found_system=yes
       		],[
	       	])
		AC_LANG_POP([C++])
	fi

	if test "$succeeded" != "yes" ; then
		if test "$_version" = "0" ; then
			AC_MSG_ERROR([[We could not detect the bam libraries (version $bam_lib_version_req_shorten or higher). If you have a staged bam library (still not installed) please specify \$BAM_ROOT in your environment and do not give a PATH to --with-bam option.  If you are sure you have bam installed, then check your version number looking in <bam/version.hpp>. See http://randspringer.de/bam for more documentation.]])
		else
			AC_MSG_NOTICE([Your bam libraries seem too old (version $_version).])
		fi
	else
		BAM_LIB="-lbam"
		AC_SUBST(BAM_CPPFLAGS)
		AC_SUBST(BAM_LDFLAGS)
		AC_SUBST(BAM_LIB)
		AC_DEFINE(HAVE_BAM,,[define if the BAM library is available])
	fi

        CPPFLAGS="$CPPFLAGS_SAVED"
       	LDFLAGS="$LDFLAGS_SAVED"
fi

])
