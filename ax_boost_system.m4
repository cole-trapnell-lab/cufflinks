# ===========================================================================
#            http://autoconf-archive.cryp.to/ax_boost_system.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BOOST_SYSTEM
#
# DESCRIPTION
#
#   Test for System library from the Boost C++ libraries. The macro requires
#   a preceding call to AX_BOOST_BASE. Further documentation is available at
#   <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_SYSTEM_LIB)
#     AC_SUBST(BOOST_SYSTEM_LIB)
#
#   And sets:
#
#     HAVE_BOOST_SYSTEM
#
# LICENSE
#
#   Copyright (c) 2009 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2009 Michael Tindal
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_BOOST_SYSTEM],
[
	AC_ARG_WITH([boost-system],
	AS_HELP_STRING([--with-boost-system@<:@=special-lib@:>@],
                   [use the System library from boost - it is possible to specify a certain library for the linker
                        e.g. --with-boost-system=boost_system-gcc-mt ]),
        [
        if test "$withval" = "no"; then
			want_boost="no"
        elif test "$withval" = "yes"; then
            want_boost="yes"
            ax_boost_user_system_lib=""
            ax_booth_user_system_lib=""
        else
		    want_boost="yes"
			echo "using $withval"
        	ax_boost_user_system_lib="$withval"
		fi
        ],
        [want_boost="yes"]
	)

	if test "x$want_boost" = "xyes"; then
        AC_REQUIRE([AC_PROG_CC])
        AC_REQUIRE([AC_CANONICAL_BUILD])
		CPPFLAGS_SAVED="$CPPFLAGS"
		CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
		export CPPFLAGS

		LDFLAGS_SAVED="$LDFLAGS"
		LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
		export LDFLAGS

        AC_CACHE_CHECK(whether the Boost::System library is available,
					   ax_cv_boost_system,
        [AC_LANG_PUSH([C++])
			 CXXFLAGS_SAVE=$CXXFLAGS
			 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[@%:@include <boost/system/system_error.hpp>]]),
                   ax_cv_boost_system=yes, ax_cv_boost_system=no)
			 CXXFLAGS=$CXXFLAGS_SAVE
             AC_LANG_POP([C++])
		])
		if test "x$ax_cv_boost_system" = "xyes"; then
            AC_SUBST(BOOST_CPPFLAGS)

			AC_DEFINE(HAVE_BOOST_SYSTEM,,[define if the Boost::System library is available])
            BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`

			LDFLAGS_SAVE=$LDFLAGS

            if test "x$ax_boost_user_system_lib" = "x"; then
                for libextension in `ls $BOOSTLIBDIR/libboost_system*.so* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_system.*\)\.so.*$;\1;'` `ls $BOOSTLIBDIR/libboost_system*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_system.*\)\.a*$;\1;'`; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_SYSTEM_LIB="-l$ax_lib"; AC_SUBST(BOOST_SYSTEM_LIB) link_system="yes"; break],
                                 [link_system="no"])
  				done
                if test "x$link_system" != "xyes"; then
                for libextension in `ls $BOOSTLIBDIR/boost_system*.dll* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_system.*\)\.dll.*$;\1;'` `ls $BOOSTLIBDIR/libboost_system*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_system.*\)\.a*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_SYSTEM_LIB="-l$ax_lib"; AC_SUBST(BOOST_SYSTEM_LIB) link_system="yes"; break],
                                 [link_system="no"])
  				done
                fi

            else
                BOOST_SYSTEM_LIB="$ax_boost_user_system_lib"; 
				AC_SUBST(BOOST_SYSTEM_LIB) 
				link_system="yes";
               

            fi
            
            if test "x$ax_boost_user_system_lib" = "x"; then
                for libextension in `ls $BOOSTLIBDIR/libboost_system*.so* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_system.*\)\.so.*$;\1;'` `ls $BOOSTLIBDIR/libboost_system*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_system.*\)\.a*$;\1;'`; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_SYSTEM_LIB="-l$ax_lib"; AC_SUBST(BOOST_SYSTEM_LIB) link_system="yes"; break],
                                 [link_system="no"])
  				done
                if test "x$link_system" != "xyes"; then
                for libextension in `ls $BOOSTLIBDIR/boost_system*.dll* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_system.*\)\.dll.*$;\1;'` `ls $BOOSTLIBDIR/libboost_system*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_system.*\)\.a*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_SYSTEM_LIB="-l$ax_lib"; AC_SUBST(BOOST_SYSTEM_LIB) link_system="yes"; break],
                                 [link_system="no"])
  				done
                fi

            else
                BOOST_SYSTEM_LIB="$ax_boost_user_system_lib"; 
				AC_SUBST(BOOST_SYSTEM_LIB) 
				link_system="yes";
               

            fi
		fi

		CPPFLAGS="$CPPFLAGS_SAVED"
    	LDFLAGS="$LDFLAGS_SAVED"
	fi
])
