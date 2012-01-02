m4_include([ax_boost_base.m4])
m4_include([ax_boost_thread.m4])
m4_include([ax_bam.m4])
m4_include([ax_check_zlib.m4])

define([svnversion], esyscmd([sh -c "svnversion|tr -d '\n'"]))dnl
AC_INIT([cufflinks], [1.3.0], [cole@cs.umd.edu])
AC_DEFINE(SVN_REVISION, "svnversion", [SVN Revision])

AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE

#AM_PATH_CPPUNIT(1.10.2)

AC_ARG_VAR(PYTHON, [python program])

 # Make sure CXXFLAGS is defined so that AC_PROG_CXX doesn't set it.
CXXFLAGS="$CXXFLAGS"
CFLAGS="$CFLAGS"

AC_LANG(C)

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL
AM_PATH_PYTHON([2.4])
AX_BOOST_BASE([1.38.0])
AX_BAM
AX_BOOST_THREAD
AX_CHECK_ZLIB()

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
#AC_FUNC_FORK
#AC_CHECK_FUNCS([floor memmove pow regcomp sqrt strchr strcspn strspn strstr])

# check the platform
AC_CANONICAL_HOST

# set CFLAGS and CXXFLAGS
user_CFLAGS=${CFLAGS}
generic_CFLAGS="-Wall -Wno-strict-aliasing -g -gdwarf-2 -Wuninitialized"
ext_CFLAGS=""
debug_CFLAGS=""
#echo "${host_cpu}-${host_os}"
case "${host_cpu}-${host_os}" in
  i*86-*linux*)
	    ext_CFLAGS="-march=i686";;
  i*86-darwin*)
	CFLAGS="-m64"
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM], [ext_CFLAGS="-arch x86_64"], []);;
  *)
    AC_MSG_CHECKING([if gcc accepts -m64])
    CFLAGS="-m64"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM], [ext_CFLAGS="-m64"; AC_MSG_RESULT([yes])],
					  [ext_CFLAGS="-D_FILE_OFFSET_BITS=64"; AC_MSG_RESULT([no])]);;
esac

AC_ARG_ENABLE(vectorize,    [  --enable-vectorize        Enable GCC auto-vectorization],
			  [ext_CFLAGS="${ext_CFLAGS} -ftree-vectorize -msse3 -ffast-math -ftree-vectorizer-verbose=99"], [])

AC_ARG_ENABLE(intel64,      [  --enable-intel64        optimize for Intel64 CPU such as Xeon and Core2],
			  [ext_CFLAGS="${ext_CFLAGS} -march=nocona"], [])

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
	                      [enable debugging info (default is no)])],
              [], [enable_debug=no])

			  
AC_ARG_ENABLE([optim],
              [AS_HELP_STRING([--enable-optim@<:@=0|1|2|3@:>@],
	                      [set optimization level (default is 3)])],
              [if test "x$enable_optim" = xyes; then enable_optim=3; fi],
              [enable_optim=3])

AS_IF([test "x$enable_optim" != xno], [ext_CFLAGS="$ext_CFLAGS -O$enable_optim"])	

AS_IF([test "x$enable_debug" = xyes],
      [debug_CFLAGS="-DDEBUG"],
      [debug_CFLAGS="-DNDEBUG"])
	  
AC_ARG_ENABLE(profiling,      [  --enable-profiling        enable profiling with google-perftools],
	  [ext_LDFLAGS="-lprofiler -ltcmalloc"], [])
	  
CFLAGS="${generic_CFLAGS} ${ext_CFLAGS} ${user_CFLAGS} ${debug_CFLAGS} ${OPENMP_CFLAGS}"
CXXFLAGS="$CFLAGS"
CXXFLAGS="$CXXFLAGS $BOOST_CPPFLAGS $BAM_CPPFLAGS"
LDFLAGS="$ext_LDFLAGS"

# Checks for structures/functions that can be used to determine system memory
AC_CHECK_MEMBERS([struct sysinfo.totalram], [], [], [#include <sys/sysinfo.h>])
AC_CHECK_DECLS([sysctl, CTL_HW, HW_PHYSMEM], [], [], [#include <sys/sysctl.h>])

AM_INIT_AUTOMAKE([-Wall -Werror tar-pax foreign])

AC_CONFIG_FILES([Makefile				 
				 src/Makefile])

AC_OUTPUT
# dump some configuration confirmations
echo \
"
-- ${PACKAGE_STRING} Configuration Results --
  C++ compiler:        ${CXX} ${CXXFLAGS} ${LDFLAGS}"

if test x"${GCC}" = x"yes" ; then
   gcc_version=`${CC} --version | head -n 1`
   echo "  GCC version:         ${gcc_version}"
else
   gcc_version=''
fi

echo \
"  Host System type:    ${host}
  Install prefix:      ${prefix}
  Install eprefix:     ${exec_prefix}

  See config.h for further configuration information.
  Email <${PACKAGE_BUGREPORT}> with questions and bug reports.
"

if test x"${PYTHON}" = x":"  ||  ! test -x "${PYTHON}"; then
    echo "WARNING! python was not found and is required to run some utility scripts"
    echo "  We recommend installing python and pointing configure to the installed location"
fi

