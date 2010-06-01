 # AC_C_OPENMP
 # -----------
 # Check which options need to be passed to the C compiler to support OpenMP.
 # Set the OPENMP_CFLAGS variable to these options.
 # The options are necessary at compile time (so the #pragmas are understood)
 # and at link time (so the appropriate library is linked with).
 # This macro takes care to not produce redundant options if $CC $CFLAGS 
already
 # supports OpenMP. It also is careful to not pass options to compilers that
 # misinterpret them; for example, most compilers accept "-openmp" and create
 # an output file called 'penmp' rather than activating OpenMP support.
 AC_DEFUN([AC_C_OPENMP],
 [
   AC_MSG_CHECKING([whether to use OpenMP])
   AC_ARG_ENABLE(openmp,
     [AS_HELP_STRING([--disable-openmp], [do not use OpenMP])],
     [],
     [enable_openmp=yes])
   AC_MSG_RESULT([$enable_openmp])
   OPENMP_CFLAGS=
   if test "$enable_openmp" = yes; then
     AC_MSG_CHECKING([for $CC option to support OpenMP])
     AC_CACHE_VAL([ac_cv_prog_cc_openmp], [
       ac_cv_prog_cc_openmp=unsupported
       AC_LINK_IFELSE([
 #ifndef _OPENMP
  choke me
 #endif
 #include <omp.h>
 int main () { return omp_get_num_threads (); }
         ], [ac_cv_prog_cc_openmp="none needed"])
       if test "$ac_cv_prog_cc_openmp" = unsupported; then
         dnl Try these flags:
         dnl   GCC >= 4.2           -fopenmp
         dnl   SunPRO C             -xopenmp
         dnl   Intel C              -openmp
         dnl   SGI C, PGI C         -mp
         dnl   Tru64 Compaq C       -omp
         dnl   IBM C (AIX, Linux)   -qsmp=omp
         for brand in GCC SunPRO Intel SGI/PGI Compaq IBM; do
           case $brand in
             GCC)
               ac_conditional='defined __GNUC__'
               ac_option='-fopenmp' ;;
             SunPRO)
               ac_conditional='defined __SUNPRO_C || defined __SUNPRO_CC'
               ac_option='-xopenmp' ;;
             Intel)
               ac_conditional='defined __INTEL_COMPILER'
               ac_option='-openmp' ;;
             SGI/PGI)
               ac_conditional='defined __sgi || defined __PGI || defined 
__PGIC__'
               ac_option='-mp' ;;
             Compaq)
               ac_conditional='defined __DECC || defined __DECCXX'
               ac_option='-omp' ;;
             IBM)
               ac_conditional='defined __xlc__ || defined __xlC__'
               ac_option='-qsmp=omp' ;;
           esac
           if test $brand = GCC; then
             if test "$GCC" = yes; then
               ac_openmp_result=yes
             else
               ac_openmp_result=no
             fi
           else
             AC_EGREP_CPP([Brand], [
               #if $ac_conditional
                Brand
               #endif
               ], [ac_openmp_result=yes], [ac_openmp_result=no])
           fi
           if test $ac_openmp_result = yes; then
             ac_save_CFLAGS=$CFLAGS
             CFLAGS="$CFLAGS $ac_option"
             AC_LINK_IFELSE([
 #ifndef _OPENMP
  choke me
 #endif
 #include <omp.h>
 int main () { return omp_get_num_threads (); }
               ], [ac_cv_prog_cc_openmp=$ac_option])
             CFLAGS=$ac_save_CFLAGS
             break
           fi
         done
       fi
       ])
     AC_MSG_RESULT([$ac_cv_prog_cc_openmp])
     case $ac_cv_prog_cc_openmp in
       "none needed" | unsupported)
         OPENMP_CFLAGS= ;;
       *)
         OPENMP_CFLAGS=$ac_cv_prog_cc_openmp ;;
     esac
   fi
   AC_SUBST([OPENMP_CFLAGS])
 ])

