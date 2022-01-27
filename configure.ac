# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([mpiSORT], [1.5.1])
AC_CONFIG_HEADERS([src/config.h])
AM_INIT_AUTOMAKE([1.13 foreign -Wall])
AC_PROG_RANLIB
AM_PROG_AR

# Checks for programs.
AC_PROG_CC([mpicc])

# Checks for libraries.
AC_CHECK_LIB([m], [log2])
AC_CHECK_LIB([z], [inflate])

# Checks for library functions.
AC_CHECK_FUNCS([malloc_trim])
AC_CHECK_FUNCS([memset]) 
AC_CHECK_FUNCS([munmap]) 
AC_CHECK_FUNCS([strchr]) 
AC_CHECK_FUNCS([strdup]) 
AC_CHECK_FUNCS([strndup])
AC_CHECK_FUNCS([strstr]) 
AC_CHECK_FUNCS([strtoull]) 
AC_CHECK_HEADERS([fcntl.h]) 
AC_CHECK_HEADERS([limits.h]) 
AC_CHECK_HEADERS([malloc.h]) 
AC_CHECK_HEADERS([stdint.h]) 
AC_CHECK_HEADERS([stdlib.h]) 
AC_CHECK_HEADERS([string.h]) 
AC_CHECK_HEADERS([unistd.h]) 
AC_CHECK_HEADER_STDBOOL 
AC_C_INLINE 
AC_FUNC_FSEEKO 
AC_FUNC_MALLOC 
AC_FUNC_MMAP 
AC_FUNC_REALLOC 
AC_TYPE_INT16_T 
AC_TYPE_INT64_T 
AC_TYPE_INT8_T 
AC_TYPE_OFF_T 
AC_TYPE_SIZE_T 
AC_TYPE_UINT16_T 
AC_TYPE_UINT32_T 
AC_TYPE_UINT8_T 

AC_ARG_WITH([libhts], [AS_HELP_STRING( [--with-libhts=DIR], [htslib installation location])], [with_libhts_path=$with_libhts] , [libhts_present="no"] )
AC_SUBST([LIBHTSDIR], $with_libhts_path)
AM_CONDITIONAL([HAVE_HTSLIB], [test "x$libhts_present" != "xno" ])
AM_COND_IF([HAVE_HTSLIB], [echo "compile with htslib"])

#Include flag for lzma, bz2, libcurl
#By Default libculr,bz2,lzma are disable
 
AC_ARG_ENABLE([curl], [AS_HELP_STRING([--enable-curl], [remove curl-based support])], [enable_curl=yes], [enable_curl=no])
AM_CONDITIONAL([ENABLE_CURL], [test "$enable_curl" = "yes"])
AC_SUBST(ENABLE_CURL)

#enable_lzma
AC_ARG_ENABLE([lzma], [AS_HELP_STRING([--enable-lzma], [remove support for LZMA-compressed CRAM files])], [enable_lzma=yes], [enable_lzma=no])
AM_CONDITIONAL([ENABLE_LZMA], [test "$enable_lzma" = "yes"])
AC_SUBST(ENABLE_LZMA)

#enable_bz2
AC_ARG_ENABLE([bz2], [AS_HELP_STRING([--enable-bz2], [remove support for BZ2-compressed CRAM files])], [enable_bz2=yes], [enable_bz2=no])
AM_CONDITIONAL([ENABLE_BZ2], [test "$enable_bz2" = "yes"])
AC_SUBST(ENABLE_BZ2)

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
