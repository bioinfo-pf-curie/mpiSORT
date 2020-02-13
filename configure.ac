# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([mpiSORT], [1.0])
AC_CONFIG_HEADERS([src/config.h])
AM_INIT_AUTOMAKE([1.13 foreign -Wall])

# Checks for programs.
AC_PROG_CC([mpicc])

# Checks for libraries.
AC_CHECK_LIB([m], [log2])
AC_CHECK_LIB([z], [inflate])

# Checks for library functions.
AC_CHECK_FUNCS([malloc_trim])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
