# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([mpisort], [1.0])
AM_INIT_AUTOMAKE([1.15 foreign -Wall])

# Checks for programs.
AC_PROG_CC([mpicc])

# Checks for libraries.
AC_CHECK_LIB([z], [inflate])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
