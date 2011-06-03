#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

m4_include([./version.m4])

AC_PREREQ([2.50])
AC_INIT([SeqTools], VERSION_NUMBER, [gb10@sanger.ac.uk])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([blixemApp/blxmain.c])
AM_INIT_AUTOMAKE([1.9 foreign])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile seqtoolsUtils/Makefile libpfetch/Makefile blixemApp/Makefile dotterApp/Makefile doc/Makefile doc/User_doc/Makefile])

if test -d "$srcdir/../.git" 
then AC_DEFINE_UNQUOTED([SEQTOOLS_VERSION], ["`git describe --dirty`"], ["SeqTools package version"])
else AC_DEFINE([SEQTOOLS_VERSION], [PACKAGE_VERSION], ["SeqTools package version"])
fi

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_ARG_VAR(GLIB_GENMARSHAL, [The glib-genmarshal executable.])
AC_CHECK_PROG(GLIB_GENMARSHAL, glib-genmarshal, glib-genmarshal)

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([arpa/inet.h netdb.h netinet/in.h stdlib.h string.h sys/socket.h unistd.h gtk.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_PID_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_CHECK_FUNCS([dup2 floor gethostbyname memset socket sqrt strcasecmp strchr strcspn strerror strncasecmp strrchr strstr strtol uname])


PKG_CHECK_MODULES(DEPS, [glib-2.0 gtk+-2.0 >= 2.10])

AC_OUTPUT


echo \
"-------------------------------------------------

 ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

 Prefix: '${prefix}'.
 Compiler: '${CC} ${CFLAGS} ${CPPFLAGS}'

 Now type 'make @<:@<target>@:>@'
   where the optional <target> is:
     all      - build all binaries
     install  - install everything

--------------------------------------------------"
