#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

m4_include([./seqtools_version.m4])

AC_PREREQ([2.50])
AC_INIT([SeqTools], VERSION_NUMBER, [gb10@sanger.ac.uk])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([blixemApp/blxmain.cpp])
AM_INIT_AUTOMAKE([1.9 foreign])
AC_CONFIG_HEADERS([config.h])

# Set up calling of Makefiles in toplevel directory and subdirectories
AC_CONFIG_FILES([Makefile
blixemApp/Makefile
dotterApp/Makefile
belvuApp/Makefile
seqtoolsUtils/Makefile
doc/Makefile
doc/User_doc/Makefile
doc/Design_notes/Makefile
doc/Design_notes/modules/Makefile
examples/Makefile
test/Makefile
test/data/Makefile
test/scripts/Makefile
test/scripts/automated/Makefile
test/scripts/automated/belvu/Makefile
test/scripts/automated/blixem/Makefile
test/scripts/automated/dotter/Makefile
test/scripts/manual/Makefile
test/scripts/manual/belvu/Makefile
test/scripts/manual/blixem/Makefile
test/scripts/manual/dotter/Makefile
])

# removed from the above: libpfetch/Makefile


# Check if this is a Git repository and, if so, include info about the current
# status in the version info. Otherwise use PACKAGE_VERSION (which will get
# set from VERSION_NUMBER in version.m4).
if test -d "$srcdir/../.git" 
then AC_DEFINE_UNQUOTED([SEQTOOLS_VERSION], ["`git describe --abbrev=1`"], ["SeqTools package version"])
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
AC_CHECK_HEADERS([arpa/inet.h netdb.h netinet/in.h stdlib.h string.h sys/socket.h unistd.h gtk.h sqlite3.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_PID_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_CHECK_FUNCS([dup2 floor gethostbyname memset socket sqrt strcasecmp strchr strcspn strerror strncasecmp strrchr strstr strtol uname])

# Check for dependencies required by all executables
PKG_CHECK_MODULES([DEPS], [glib-2.0 gtk+-2.0 >= 2.10])

# Check for dependencies required by sqlite code
PKG_CHECK_MODULES([DEPS_SQLITE3], [sqlite3], [HAVE_SQLITE3=1], [HAVE_SQLITE3=0])

# Check if the gbtools library exists as a subdirectory. This should be the case for
# the dist. It doesn't exist in the git repository though, so for development it
# either needs to be copied in or installed locally.
AM_CONDITIONAL([USE_GBTOOLS], [test -e "$srcdir/gbtools/configure.ac"])
AM_COND_IF([USE_GBTOOLS], [AC_CONFIG_SUBDIRS([gbtools])])

# Check for the curl library. This is only required by blixemh, so don't fail
# if we can't find it, but set the USE_LIBCURL to false. Note: I think we should 
# probably be using pkg-config to find the curl dependencies, but that doesn't work
# on my mac (perhaps because it's in /usr instead of /opt/local), so I just check for
# the header and set the dependencies manually in Makefile.am.
#PKG_CHECK_MODULES([LIBCURL_DEPS], [curl], [HAVE_LIBCURL=1], [HAVE_LIBCURL=0])
AC_CHECK_HEADERS([curl/curl.h], [HAVE_LIBCURL=1], [HAVE_LIBCURL=0])
AM_CONDITIONAL([USE_LIBCURL], [test "$HAVE_LIBCURL" -eq 1])

# Similar check for sqlite
AM_CONDITIONAL([USE_SQLITE3], [test "$HAVE_SQLITE3" -eq 1])

# Check for execinfo.h. This is used to provide backtraces but is not available
# on all systems.
AC_CHECK_HEADERS([execinfo.h])

AC_OUTPUT


echo \
"-------------------------------------------------

 ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

 Prefix: '${prefix}'
 Compiler: '${CC} ${CFLAGS} ${CPPFLAGS}'

 Now type 'make @<:@<target>@:>@'
   where the optional @<:@<target>@:>@ is:
     all       - build all binaries (default)
     install   - build and install everything
     <program> - build specific program
"

if test "$HAVE_LIBCURL" -eq 0; then
echo " Warning: libcurl not found; excluding blixemh
"
fi

if test "$HAVE_SQLITE3" -eq 0; then
echo " Warning: sqlite3 not found
"
fi

echo "-------------------------------------------------
"

