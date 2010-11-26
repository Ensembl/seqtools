# The main SeqTools makefile.am
# Process this file with automake to create Makefile.in.

SUBDIRS = seqtoolsUtils dotterApp blixemApp libpfetch

bin_PROGRAMS = dotter blixem blixemh

# blixem and blixemh have the same source files, most of which are compiled in libBlixemApp.a. However, there are a 
# few that need to be compiled with different flags/libraries, and those are specified here
blixem_SOURCES = $(top_srcdir)/blixemApp/blxmain.c $(top_srcdir)/blixemApp/blxview.c $(top_srcdir)/blixemApp/blxFetch.c
blixemh_SOURCES = $(blixem_SOURCES)
dotter_SOURCES = $(top_srcdir)/dotterApp/dotterMain.c

AM_CPPFLAGS = $(DEPS_CFLAGS) -I$(top_srcdir)
blixemh_CPPFLAGS = $(AM_CPPFLAGS) -DPFETCH_HTML

COMMON_LIBS = $(DEPS_LIBS) $(top_builddir)/seqtoolsUtils/libSeqtoolsUtils.a

blixem_LDADD = $(COMMON_LIBS) $(top_builddir)/blixemApp/libBlixemApp.a
blixemh_LDADD = $(COMMON_LIBS) $(top_builddir)/blixemApp/libBlixemApp.a $(top_builddir)/libpfetch/libpfetch.a -lcurl
dotter_LDADD = $(COMMON_LIBS) $(top_builddir)/dotterApp/libDotterApp.a

