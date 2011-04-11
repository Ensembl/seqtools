# The main SeqTools makefile.am
# Process this file with automake to create Makefile.in.

SUBDIRS = seqtoolsUtils dotterApp blixemApp belvuApp libpfetch doc

bin_PROGRAMS = dotter blixem blixemh belvu

# blixem and blixemh have the same source files, most of which are compiled in libBlixemApp.a. However, there are a 
# few that need to be compiled with different flags/libraries, and those are specified here
blixem_SOURCES = $(top_srcdir)/blixemApp/blxmain.c $(top_srcdir)/blixemApp/blxview.c $(top_srcdir)/blixemApp/blxFetch.c
blixemh_SOURCES = $(blixem_SOURCES)
dotter_SOURCES = $(top_srcdir)/dotterApp/dotterMain.c
belvu_SOURCES = $(top_srcdir)/belvuApp/belvuMain.c

AM_CPPFLAGS = $(DEPS_CFLAGS) -I$(top_srcdir)
blixemh_CPPFLAGS = $(AM_CPPFLAGS) -DPFETCH_HTML

AM_LDFLAGS = -L$(top_builddir)/seqtoolsUtils
blixem_LDFLAGS = $(AM_LDFLAGS) -L$(top_builddir)/blixemApp
blixemh_LDFLAGS = $(blixem_LDFLAGS)
dotter_LDFLAGS = $(AM_LDFLAGS) -L$(top_builddir)/dotterApp
belvu_LDFLAGS = $(AM_LDFLAGS) -L$(top_builddir)/belvuApp

COMMON_LIBS = $(DEPS_LIBS) $(top_builddir)/seqtoolsUtils/libSeqtoolsUtils.a

blixem_LDADD = $(COMMON_LIBS) $(top_builddir)/blixemApp/libBlixemApp.a -lSeqtoolsUtils -lBlixemApp
blixemh_LDADD = $(COMMON_LIBS) $(top_builddir)/blixemApp/libBlixemApp.a $(top_builddir)/libpfetch/libpfetch.a -lcurl -lSeqtoolsUtils -lBlixemApp
dotter_LDADD = $(COMMON_LIBS) $(top_builddir)/dotterApp/libDotterApp.a -lSeqtoolsUtils -lDotterApp
belvu_LDADD = $(COMMON_LIBS) $(top_builddir)/belvuApp/libBelvuApp.a -lSeqtoolsUtils -lBelvuApp


