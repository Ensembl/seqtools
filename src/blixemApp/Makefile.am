
SUBDIRS = .

noinst_LIBRARIES = libBlixemApp.a

libBlixemApp_a_SOURCES = sequencecellrenderer.c bigpicture.c bigpicturegrid.c detailview.c detailviewtree.c blxwindow.c exonview.c coverageview.c blxdotter.c blxFetchDb.c blxview.h blixem_.h detailview.h detailviewtree.h sequencecellrenderer.h bigpicture.h bigpicturegrid.h blxwindow.h exonview.h coverageview.h blxdotter.h 
libBlixemApp_a_LIBADD  = 
libBlixemApp_a_CFLAGS  = -I$(top_srcdir)

# there's a problem on ubuntu 14 with finding X11 so this is a hopefully temporary hack...
X_LIB = -lX11 -lm

# use sqlite for blixem if it's available
if USE_SQLITE3
AM_CPPFLAGS = -g -Wall $(DEPS_CFLAGS) $(DEPS_SQLITE3) -I$(top_srcdir)
BLX_LIBS = $(DEPS_LIBS) $(DEPS_SQLITE3_LIBS)
blixem_CPPFLAGS = $(AM_CPPFLAGS) -DSQLITE3
else
AM_CPPFLAGS = -g -Wall $(DEPS_CFLAGS) -I$(top_srcdir)
BLX_LIBS = $(DEPS_LIBS)
blixem_CPPFLAGS = $(AM_CPPFLAGS)
endif

if USE_LIBCURL
bin_PROGRAMS = blixem blixemh
else
bin_PROGRAMS = blixem
endif

AM_LDFLAGS = -L$(top_builddir)/seqtoolsUtils -L$(top_builddir)/blixemApp

blixem_SOURCES = blxmain.c blxview.c blxFetch.c
blixem_LDADD = $(top_builddir)/seqtoolsUtils/libSeqtoolsUtils.a libBlixemApp.a -lSeqtoolsUtils -lBlixemApp $(BLX_LIBS) $(X_LIB)

# Only compile the blixemh target if we have the libcurl library
if USE_LIBCURL
blixemh_SOURCES = $(blixem_SOURCES)
blixemh_CPPFLAGS = $(blixem_CPPFLAGS) -DPFETCH_HTML
blixemh_LDADD = $(top_builddir)/seqtoolsUtils/libSeqtoolsUtils.a libBlixemApp.a $(top_builddir)/libpfetch/libpfetch.a $(BLX_LIBS) -lcurl -lSeqtoolsUtils -lBlixemApp  $(X_LIB)
endif


# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/blixemApp/Makefile.in
