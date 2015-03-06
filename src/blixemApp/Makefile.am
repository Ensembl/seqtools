
SUBDIRS = .

noinst_LIBRARIES = libBlixemApp.a

AM_CPPFLAGS = -g -Wall $(DEPS_CFLAGS) -I$(top_srcdir) -I$(top_srcdir)/gbtools/src/include
AM_LDFLAGS = -L$(top_builddir)/gbtools -L$(top_builddir)/seqtoolsUtils -L$(top_builddir)/blixemApp

libBlixemApp_a_SOURCES = sequencecellrenderer.cpp bigpicture.cpp bigpicturegrid.cpp detailview.cpp detailviewtree.cpp blxwindow.cpp exonview.cpp coverageview.cpp blxdotter.cpp blxFetchDb.cpp blxview.hpp blixem_.hpp detailview.hpp detailviewtree.hpp sequencecellrenderer.hpp bigpicture.hpp bigpicturegrid.hpp blxwindow.hpp exonview.hpp coverageview.hpp blxdotter.hpp 
libBlixemApp_a_LIBADD  = 
libBlixemApp_a_CFLAGS  = -I$(top_srcdir)

# there's a problem on ubuntu 14 with finding X11 so this is a hopefully temporary hack...
X_LIB = -lX11 -lm

# use sqlite for blixem if it's available
BLX_LIBS = $(DEPS_LIBS)
blixem_CPPFLAGS = $(AM_CPPFLAGS)

if USE_SQLITE3
AM_CPPFLAGS += $(DEPS_SQLITE3)
BLX_LIBS += $(DEPS_SQLITE3_LIBS)
blixem_CPPFLAGS += -DSQLITE3
endif

if USE_LIBCURL
bin_PROGRAMS = blixem blixemh
else
bin_PROGRAMS = blixem
endif

blixem_SOURCES = blxmain.cpp blxview.cpp blxFetch.cpp
blixem_LDADD = -lSeqtoolsUtils -lBlixemApp $(BLX_LIBS) $(X_LIB)

# Only compile the blixemh target if we have the libcurl library
if USE_LIBCURL
blixemh_SOURCES = $(blixem_SOURCES)
blixemh_CPPFLAGS = $(blixem_CPPFLAGS) -DPFETCH_HTML
blixemh_LDADD = $(top_builddir)/libpfetch/libpfetch.a $(BLX_LIBS) -lcurl -lSeqtoolsUtils -lBlixemApp  $(X_LIB)
endif


# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/blixemApp/Makefile.in
