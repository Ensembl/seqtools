
SUBDIRS = .

include $(top_srcdir)/Makefile.am.common

# there's a problem on ubuntu 14 with finding X11 so this is a hopefully temporary hack...
X_LIB = -lX11 -lm

blixem_CPPFLAGS = $(AM_CPPFLAGS)
BLX_LIBS = $(top_builddir)/seqtoolsUtils/libSeqtoolsUtils.a

# If gbtools is in a subdirectory, add it; otherwise look for a local installation
if USE_GBTOOLS
BLX_LIBS += $(top_builddir)/gbtools/.libs/libgbtools.a
else
BLX_LIBS += -lgbtools
endif

# the gtk deps etc. must go at the end so that gbtools can pick them up
BLX_LIBS += $(DEPS_LIBS) $(X_LIB)

# use sqlite for blixem if it's available
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

blixem_SOURCES = blxmain.cpp blxview.cpp blxFetch.cpp sequencecellrenderer.cpp bigpicture.cpp bigpicturegrid.cpp detailview.cpp detailviewtree.cpp blxwindow.cpp exonview.cpp coverageview.cpp blxdotter.cpp blxFetchDb.cpp blxview.hpp blxcontext.hpp blixem_.hpp detailview.hpp detailviewtree.hpp sequencecellrenderer.hpp bigpicture.hpp bigpicturegrid.hpp blxwindow.hpp exonview.hpp coverageview.hpp blxdotter.hpp 
blixem_LDADD = $(BLX_LIBS)

# Only compile the blixemh target if we have the libcurl library
if USE_LIBCURL
blixemh_SOURCES = $(blixem_SOURCES)
blixemh_CPPFLAGS = $(blixem_CPPFLAGS) -DPFETCH_HTML

# commented out for now....
#blixemh_LDADD = $(top_builddir)/libpfetch/libpfetch.a $(BLX_LIBS) -lcurl

# but I've commented it back in ....need to try a compile....
blixemh_LDADD =  $(BLX_LIBS) -lcurl

endif


# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/blixemApp/Makefile.in
