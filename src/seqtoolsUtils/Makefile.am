
SUBDIRS = .

AM_CPPFLAGS = -g -Wall $(DEPS_CFLAGS) -I$(top_srcdir)

# If gbtools is in a subdirectory, add it to the include path
if USE_GBTOOLS
AM_CPPFLAGS += -I$(top_srcdir)/gbtools/src/include
endif

if USE_LIBCURL
AM_CPPFLAGS += -DPFETCH_HTML
endif

noinst_LIBRARIES = libSeqtoolsUtils.a

libSeqtoolsUtils_a_SOURCES = iupac.hpp version.hpp utilities.hpp utilities.cpp blxmsp.hpp blxmsp.cpp translate.cpp seqtoolsWebBrowser.cpp blxGff3Parser.hpp blxGff3Parser.cpp blxparser.hpp blxparser.cpp seqtoolsFetch.hpp seqtoolsFetch.cpp
libSeqtoolsUtils_a_LIBADD  = 
libSeqtoolsUtils_a_CFLAGS  =

# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/seqtoolsUtils/Makefile.in
