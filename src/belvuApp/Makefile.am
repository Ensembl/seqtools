
SUBDIRS = .

noinst_LIBRARIES = libBelvuApp.a

AM_CPPFLAGS = -g -Wall $(DEPS_CFLAGS) -I$(top_srcdir)

# If gbtools is in a subdirectory, add it to the include path
if USE_GBTOOLS
AM_CPPFLAGS += -I$(top_srcdir)/gbtools/src/include
endif

libBelvuApp_a_SOURCES = belvu_.hpp belvu.cpp belvuWindow.hpp belvuWindow.cpp belvuAlignment.hpp belvuAlignment.cpp belvuTree.hpp belvuTree.cpp belvuConsPlot.hpp belvuConsPlot.cpp
libBelvuApp_a_LIBADD  = 
libBelvuApp_a_CFLAGS  = 

# there's a problem on ubuntu 14 with finding X11 so this is a hopefully temporary hack...
X_LIB = -lX11 -lm

bin_PROGRAMS = belvu
belvu_SOURCES = belvuMain.cpp
belvu_LDADD = $(top_builddir)/belvuApp/libBelvuApp.a $(top_builddir)/seqtoolsUtils/libSeqtoolsUtils.a

# If gbtools is in a subdirectory, add it; otherwise look for a local installation
if USE_GBTOOLS
belvu_LDADD += $(top_builddir)/gbtools/libgbtools.a        
else
belvu_LDADD += -lgbtools
endif

# the gtk deps etc. must go at the end so that gbtools can pick them up
belvu_LDADD += $(DEPS_LIBS) $(X_LIB)

# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/belvuApp/Makefile.in
