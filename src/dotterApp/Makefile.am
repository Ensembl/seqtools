
SUBDIRS = .

noinst_LIBRARIES = libDotterApp.a

AM_CPPFLAGS = -g -Wall $(DEPS_CFLAGS) -I$(top_srcdir)

# If gbtools is in a subdirectory, add it to the include path
if USE_GBTOOLS
AM_CPPFLAGS += -I$(top_srcdir)/gbtools/src/include
endif

libDotterApp_a_SOURCES = greyramptool.cpp alignmenttool.cpp dotplot.cpp dotter.cpp dotterKarlin.cpp seqtoolsExonView.cpp dotter_.hpp dotter.hpp seqtoolsExonView.hpp
libDotterApp_a_LIBADD  = 
libDotterApp_a_CFLAGS  = 

# there's a problem on ubuntu 14 with finding X11 so this is a hopefully temporary hack...
X_LIB = -lX11 -lm

bin_PROGRAMS = dotter
dotter_SOURCES = dotterMain.cpp
dotter_LDADD = $(top_builddir)/dotterApp/libDotterApp.a $(top_builddir)/seqtoolsUtils/libSeqtoolsUtils.a $(DEPS_LIBS)  $(X_LIB)

# If gbtools is in a subdirectory, add it; otherwise look for a local installation
if USE_GBTOOLS
dotter_LDADD += $(top_builddir)/gbtools/libgbtools.a            
else
dotter_LDADD += -lgbtools
endif

# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/dotterApp/Makefile.in
