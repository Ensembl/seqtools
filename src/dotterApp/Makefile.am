
SUBDIRS = .

noinst_LIBRARIES = libDotterApp.a

AM_CPPFLAGS = -g -Wall $(DEPS_CFLAGS) -I$(top_srcdir) -I$(top_srcdir)/gbtools/src/include
AM_LDFLAGS = -L$(top_builddir)/gbtools -L$(top_builddir)/seqtoolsUtils -L$(top_builddir)/dotterApp 

libDotterApp_a_SOURCES = greyramptool.cpp alignmenttool.cpp dotplot.cpp dotter.cpp dotterKarlin.cpp seqtoolsExonView.cpp dotter_.hpp dotter.hpp seqtoolsExonView.hpp
libDotterApp_a_LIBADD  = 
libDotterApp_a_CFLAGS  = 

# there's a problem on ubuntu 14 with finding X11 so this is a hopefully temporary hack...
X_LIB = -lX11 -lm

bin_PROGRAMS = dotter
dotter_SOURCES = dotterMain.cpp
dotter_LDADD = -lgbtools -lSeqtoolsUtils -lDotterApp $(DEPS_LIBS)  $(X_LIB)

# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/dotterApp/Makefile.in
