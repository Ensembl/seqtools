
SUBDIRS = .

noinst_LIBRARIES = libBelvuApp.a

AM_CPPFLAGS = -g -Wall $(DEPS_CFLAGS) -I$(top_srcdir)

libBelvuApp_a_SOURCES = belvu_.hpp belvu.cpp belvuWindow.hpp belvuWindow.cpp belvuAlignment.hpp belvuAlignment.cpp belvuTree.hpp belvuTree.cpp belvuConsPlot.hpp belvuConsPlot.cpp
libBelvuApp_a_LIBADD  = 
libBelvuApp_a_CFLAGS  = 

# there's a problem on ubuntu 14 with finding X11 so this is a hopefully temporary hack...
X_LIB = -lX11 -lm

AM_LDFLAGS = -L$(top_builddir)/seqtoolsUtils -L$(top_builddir)/belvuApp

bin_PROGRAMS = belvu
belvu_SOURCES = belvuMain.cpp
belvu_LDADD = $(top_builddir)/seqtoolsUtils/libSeqtoolsUtils.a libBelvuApp.a -lSeqtoolsUtils -lBelvuApp $(DEPS_LIBS) $(X_LIB)

# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/belvuApp/Makefile.in
