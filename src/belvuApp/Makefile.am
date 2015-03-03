
SUBDIRS = .

noinst_LIBRARIES = libBelvuApp.a

AM_CPPFLAGS = -g -Wall $(DEPS_CFLAGS)

libBelvuApp_a_SOURCES = belvu_.h belvu.c belvuWindow.h belvuWindow.c belvuAlignment.h belvuAlignment.c belvuTree.h belvuTree.c belvuConsPlot.h belvuConsPlot.c
libBelvuApp_a_LIBADD  = 
libBelvuApp_a_CFLAGS  = 

# there's a problem on ubuntu 14 with finding X11 so this is a hopefully temporary hack...
X_LIB = -lX11 -lm

AM_LDFLAGS = -L$(top_builddir)

bin_PROGRAMS = belvu
belvu_SOURCES = belvuMain.c
belvu_LDADD = $(top_builddir)/seqtoolsUtils/libSeqtoolsUtils.a libBelvuApp.a -lSeqtoolsUtils -lBelvuApp $(DEPS_LIBS) $(X_LIB)

# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/belvuApp/Makefile.in
