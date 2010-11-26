
noinst_LIBRARIES = libDotterApp.a

libDotterApp_a_SOURCES = greyramptool.c alignmenttool.c dotplot.c dotter.c dotterKarlin.c dotterMain.c seqtoolsExonView.c dotter_.h dotter.h seqtoolsExonView.h

AM_CPPFLAGS = $(DEPS_CFLAGS) -I$(top_srcdir)

libDotterApp_a_LIBADD = $(top_builddir)/seqtoolsUtils/libSeqtoolsUtils.a
