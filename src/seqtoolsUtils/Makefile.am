
SUBDIRS = .

noinst_LIBRARIES = libSeqtoolsUtils.a

AM_CPPFLAGS = -g -Wall $(DEPS_CFLAGS)

libSeqtoolsUtils_a_SOURCES = iupac.h version.h utilities.h utilities.c blxmsp.h blxmsp.c translate.c seqtoolsWebBrowser.c blxGff3Parser.h blxGff3Parser.c blxparser.h blxparser.c
libSeqtoolsUtils_a_LIBADD  = 
libSeqtoolsUtils_a_CFLAGS  = 

# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/seqtoolsUtils/Makefile.in
