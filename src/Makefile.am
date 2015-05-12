#
# The main SeqTools makefile.am
# Process this file with automake to create Makefile.in.
#

SUBDIRS = .

if USE_GBTOOLS
SUBDIRS += gbtools
endif

if USE_LIBCURL
SUBDIRS += libpfetch
endif

SUBDIRS += seqtoolsUtils dotterApp blixemApp belvuApp doc test examples
# these are included in the distribution, but are not installed
EXTRA_DIST = INSTALL.windows 

