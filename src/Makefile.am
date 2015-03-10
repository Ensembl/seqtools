#
# The main SeqTools makefile.am
# Process this file with automake to create Makefile.in.
#

# libpfetch is only used if we have libcurl
if USE_LIBCURL
SUBDIRS = libpfetch seqtoolsUtils dotterApp blixemApp belvuApp doc
else
SUBDIRS = seqtoolsUtils dotterApp blixemApp belvuApp doc
endif

# these are included in the distribution, but are not installed
EXTRA_DIST = INSTALL.windows 

