#
# The main SeqTools makefile.am
# Process this file with automake to create Makefile.in.
#

SUBDIRS = .

if USE_GBTOOLS
SUBDIRS += gbtools
endif

SUBDIRS += seqtoolsUtils dotterApp blixemApp belvuApp doc test examples

# these are included in the distribution, but are not installed
EXTRA_DIST = INSTALL.windows INSTALL.linux INSTALL.mac


# Extra files to remove for the maintainer-clean target.
# Note you cannot use this target to remove directories,
# hence the extra "local" target.
#
MAINTAINERCLEANFILES = \
$(top_srcdir)/Makefile.in \
$(top_srcdir)/aclocal.m4 \
$(top_srcdir)/config.h.in \
$(top_srcdir)/config.h.in~ \
$(top_srcdir)/configure \
$(top_srcdir)/seqtools_version.m4

GENERATED_FILES = $(top_srcdir)/config $(top_srcdir)/autom4te.cache

maintainer-clean-local:
	rm -rf $(GENERATED_FILES)
	mkdir $(top_srcdir)/config
	rm -rf $(top_srcdir)/build/*
	git checkout $(top_srcdir)/build

