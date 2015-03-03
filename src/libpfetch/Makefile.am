
SUBDIRS = .

noinst_LIBRARIES = libpfetch.a

libpfetch_a_SOURCES = libpfetch-utils.c \
libpfetch.c \
libpfetch.h \
libpfetch_I.h \
libpfetch_P.h \
libcurlobject.c \
libcurlobject.h \
libcurlobject_I.h \
$(NULL)  

AM_CPPFLAGS = -g -Wall $(DEPS_CFLAGS)

libpfetch_a_CFLAGS = 

libpfetch_a_marshal = libpfetch-cmarshal.list

libpfetch_a_DEPENDENCIES = $(nodist_libpfetch_a_SOURCES)

libpfetch_a_extra_sources = libpfetch-cmarshal.h \
libpfetch-cmarshal.c \
$(NULL)

nodist_libpfetch_a_SOURCES = $(libpfetch_a_extra_sources)

BUILT_SOURCES = $(libpfetch_a_extra_sources)
CLEANFILES = $(BUILT_SOURCES)
EXTRA_DIST = $(libpfetch_a_marshal)   
#PFETCH_EXTRA_DIST = $(libpfetch_a_marshal)


GENMARSHAL_COMMAND = $(GLIB_GENMARSHAL) --prefix=libpfetch_cmarshal


# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/libpfetch/Makefile.in



.list.h:
	$(GENMARSHAL_COMMAND) --header $< >$@

.list.c:
	$(GENMARSHAL_COMMAND) --body $< >$@

