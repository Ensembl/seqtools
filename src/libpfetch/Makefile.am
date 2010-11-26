
noinst_LIBRARIES = libpfetch.a

libpfetch_a_SOURCES = libpfetch-utils.c libpfetch.c libpfetch.h libpfetch-cmarshal.c libcurlobject.c libcurlobject.h libpfetch_I.h libpfetch_P.h libcurlobject_I.h

libpfetch_a_marshal = libpfetch-cmarshal.list

libpfetch_a_DEPENDENCIES = $(nodist_libpfetch_a_SOURCES)

libpfetch_a_extra_sources = libpfetch-cmarshal.h \
libpfetch-cmarshal.c \
$(NULL)

nodist_libpfetch_a_SOURCES = $(libpfetch_a_extra_sources)

BUILT_SOURCES = $(libpfetch_a_extra_sources)
CLEANFILES = $(BUILT_SOURCES)
EXTRA_DIST = $(libpfetch_a_marshal)   

AM_CPPFLAGS = $(DEPS_CFLAGS) -I$(top_srcdir)

GENMARSHAL_COMMAND = $(GLIB_GENMARSHAL) --prefix=libpfetch_cmarshal

.list.h:
	$(GENMARSHAL_COMMAND) --header $< >$@

.list.c:
	$(GENMARSHAL_COMMAND) --body $< >$@

