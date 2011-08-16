#
# The main SeqTools makefile.am
# Process this file with automake to create Makefile.in.
#

# Define some shortcuts to source-code subdirectories
UTILS = $(top_srcdir)/seqtoolsUtils
PFETCH = $(top_srcdir)/libpfetch
BLIXEM = $(top_srcdir)/blixemApp
DOTTER = $(top_srcdir)/dotterApp
BELVU = $(top_srcdir)/belvuApp
DOC = $(top_srcdir)/doc
USERDOC = $(top_srcdir)/doc/User_doc


#
# Convenience libraries
#
# These get compiled into the top-level directory to make finding
# make targent dependencies easier, even though the source-code 
# lives in subdirectories
#

noinst_LIBRARIES = libSeqtoolsUtils.a libpfetch.a libBlixemApp.a libDotterApp.a libBelvuApp.a

# Define the sources for the convenience libraries
libSeqtoolsUtils_a_SOURCES = $(UTILS)/iupac.h $(UTILS)/version.h $(UTILS)/utilities.h $(UTILS)/utilities.c $(UTILS)/blxmsp.h $(UTILS)/blxmsp.c $(UTILS)/translate.c $(UTILS)/seqtoolsWebBrowser.c $(UTILS)/blxGff3Parser.h $(UTILS)/blxGff3Parser.c $(UTILS)/blxparser.h $(UTILS)/blxparser.c
libpfetch_a_SOURCES = $(PFETCH)/libpfetch-utils.c $(PFETCH)/libpfetch.c $(PFETCH)/libpfetch.h $(PFETCH)/libpfetch-cmarshal.c $(PFETCH)/libcurlobject.c $(PFETCH)/libcurlobject.h $(PFETCH)/libpfetch_I.h $(PFETCH)/libpfetch_P.h $(PFETCH)/libcurlobject_I.h
libBlixemApp_a_SOURCES = $(BLIXEM)/sequencecellrenderer.c $(BLIXEM)/bigpicture.c $(BLIXEM)/bigpicturegrid.c $(BLIXEM)/detailview.c $(BLIXEM)/detailviewtree.c $(BLIXEM)/blxwindow.c $(BLIXEM)/exonview.c $(BLIXEM)/coverageview.c $(BLIXEM)/blxdotter.c $(BLIXEM)/blxview.h $(BLIXEM)/blixem_.h $(BLIXEM)/detailview.h $(BLIXEM)/detailviewtree.h $(BLIXEM)/sequencecellrenderer.h $(BLIXEM)/bigpicture.h $(BLIXEM)/bigpicturegrid.h $(BLIXEM)/blxwindow.h $(BLIXEM)/exonview.h $(BLIXEM)/coverageview.h $(BLIXEM)/blxdotter.h
libDotterApp_a_SOURCES = $(DOTTER)/greyramptool.c $(DOTTER)/alignmenttool.c $(DOTTER)/dotplot.c $(DOTTER)/dotter.c $(DOTTER)/dotterKarlin.c $(DOTTER)/seqtoolsExonView.c $(DOTTER)/dotter_.h $(DOTTER)/dotter.h $(DOTTER)/seqtoolsExonView.h
libBelvuApp_a_SOURCES = $(BELVU)/belvu_.h $(BELVU)/belvu.c $(BELVU)/belvuWindow.h $(BELVU)/belvuWindow.c $(BELVU)/belvuAlignment.h $(BELVU)/belvuAlignment.c $(BELVU)/belvuTree.h $(BELVU)/belvuTree.c $(BELVU)/belvuConsPlot.h $(BELVU)/belvuConsPlot.c

# libpfetch marshal stuff
libpfetch_a_marshal = $(PFETCH)/libpfetch-cmarshal.list
libpfetch_a_DEPENDENCIES = $(nodist_libpfetch_a_SOURCES)
libpfetch_a_extra_sources = $(PFETCH)/libpfetch-cmarshal.h $(PFETCH)/libpfetch-cmarshal.c $(NULL)
nodist_libpfetch_a_SOURCES = $(libpfetch_a_extra_sources)
BUILT_SOURCES = $(libpfetch_a_extra_sources)
CLEANFILES = $(BUILT_SOURCES)
PFETCH_EXTRA_DIST = $(libpfetch_a_marshal)
GENMARSHAL_COMMAND = $(GLIB_GENMARSHAL) --prefix=libpfetch_cmarshal
.list.h:
	$(GENMARSHAL_COMMAND) --header $< >$@
.list.c:
	$(GENMARSHAL_COMMAND) --body $< >$@


#
# Documentation
#

# these will be installed alongside the programs
dist_doc_DATA = $(USERDOC)/todo.html \
                $(USERDOC)/belvu_quick_start.html $(USERDOC)/blixem_quick_start.html $(USERDOC)/dotter_quick_start.html \
                $(USERDOC)/belvu_usage.txt $(USERDOC)/blixem_usage.txt $(USERDOC)/dotter_usage.txt \
                $(USERDOC)/belvu_news.html $(USERDOC)/blixem_news.html $(USERDOC)/dotter_news.html \
                $(USERDOC)/Blixem_manual.pdf $(USERDOC)/Dotter_manual.pdf

# these are included in the distribution, but are not installed
EXTRA_DIST = $(PFETCH_EXTRA_DIST) $(DOC)/Design_notes $(DOC)/index.html $(USERDOC)/index.html $(USERDOC)/manuals.html Blixem_$(USERDOC)/manual.odt Dotter_$(USERDOC)/manual.odt


#
# Executables
#
bin_PROGRAMS = dotter blixem blixemh belvu

# blixem and blixemh have the same source files, most of which are compiled into libBlixemApp.a. However, there are a 
# few that need to be compiled with different flags/libraries, and those are specified here
blixem_SOURCES = $(top_srcdir)/blixemApp/blxmain.c $(top_srcdir)/blixemApp/blxview.c $(top_srcdir)/blixemApp/blxFetch.c
blixemh_SOURCES = $(blixem_SOURCES)
dotter_SOURCES = $(top_srcdir)/dotterApp/dotterMain.c
belvu_SOURCES = $(top_srcdir)/belvuApp/belvuMain.c

AM_CPPFLAGS = $(DEPS_CFLAGS) -I$(top_srcdir)
blixemh_CPPFLAGS = $(AM_CPPFLAGS) -DPFETCH_HTML

AM_LDFLAGS = -L$(top_builddir)/seqtoolsUtils
blixem_LDFLAGS = $(AM_LDFLAGS) -L$(top_builddir)/blixemApp
blixemh_LDFLAGS = $(blixem_LDFLAGS)
dotter_LDFLAGS = $(AM_LDFLAGS) -L$(top_builddir)/dotterApp
belvu_LDFLAGS = $(AM_LDFLAGS) -L$(top_builddir)/belvuApp

COMMON_LIBS = $(DEPS_LIBS) $(top_builddir)/libSeqtoolsUtils.a

blixem_LDADD = $(COMMON_LIBS) $(top_builddir)/libBlixemApp.a
blixemh_LDADD = $(COMMON_LIBS) $(top_builddir)/libBlixemApp.a $(top_builddir)/libpfetch.a -lcurl
dotter_LDADD = $(COMMON_LIBS) $(top_builddir)/libDotterApp.a
belvu_LDADD = $(COMMON_LIBS) $(top_builddir)/libBelvuApp.a


