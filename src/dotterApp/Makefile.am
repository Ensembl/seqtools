
SUBDIRS = .

include $(top_srcdir)/Makefile.am.common

# there's a problem on ubuntu 14 with finding X11 so this is a hopefully temporary hack...
#X_LIB = -lX11 -lm

bin_PROGRAMS = dotter
dotter_SOURCES = dotterMain.cpp greyramptool.cpp alignmenttool.cpp dotplot.cpp dotter.cpp dotterKarlin.cpp seqtoolsExonView.cpp dotter_.hpp dotter.hpp seqtoolsExonView.hpp 
dotter_LDADD = $(top_builddir)/seqtoolsUtils/libSeqtoolsUtils.a 

# If gbtools is in a subdirectory, add it; otherwise look for a local installation
if USE_GBTOOLS
dotter_LDADD += $(top_builddir)/gbtools/.libs/libgbtools.a            
else
dotter_LDADD += -lgbtools
endif

# the gtk deps etc. must go at the end so that gbtools can pick them up
dotter_LDADD += $(DEPS_LIBS) $(X_LIB)

# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/dotterApp/Makefile.in
