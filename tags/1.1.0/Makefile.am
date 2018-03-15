ACLOCAL_AMFLAGS=-I ./make/m4

SUBDIRS=src/stat_calc src/data_manage src/annie src/galib247 src/libGE-0.26 src/athena

EXTRA_DIST=src/libGE-0.26/GE doc CHANGELOG sample_sets samples scripts README

CLEANFILES=src/libGE-0.26/GE doc CHANGELOG example.gram sample_sets makeathenaconfig.pl

# This is here to prevent installation of LOKI during a distcheck
DISTCHECK_CONFIGURE_FLAGS=--disable-loki

