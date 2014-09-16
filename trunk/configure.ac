# Ritchie Lab configure.ac template

AC_PREREQ(2.59)
AC_INIT(athena, 1.1.0, software@ritchielab.psu.edu)

# Checks for --enable-debug added
AX_CHECK_ENABLE_DEBUG

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX

RL_CONFIGURE()

AM_INIT_AUTOMAKE([1.8 foreign])

AC_CONFIG_SRCDIR([src/athena/athena.cpp])
AC_CONFIG_HEADER([config.h])


AC_CONFIG_FILES([Makefile
                 src/data_manage/Makefile
                 src/stat_calc/Makefile
                 src/annie/Makefile
                 src/galib247/Makefile
                 src/libGE-0.26/Makefile
                 src/athena/Makefile
                ])
AC_OUTPUT
