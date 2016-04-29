# Ritchie Lab configure.ac template

AC_PREREQ(2.59)
AC_INIT(athena, 1.1.0, software@ritchielab.psu.edu)

# Checks for --enable-debug added
AX_CHECK_ENABLE_DEBUG
#AC_PROG_CC
#AC_PROG_CPP
#AC_PROG_CXX

AC_CANONICAL_HOST
case $host_os in
  darwin*)
    CPPFLAGS="$CPPFLAGS --stdlib=libstdc++"
    LDFLAGS="$LDFLAGS --stdlib=libstdc++"
  ;;
  *)
  ;;
esac

AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
[compile with MPI (parallelization) support. If none is found, MPI is not used. Default: auto]) ],
,[with_mpi=auto])

AX_PROG_CXX_MPI([test x"$with_mpi" != xno],[AC_DEFINE([HAVE_CXX_MPI],[],[CXX MPI Found])],[
   AC_DEFINE([NO_CXX_MPI],[],[No C MPI Found])
   if test x"$with_mpi" = xyes; then
     AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
   else
     AC_MSG_WARN([No MPI compiler found, won't use MPI.])
   fi
])

AX_PROG_CC_MPI([test x"$with_mpi" != xno],[AC_DEFINE([HAVE_CC_MPI],[],[C MPI Found])],[
   AC_DEFINE([NO_CC_MPI],[],[No C MPI Found])
   if test x"$with_mpi" = xyes; then
     AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
   else
     AC_MSG_WARN([No MPI compiler found, won't use MPI.])
   fi
])

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
