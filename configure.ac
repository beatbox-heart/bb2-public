AC_INIT([BeatBox],[2.0],[V.N.Biktashev@exeter.ac.uk])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE([-Wall gnu])
myhost=`uname -n`
case $myhost in
   hector*) 
     LIBS=" -lm -lX11 -lxcb -lxcb-xlib -ldl -lXau";;
   *) 
   ;;
   esac
myos=`uname`
case $myos in
   Darwin*)
     LIBS="-L/opt/X11/lib"
     CFLAGS="$CFLAGS -I/opt/X11/include"
   ;;
   *)
   ;;
   esac

# Checks for graphics libraries
AC_SEARCH_LIBS(XOpenDisplay,
  [X11],
  [AC_SUBST(X11_AVAIL,yes)],
  [AC_MSG_WARN([X11 libraries not found, Beatbox will be built without graphical displays.  If you need graphical displays you will need to configure Beatbox with X11, please contact the Beatbox developers for further information.])
    AC_SUBST(X11_AVAIL,no)],
  [])
if test "x$X11_AVAIL" = "xno"; then
  AC_SUBST(GL_AVAIL,no)
else
  AC_SEARCH_LIBS(glXQueryExtension,
    [GL],
    [AC_SUBST(GL_AVAIL,yes)],
    [AC_MSG_WARN([OpenGL libraries not found, Beatbox will be built without GL graphical extention. If you need graphical displays you will need to configure Beatbox with OpenGL please contact the Beatbox developers for further information.])
      AC_SUBST(GL_AVAIL,no)],
  [])
fi

AC_SEARCH_LIBS(pow, [m], [], [])
AC_PROG_CC(mpicc gcc cl cc)
AM_PROG_CC_C_O
AC_CONFIG_HEADERS([config.h])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include "mpi.h"]], [[MPI_Init(0, 0);]])],[mpi_link="yes"],[mpi_link="no"])
AC_CONFIG_FILES([
	Makefile
	src/Makefile
])

if test "x$X11_AVAIL" = "xno"; then
  echo "***************"
  echo "Reminder: X11 libraries not found so the graphical part of Beatbox is not being built.  If you want graphical display then please contact the Beatbox developers for help building with X11 libraries."
  echo "***************"
  CFLAGS="$CFLAGS -D NOX11"
else
  CFLAGS="$CFLAGS -D X11"
fi

if test "$GL_AVAIL" = "no"; then
  echo "***************"
  echo "Reminder: OpenGL libraries not found, Beatbox will be built without GL graphical extention. If you need graphical displays you will need to configure Beatbox with X11, please contact the Beatbox developers for further information."
  echo "***************"
  CFLAGS="$CFLAGS -D NOGL"
else
  CFLAGS="$CFLAGS -D GL"
fi

LT_INIT

AC_CONFIG_SUBDIRS([src/gsl-1.16.extract])
AC_OUTPUT
