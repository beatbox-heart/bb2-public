AC_INIT([BeatBox],[2.2],[V.N.Biktashev@exeter.ac.uk])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE([-Wall gnu])

# Initialize pkg-config (needed for PKG_CHECK_MODULES macros)
m4_pattern_forbid([^PKG_CHECK_MODULES$])
m4_include([m4/pkg.m4])

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

# -------------------------------
# Option: allow forcing Mesa/Zink
# -------------------------------
AC_ARG_WITH([mesa],
  AS_HELP_STRING([--with-mesa], [Use Mesa pkg-config instead of system OpenGL]),
  [use_mesa=$withval], [use_mesa=no])

# -------------------------------
# Checks for graphics libraries
# -------------------------------
AC_SEARCH_LIBS(XOpenDisplay,
  [X11],
  [AC_SUBST(X11_AVAIL,yes)],
  [AC_MSG_WARN([X11 libraries not found, Beatbox will be built without graphical displays.  If you need graphical displays you will need to configure Beatbox with X11, please contact the Beatbox developers for further information.])
    AC_SUBST(X11_AVAIL,no)],
  [])

# -------------------------------
# OpenGL detection (cross-platform)
# -------------------------------
GL_CFLAGS=""
GL_LIBS=""

if test "x$X11_AVAIL" = "xno"; then
  AC_SUBST(GL_AVAIL,no)
else
  case "$myos" in
    Darwin*)
      if test "x$use_mesa" = "xyes"; then
        PKG_CHECK_MODULES([GL], [gl], [], 
          [AC_MSG_WARN([Mesa not found, install with `brew install mesa`])])
        GL_CFLAGS="$GL_CFLAGS $GL_CFLAGS"
        GL_LIBS="$GL_LIBS $GL_LIBS"
        AC_SUBST(GL_AVAIL,yes)
      else
        AC_CHECK_HEADER([OpenGL/gl.h], [],
          [AC_MSG_WARN([MacOS OpenGL headers not found (install Xcode Command Line Tools)])])
	GL_LIBS="-L/opt/X11/lib -lGL -lGLU -lX11"
        AC_SUBST(GL_AVAIL,yes)

        AS_WARN([
**********************************************************************
* You are building on macOS using Apple's OpenGL framework.          *
* This implementation is deprecated and limited to OpenGL 4.1 Core.  *
* For newer OpenGL (4.6+), consider --with-mesa and Homebrew Mesa.   *
**********************************************************************
])
      fi
      ;;
    *)
      # Linux and others
      AC_SEARCH_LIBS(glXQueryExtension,
        [GL],
        [AC_SUBST(GL_AVAIL,yes) GL_LIBS="-lGL -lGLU"],
        [AC_MSG_WARN([OpenGL libraries not found, Beatbox will be built without GL graphical extension.])
          AC_SUBST(GL_AVAIL,no)],
      [])
      ;;
  esac
fi

AC_SUBST(GL_CFLAGS)
AC_SUBST(GL_LIBS)

# -------------------------------
# CFLAGS adjustments for X11/GL
# -------------------------------
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
  LIBS="$LIBS $GL_LIBS"
fi

# -------------------------------
# Other checks
# -------------------------------
AC_SEARCH_LIBS(pow, [m], [], [])
AC_PROG_CC(mpicc gcc cl cc)
AM_PROG_CC_C_O
AC_CONFIG_HEADERS([config.h])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include "mpi.h"]], [[MPI_Init(0, 0);]])],[mpi_link="yes"],[mpi_link="no"])

# -------------------------------
# Libtool
# -------------------------------
LT_INIT

# -------------------------------
# Subdirectories and output
# -------------------------------
AC_CONFIG_SUBDIRS([src/gsl-1.16.extract])
AC_CONFIG_FILES([
	Makefile
	src/Makefile
])

AC_OUTPUT
