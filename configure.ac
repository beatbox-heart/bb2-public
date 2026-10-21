AC_INIT([BeatBox],[2.2],[V.N.Biktashev@exeter.ac.uk])
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE([-Wall gnu])

# Generated files (configure, Makefile.in) are distributed, and their
# timestamps are not preserved by git, so by default do not try to
# re-run autotools from make.  Developers who edit configure.ac or
# Makefile.am should either run "autoreconf -fi" by hand, or configure
# with --enable-maintainer-mode to have make do it automatically.
AM_MAINTAINER_MODE([disable])

# -------------------------------
# Compiler
# -------------------------------
# Prefer the MPI wrapper; a user-supplied CC=... takes precedence.
# This must come before anything modifies CFLAGS, otherwise the
# default "-g -O2" is lost.
AC_PROG_CC([mpicc gcc cc])
AM_PROG_CC_C_O
LT_INIT

# -------------------------------
# Platform-specific paths
# -------------------------------
case `uname` in
   Darwin*)
     # XQuartz
     CPPFLAGS="$CPPFLAGS -I/opt/X11/include"
     LDFLAGS="$LDFLAGS -L/opt/X11/lib"
   ;;
   *) ;;
esac

AC_SEARCH_LIBS(pow, [m], [], [])

# -------------------------------
# MPI: build bb2p only if the compiler can actually link an MPI program
# -------------------------------
AC_MSG_CHECKING([whether $CC can build MPI programs])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]], [[MPI_Init(0, 0);]])],
  [mpi_ok=yes], [mpi_ok=no])
AC_MSG_RESULT([$mpi_ok])
if test "x$mpi_ok" = "xno"; then
  AC_MSG_WARN([MPI not available, the parallel version bb2p will not be built. To build it, install MPI and make sure mpicc is in PATH, or set CC to an MPI compiler wrapper.])
fi
AM_CONDITIONAL([MPI_AVAILABLE], [test "x$mpi_ok" = "xyes"])

# -------------------------------
# Graphics
# -------------------------------
AC_ARG_WITH([mesa],
  AS_HELP_STRING([--with-mesa], [find OpenGL via pkg-config (e.g. Homebrew Mesa) instead of searching for libGL]),
  [use_mesa=$withval], [use_mesa=no])

AC_SEARCH_LIBS(XOpenDisplay, [X11],
  [X11_AVAIL=yes],
  [AC_MSG_WARN([X11 libraries not found, Beatbox will be built without graphical displays.])
   X11_AVAIL=no])

GL_AVAIL=no
if test "x$X11_AVAIL" = "xyes"; then
  AS_IF([test "x$use_mesa" = "xyes"],
    # pkg-config is called directly: PKG_CHECK_MODULES would make
    # pkg-config a requirement even when --with-mesa is not given.
    [AC_PATH_TOOL([PKG_CONFIG], [pkg-config])
     AS_IF([test -n "$PKG_CONFIG" && $PKG_CONFIG --exists gl],
       [GL_AVAIL=yes
        CPPFLAGS="$CPPFLAGS `$PKG_CONFIG --cflags gl`"
        LIBS="`$PKG_CONFIG --libs gl` $LIBS"],
       [AC_MSG_WARN([--with-mesa given but pkg-config cannot find gl (on macOS: brew install mesa)])])],
    [AC_SEARCH_LIBS(glXQueryExtension, [GL], [GL_AVAIL=yes])])
  # The GL sources include <GL/gl.h> and <GL/glx.h>
  if test "x$GL_AVAIL" = "xyes"; then
    AC_CHECK_HEADER([GL/glx.h], [], [GL_AVAIL=no])
  fi
  if test "x$GL_AVAIL" = "xno"; then
    AC_MSG_WARN([OpenGL libraries not found, Beatbox will be built without GL graphical extension.])
  fi
fi

# Passed to the compiler via src/Makefile.am rather than appended to
# CFLAGS, so that they survive "make CFLAGS=..." and do not leak into
# the GSL sub-configure.
GRAPHICS_CPPFLAGS=""
if test "x$X11_AVAIL" = "xyes"; then
  GRAPHICS_CPPFLAGS="$GRAPHICS_CPPFLAGS -D X11"
else
  GRAPHICS_CPPFLAGS="$GRAPHICS_CPPFLAGS -D NOX11"
fi
if test "x$GL_AVAIL" = "xyes"; then
  GRAPHICS_CPPFLAGS="$GRAPHICS_CPPFLAGS -D GL"
else
  GRAPHICS_CPPFLAGS="$GRAPHICS_CPPFLAGS -D NOGL"
fi
AC_SUBST([GRAPHICS_CPPFLAGS])

# -------------------------------
# Output
# -------------------------------
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
	Makefile
	src/Makefile
])
AC_CONFIG_SUBDIRS([src/gsl-1.16.extract])
AC_OUTPUT

# Printed last, so it is not buried under the GSL sub-configure output
echo ""
echo "================= Configuration Summary ================="
echo " C compiler:       $CC"
echo " MPI (bb2p):       $mpi_ok"
echo " X11 support:      $X11_AVAIL"
echo " OpenGL support:   $GL_AVAIL"
echo "========================================================="
