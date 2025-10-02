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
  [AC_MSG_WARN([X11 libraries not found, Beatbox will be built without graphical displays.])
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
# GLX detection (requires X11 + GL)
# -------------------------------
GLX_CFLAGS=""
GLX_LIBS=""

if test "x$X11_AVAIL" = "xno" || test "x$GL_AVAIL" = "xno"; then
  AC_MSG_WARN([Skipping GLX check because X11 or GL not available])
  AC_SUBST(GLX_AVAIL,no)
else
  AC_CHECK_HEADER([GL/glx.h],
    [AC_SEARCH_LIBS([glXQueryExtension],
       [GL GLX],
       [AC_SUBST(GLX_AVAIL,yes)
        GLX_LIBS="-lGL"
        AC_DEFINE([HAVE_GLX], [1], [Define if GLX is available])],
       [AC_MSG_WARN([GLX libraries not found, building without GLX support])
        AC_SUBST(GLX_AVAIL,no)])],
    [AC_MSG_WARN([GLX header (GL/glx.h) not found, building without GLX support])
     AC_SUBST(GLX_AVAIL,no)])
fi

AC_SUBST(GLX_CFLAGS)
AC_SUBST(GLX_LIBS)

# -------------------------------
# CFLAGS adjustments for X11/GL/GLX
# -------------------------------
if test "x$X11_AVAIL" = "xno"; then
  echo "***************"
  echo "Reminder: X11 libraries not found so the graphical part of Beatbox is not being built."
  echo "***************"
  CFLAGS="$CFLAGS -D NOX11"
else
  CFLAGS="$CFLAGS -D X11"
fi

if test "$GL_AVAIL" = "no"; then
  echo "***************"
  echo "Reminder: OpenGL libraries not found, Beatbox will be built without GL graphical extension."
  echo "***************"
  CFLAGS="$CFLAGS -D NOGL"
else
  CFLAGS="$CFLAGS -D GL"
  LIBS="$LIBS $GL_LIBS"
fi

if test "$GLX_AVAIL" = "no"; then
  echo "***************"
  echo "Reminder: GLX not found, Beatbox will be built without GLX support."
  echo "***************"
  CFLAGS="$CFLAGS -D NOGLX"
else
  CFLAGS="$CFLAGS -D GLX"
fi

# -------------------------------
# Other checks
# -------------------------------
AC_SEARCH_LIBS(pow, [m], [], [])

# -------------------------------
# Compiler detection
# -------------------------------
AC_PROG_CC
AM_PROG_CC_C_O

# -------------------------------
# MPI detection (override CC if mpicc found)
# -------------------------------
AC_PATH_PROG([MPICC], [mpicc], [no])

if test "x$MPICC" != "xno"; then
   AC_MSG_RESULT([MPI compiler found: $MPICC])
   CC="$MPICC"
   AC_DEFINE([HAVE_MPI], [1], [Define if MPI is available])
   AM_CONDITIONAL([MPI_AVAILABLE], [true])
else
   AC_MSG_WARN([MPI compiler not found, building without MPI support])
   AM_CONDITIONAL([MPI_AVAILABLE], [false])
fi

# -------------------------------
# GLX availability summary (double-check)
# -------------------------------
if test "x$X11_AVAIL" = "xyes" -a "x$GL_AVAIL" = "xyes"; then
  AC_SEARCH_LIBS([glXQueryExtension], [GL], 
    [GLX_AVAIL=yes], 
    [GLX_AVAIL=no])
else
  GLX_AVAIL=no
fi
AC_SUBST(GLX_AVAIL)

# -------------------------------
# Summary
# -------------------------------
echo ""
echo "================= Configuration Summary ================="
echo " C compiler:       $CC"
if test "x$MPICC" != "xno"; then
  echo " MPI support:      enabled ($MPICC)"
else
  echo " MPI support:      disabled"
fi
echo " X11 support:      $X11_AVAIL"
echo " OpenGL support:   $GL_AVAIL"
echo " GLX support:      $GLX_AVAIL"
echo "========================================================="

# -------------------------------
# Headers
# -------------------------------
AC_CONFIG_HEADERS([config.h])

# -------------------------------
# Libtool
# -------------------------------
LT_INIT

# -------------------------------
# Subdirectories and output
# -------------------------------
AC_CONFIG_SUBDIRS([src/gsl-1.16.extract], [CC="$CC"])
AC_CONFIG_FILES([
	Makefile
	src/Makefile
])

AC_OUTPUT
