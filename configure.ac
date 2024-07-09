AC_INIT([abn],[3.1.1])

## Use gsl-config to find arguments for compiler and linker flags
##
## Check for non-standard programs: gsl-config(1)
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
## If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then
    # Use gsl-config for header and linker arguments
    GSL_CFLAGS=`${GSL_CONFIG} --cflags`
    GSL_LIBS=`${GSL_CONFIG} --libs`
else
    AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi

# Now substitute these variables in src/Makevars.in to create src/Makevars
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

## Check for non-standard programs: JAGS
AC_PATH_PROG([JAGS], [jags])
if test "${JAGS}" != ""; then
		## Try to fetch the paths for jags with pkg-config(1)
		if test "${PKG_CONFIG} --modversion ${JAGS}" != ""; then
          JAGS_CFLAGS=`pkg-config --cflags --static jags`
          JAGS_LIBS=`pkg-config --libs --static jags`
		fi
else
    AC_MSG_ERROR([jags not found, is JAGS installed?])
fi
#
# Now substitute these variables in src/Makevars.in to create src/Makevars
AC_SUBST(JAGS_CFLAGS)
AC_SUBST(JAGS_LIBS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
