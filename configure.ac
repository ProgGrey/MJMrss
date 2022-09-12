# Checking:
AC_INIT([MJMrss], 1.0.1)
AC_LANG(C++)
: ${CXXFLAGS=""}
AC_OPENMP
AX_CHECK_COMPILE_FLAG([-O2], [CXXFLAGS="$CXXFLAGS -O2"])
m4_define([LTO_CXXFLAGS], "")
AX_CHECK_COMPILE_FLAG([-flto], [[LTO_CXXFLAGS]="-flto"])
AX_CHECK_COMPILE_FLAG([-flto=thin], [[LTO_CXXFLAGS]="-flto=thin"])
AX_CHECK_COMPILE_FLAG([-ffat-lto-objects], [[LTO_CXXFLAGS]+=" -ffat-lto-objects"])
AX_CHECK_COMPILE_FLAG([-ffp-contract=fast], [CXXFLAGS="$CXXFLAGS -ffp-contract=fast"])
AX_CXX_COMPILE_STDCXX([14], [ext], [mandatory])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
AC_SUBST([LTO_CXXFLAGS])