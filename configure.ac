AC_INIT([MJMrss], 1.0.3)
AC_LANG(C++)
: ${CXXFLAGS=""}
# Check openmp
AC_OPENMP
# Check optimiztion flags
AX_CHECK_COMPILE_FLAG([-O2], [CXXFLAGS="$CXXFLAGS -O2"])
m4_define([LTO_CXXFLAGS], "")
AX_CHECK_COMPILE_FLAG([-flto], [[LTO_CXXFLAGS]="-flto"])
AX_CHECK_COMPILE_FLAG([-flto=thin], [[LTO_CXXFLAGS]="-flto=thin"])
AX_CHECK_COMPILE_FLAG([-ffat-lto-objects], [[LTO_CXXFLAGS]+=" -ffat-lto-objects"])
AX_CHECK_COMPILE_FLAG([-ffp-contract=fast], [CXXFLAGS="$CXXFLAGS -ffp-contract=fast"])
# Check system Eigen3
AC_CHECK_HEADERS([eigen3/Eigen/Dense], [have_eigen_in_system=yes], [have_eigen_in_system=no])
AS_IF([have_eigen_in_system=yes], 
[
AC_RUN_IFELSE([AC_LANG_PROGRAM([
#include <eigen3/Eigen/Eigen>
],[
#if EIGEN_VERSION_AT_LEAST(3,3,10)
exit(0);
#else
exit(1);
#endif
])], [have_eigen_in_system=1], [have_eigen_in_system=0])])
# Ouput files
AC_CONFIG_FILES([src/Makevars])
AC_CONFIG_FILES([src/eigen_wrap.h])
#AC_CONFIG_HEADERS([src/eigen_wrap.h])
AC_OUTPUT
AC_SUBST([LTO_CXXFLAGS])
AC_SUBST([have_eigen_in_system])