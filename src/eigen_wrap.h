#define HAVE_EIGEN_IN_SYSTEM 0

#if HAVE_EIGEN_IN_SYSTEM == 1
// Use system Eigen (3.4.x) instead RcppEigen
#include <eigen3/Eigen/Dense>

#else
// Use RcppEigen (3.3.9)
#include <Eigen/Dense>

#if (EIGEN_WORLD_VERSION == 3) && (EIGEN_MAJOR_VERSION == 3)
namespace Eigen{
    template <typename Type>
    using VectorX = Matrix<Type, Dynamic, 1>;
};
#endif

#endif
