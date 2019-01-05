#ifndef COLORS_TYPE_H_
#define COLORS_TYPE_H_

#include <Eigen/Core>
#include <Eigen/SVD>

namespace colors {

typedef double FloatType;

typedef Eigen::Matrix<FloatType, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Array<FloatType, Eigen::Dynamic, Eigen::Dynamic> Array;
typedef Eigen::Matrix<FloatType, Eigen::Dynamic, 1> Vector;

}  // namespace colors
#endif
