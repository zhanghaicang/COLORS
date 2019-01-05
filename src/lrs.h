#ifndef COLORS_LRS_H_
#define COLORS_LRS_H_

#include <vector>

#include <Eigen/Core>
#include <Eigen/SVD>
#include "type.h"

namespace colors {

inline int count_larger_than(const Vector& v, double value);

int calc_lrs(const double lambda, const Matrix& M, Matrix& L, Matrix& S);

int calc_lrs(const double lambda, const Matrix& M, const char* output_path);

int save_lrs(const Matrix& M, const Matrix& S, const Matrix& apc,
             const char* output_path);

}  // namespace colors
#endif
