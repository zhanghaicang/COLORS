#ifndef COLORS_SCORE_H_
#define COLORS_SCORE_H_

#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SVD>

#include "type.h"

namespace colors {

const int ALPHA = 21;

int calc_matrix(const std::vector<std::string>& msa, const Vector& weight,
                Matrix& mi, Matrix& omes, Matrix& cov);

double calc_mi(const Matrix& pab, const Vector& pa, const Vector& pb);
double calc_cov(const Matrix& pab, const Vector& pa, const Vector& pb);
double calc_omes(const Matrix& pab, const Vector& pa, const Vector& pb);

int calc_apc(const Matrix& m, Matrix& apc);

}  // namespace colors

#endif
